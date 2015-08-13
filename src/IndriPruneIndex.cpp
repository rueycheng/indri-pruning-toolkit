#include "IndriPruneIndex.h"

#include <algorithm>
#include <string>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <boost/algorithm/string.hpp>
#include <boost/bind.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/function.hpp>

// FIXME: check if this quick fix is appropriate
#define __atomic_add __gnu_cxx::__atomic_add

#include "indri/CompressedCollection.hpp"
#include "indri/DeletedDocumentList.hpp"
#include "indri/DiskIndex.hpp"
#include "indri/IndexWriter.hpp"
#include "indri/Parameters.hpp"
#include "indri/Path.hpp"
#include "indri/Repository.hpp"
#include "indri/TermScoreFunctionFactory.hpp"

#include "BinnedFilter.h"
#include "DocumentCentricTermScoreFunction.h"
#include "FilteredIndex.h"
#include "GainFunction.h"
#include "QuantileSampler.h"
#include "RunningStat.h"

indri_contrib::DocumentCentricTermScoreFunction*
getDocumentCentricTermScoreFunction(const std::string& stringSpec, 
				    double documentLength, double collectionLength, int documentCount)
{
    // Notation:
    //   collectionLength    := N  (collection length, in # of words)
    //   documentCount       := D  (total # of documents)
    //   documentLength      := DL (document length, in # of words)
    //
    // When invoked, pass 3 parameters:
    //   termCount           := TF (term frequency in some document)
    //   collectionTermCount := CF (collection term frequency)
    //   documentFrequency   := DF (document frequency, in # of docs)
    indri::api::Parameters spec;
    parseSpec(spec, stringSpec);
    std::string method = spec.get("method", "dirichlet");

    if (method == "dirichlet" || method == "d" || method == "dir") {
	// dirichlet -- takes parameter "mu"
	double mu = spec.get("mu", 2500);
	return new indri_contrib::DirichletDocumentCentricTermScoreFunction(
	    mu, documentLength, collectionLength);
    } 
    else if (method == "linear" || method == "jm" || method == "jelinek-mercer") {
	// jelinek-mercer -- can take parameters collectionLambda (or just lambda) and documentLambda
	double collectionLambda = spec.exists("collectionLambda")?
	    spec.get("collectionLambda", 0.4): spec.get("lambda", 0.4);
	return new indri_contrib::JelinekMercerDocumentCentricTermScoreFunction(
	    collectionLambda, documentLength, collectionLength);
    }
    else if (method == "tfidf" || method == "okapi") {
	double k1 = spec.get("k1", 1.2);
	double b = spec.get("b", 0.75);
	double k3 = spec.get("k3", 7);
	int qtf = spec.get("qtf", 1);

	return new indri_contrib::TFIDFDocumentCentricTermScoreFunction(
	    qtf, k1, b, k3, (method == "okapi"), documentLength, collectionLength, documentCount);
    } 
    else if (method == "kld") {
	double collectionLambda = spec.exists("collectionLambda")?
	    spec.get("collectionLambda", 0.0): spec.get("lambda", 0.0);
	return new indri_contrib::KLDivergenceDocumentCentricTermScoreFunction(
	    collectionLambda, documentLength, collectionLength);
    }
    else die(boost::format("No support for method '%s'") % method);
}

// Weight
//
typedef std::vector<double> Weight;

// LazyVectorLoader
//
class LazyVectorLoader {
    std::ifstream _in;
    std::string _buf;
    std::istringstream _linein;

public:
    LazyVectorLoader(std::string guardPath): _in(guardPath) { }

    bool hasNext() const { return _in; }

    void next(std::unordered_set<int>& v) {
	v.clear();

	if (getline(_in, _buf)) {
	    _linein.str(_buf);
	    _linein.clear();

	    int n;
	    while (_linein >> n) v.insert(n);
	}
    }
};


// Generalized uniform pruning
//
template<class GainFunction>
void generalized_uniform_pruning(indri::index::Index& index, 
				 indri::api::Parameters& parameters, 
				 Weight& weight,
				 indri_contrib::QuantileSampler* sampler,
				 indri_contrib::BinnedFilter* filter)
{
    string rule = parameters.get("rule", "method:dir");
    double epsilon = parameters.get("epsilon", 0.0);
    double scale = parameters.get("scale", 1.0);
    bool softmax = parameters.get("softmax", 0);

    string partial_gain = parameters.get("partialGain", "");
    bool point_eval_exclusive = (partial_gain == "exclusive");
    bool point_eval_inclusive = (partial_gain == "inclusive");

    string rank_transform = parameters.get("rankTransform", "");
    bool step_transform = (rank_transform == "step");
    bool uniform_transform = (rank_transform == "uniform");

    unique_ptr<LazyVectorLoader> guardLoader;
    if (parameters.exists("guard")) {
	string guardPath = parameters["guard"];
	guardLoader = unique_ptr<LazyVectorLoader>(new LazyVectorLoader(guardPath));
    }

    // Initialize the gain function
    GainFunction gain;
    gain.loadParameters(parameters);

    // There are totally N term occurrences in D documents
    long N = index.termCount();
    int D = index.documentCount();

    vector<int> CF_cache;
    CF_cache.push_back(0);
    vector<int> DF_cache;
    DF_cache.push_back(0);

    int T = index.uniqueTermCount();
    vector<int> term2seq(T + 1); // map termID to seq # in inverted list
    term2seq[0] = 0;
    // unordered_map<string, int> pos;

    // Pass 1: Term-based scan, for setting up (empty) bitmaps
    {
	unique_ptr<indri::index::DocListFileIterator>
	    docListFileIterator(index.docListFileIterator());
	int seqID = 0;

	for (docListFileIterator->startIteration();
	     !docListFileIterator->finished();
	     docListFileIterator->nextEntry())
	{
	    indri::index::DocListFileIterator::DocListData* docListData = docListFileIterator->currentEntry();
	    assert(docListData);

	    indri::index::TermData* termData = docListData->termData;
	    assert(termData);

	    indri::index::DocListIterator* docListIterator = docListData->iterator;
	    assert(docListIterator);

	    ++seqID;
	    term2seq[index.term(termData->term)] = seqID;

	    int CF = termData->corpus.totalCount; // collection term frequency
	    CF_cache.push_back(CF);

	    int DF = termData->corpus.documentCount; // document frequency
	    DF_cache.push_back(DF);

	    while (!docListIterator->finished())
		docListIterator->nextEntry(D + 1);

	    if (filter) {
		boost::dynamic_bitset<unsigned long> bitmap; 
		filter->add(bitmap);
	    }
	}
    }

    // Pass 2: Document-based scan, to compute scores
    {
	unique_ptr<indri::index::TermListFileIterator>
	    termListFileIterator(index.termListFileIterator());

	unordered_map<lemur::api::TERMID_T, int> bag;
	vector<pair<lemur::api::TERMID_T, double>> prob;
	lemur::api::DOCID_T docID = 0;

	unordered_set<int> guards;
	
	for (termListFileIterator->startIteration();
	     !termListFileIterator->finished();
	     termListFileIterator->nextEntry())
	{
	    ++docID;

	    indri::index::TermList* termList = termListFileIterator->currentEntry();
	    assert(termList);
	
	    bag.clear();
	    const indri::utility::greedy_vector<lemur::api::TERMID_T>& vec = termList->terms();
	    for (size_t i = 0; i < vec.size(); ++i) ++bag[vec[i]];

	    int DL = vec.size();
	    unique_ptr<indri_contrib::DocumentCentricTermScoreFunction> scoreFunction(
		getDocumentCentricTermScoreFunction(rule, DL, N, D));

	    prob.clear();
	    for (unordered_map<lemur::api::TERMID_T, int>::const_iterator iter = 
		 bag.begin(); iter != bag.end(); ++iter)
	    {
		lemur::api::TERMID_T seqID = term2seq[iter->first];
		if (seqID == 0) continue; // FIXME: a quick hack for removing OOVs

		int TF = iter->second;
		double probScore = scoreFunction->compute(TF, CF_cache[seqID], DF_cache[seqID]);
		prob.push_back(make_pair(seqID, probScore));
	    }

	    if (!scoreFunction->isProbability()) {
		double sum = 0.0;
		for (size_t i = 0; i < prob.size(); ++i) {
		    if (prob[i].second < 0) prob[i].second = 0;
		    sum += prob[i].second;
		}

		for (size_t i = 0; i < prob.size(); ++i)
		    prob[i].second = scale * prob[i].second / sum;

		if (softmax) {
		    sum = 0.0;
		    for (size_t i = 0; i < prob.size(); ++i) {
			prob[i].second = std::exp(prob[i].second);
			sum += prob[i].second;
		    }

		    for (size_t i = 0; i < prob.size(); ++i) 
			prob[i].second /= sum;
		}
	    }

	    if (guardLoader)
		guardLoader->next(guards);

	    if (!prob.empty()) {
		double cumprob = 0.0;
		int offset = 0;

		if (!guards.empty()) {
		    for (size_t i = 0; i < prob.size(); ++i) {
			if (guards.find(prob[i].first) != guards.end()) {
			    swap(prob[offset], prob[i]);
			    cumprob += prob[offset++].second;
			}
		    }
		}

		stable_sort(prob.begin() + offset, prob.end(), 
			    greater_value<lemur::api::TERMID_T, double>());

		if (step_transform) {
		    double step = 1.0 / prob.size();
		    double sum = 0.5 * (prob.size() + 1);
		    for (size_t i = 0; i < prob.size(); ++i)
			prob[i].second = (1.0 - i * step) / sum;
		    cumprob = 0.0;
		    offset = 0;
		} 
		else if (uniform_transform) {
		    double step = 1.0 / prob.size();
		    for (size_t i = 0; i < prob.size(); ++i)
			prob[i].second = step;
		    cumprob = 0.0;
		    offset = 0;
		}

		double docWeight = weight[docID];
		double cumgain0 = docWeight * gain(0.0);
		double cumgain = docWeight * gain(cumprob); // save for later

		for (size_t i = offset; i < prob.size(); ++i) {
		    cumprob += prob[i].second;
		    prob[i].second = docWeight * gain(cumprob);
		}

		if (point_eval_inclusive) {
		    double shift = gain(1.0);
		    for (size_t i = 0; i < prob.size(); ++i)
			prob[i].second = shift - prob[i].second;
		}
		else if (point_eval_exclusive) {
		    double shift = gain(1.0);
		    for (size_t i = prob.size() - 1; i > offset; --i)
			prob[i].second = shift - prob[i - 1].second;
		    prob[offset].second = shift - cumgain;
		}
		else {
		    for (size_t i = prob.size() - 1; i > offset; --i)
			prob[i].second -= prob[i - 1].second;
		    prob[offset].second -= cumgain; 
		}

		if (sampler)
		    for (size_t i = 0; i < prob.size(); ++i) 
			sampler->add(prob[i].second);

		if (filter) {
		    for (size_t i = 0; i < offset; ++i) {
			lemur::api::TERMID_T seqID = prob[i].first;
			filter->operator[](seqID - 1).push_back(false);
		    }

		    for (size_t i = offset; i < prob.size(); ++i) {
			lemur::api::TERMID_T seqID = prob[i].first;
			double weighted_gain = prob[i].second;

			if (seqID)
			    filter->operator[](seqID - 1).push_back(weighted_gain < epsilon);
		    }
		}
	    }
	}
    }
}

// PowerFunction
//
struct NegativeExponentPower: std::unary_function<double, double> {
    int minusN;
    NegativeExponentPower(int n): minusN(-n) {}
    double operator()(double p) { return std::pow(p, minusN); }
};

struct NegativeOnePower: std::unary_function<double, double> {
    double operator()(double p) { return 1.0 / p; }
};

// Renyi divergence of infinite order
//
void renyi_divergence_infinite_order(indri::index::Index& index, 
				     indri::api::Parameters& parameters, 
				     Weight& weight,
				     indri_contrib::QuantileSampler* sampler,
				     indri_contrib::BinnedFilter* filter)
{
    string rule = parameters.get("rule", "method:dir");
    double epsilon = parameters.get("epsilon", 0.0);
    double scale = parameters.get("scale", 1.0);
    bool softmax = parameters.get("softmax", 0);
    int cardinality = parameters.get("cardinality", 1);

    boost::function<double(double)> power_function;
    if (cardinality == 1)
	power_function = NegativeOnePower();
    else
	power_function = NegativeExponentPower(cardinality);

    // string partial_gain = parameters.get("partialGain", "");
    // bool point_eval_exclusive = (partial_gain == "exclusive");
    // bool point_eval_inclusive = (partial_gain == "inclusive");

    string rank_transform = parameters.get("rankTransform", "");
    bool step_transform = (rank_transform == "step");
    bool uniform_transform = (rank_transform == "uniform");

    unique_ptr<LazyVectorLoader> guardLoader;
    if (parameters.exists("guard")) {
	string guardPath = parameters["guard"];
	guardLoader = unique_ptr<LazyVectorLoader>(new LazyVectorLoader(guardPath));
    }

    // There are totally N term occurrences in D documents
    long N = index.termCount();
    int D = index.documentCount();

    vector<int> CF_cache;
    CF_cache.push_back(0);
    vector<int> DF_cache;
    DF_cache.push_back(0);

    int T = index.uniqueTermCount();
    vector<int> term2seq(T + 1); // map termID to seq # in inverted list
    term2seq[0] = 0;
    // unordered_map<string, int> pos;

    // Pass 1: Term-based scan, for setting up (empty) bitmaps
    {
	unique_ptr<indri::index::DocListFileIterator>
	    docListFileIterator(index.docListFileIterator());
	int seqID = 0;

	for (docListFileIterator->startIteration();
	     !docListFileIterator->finished();
	     docListFileIterator->nextEntry())
	{
	    indri::index::DocListFileIterator::DocListData* docListData = docListFileIterator->currentEntry();
	    assert(docListData);

	    indri::index::TermData* termData = docListData->termData;
	    assert(termData);

	    indri::index::DocListIterator* docListIterator = docListData->iterator;
	    assert(docListIterator);

	    ++seqID;
	    term2seq[index.term(termData->term)] = seqID;

	    int CF = termData->corpus.totalCount; // collection term frequency
	    CF_cache.push_back(CF);

	    int DF = termData->corpus.documentCount; // document frequency
	    DF_cache.push_back(DF);

	    while (!docListIterator->finished())
		docListIterator->nextEntry(D + 1);

	    if (filter) {
		boost::dynamic_bitset<unsigned long> bitmap; 
		filter->add(bitmap);
	    }
	}
    }

    // Pass 2: Document-based scan, to compute scores
    {
	unique_ptr<indri::index::TermListFileIterator>
	    termListFileIterator(index.termListFileIterator());

	unordered_map<lemur::api::TERMID_T, int> bag;
	vector<pair<lemur::api::TERMID_T, double>> prob;
	lemur::api::DOCID_T docID = 0;

	unordered_set<int> guards;
	
	for (termListFileIterator->startIteration();
	     !termListFileIterator->finished();
	     termListFileIterator->nextEntry())
	{
	    ++docID;

	    indri::index::TermList* termList = termListFileIterator->currentEntry();
	    assert(termList);
	
	    bag.clear();
	    const indri::utility::greedy_vector<lemur::api::TERMID_T>& vec = termList->terms();
	    for (size_t i = 0; i < vec.size(); ++i) ++bag[vec[i]];

	    int DL = vec.size();
	    unique_ptr<indri_contrib::DocumentCentricTermScoreFunction> scoreFunction(
		getDocumentCentricTermScoreFunction(rule, DL, N, D));

	    prob.clear();
	    for (unordered_map<lemur::api::TERMID_T, int>::const_iterator iter = 
		 bag.begin(); iter != bag.end(); ++iter)
	    {
		lemur::api::TERMID_T seqID = term2seq[iter->first];
		if (seqID == 0) continue; // FIXME: a quick hack for removing OOVs

		int TF = iter->second;
		double probScore = scoreFunction->compute(TF, CF_cache[seqID], DF_cache[seqID]);
		prob.push_back(make_pair(seqID, probScore));
	    }

	    if (!scoreFunction->isProbability()) {
		double sum = 0.0;
		for (size_t i = 0; i < prob.size(); ++i) {
		    if (prob[i].second < 0) prob[i].second = 0;
		    sum += prob[i].second;
		}

		for (size_t i = 0; i < prob.size(); ++i)
		    prob[i].second = scale * prob[i].second / sum;

		if (softmax) {
		    sum = 0.0;
		    for (size_t i = 0; i < prob.size(); ++i) {
			prob[i].second = std::exp(prob[i].second);
			sum += prob[i].second;
		    }

		    for (size_t i = 0; i < prob.size(); ++i) 
			prob[i].second /= sum;
		}
	    }

	    if (guardLoader)
		guardLoader->next(guards);

	    if (!prob.empty()) {
		double cumprob = 0.0;
		int offset = 0;

		if (!guards.empty()) {
		    for (size_t i = 0; i < prob.size(); ++i) {
			if (guards.find(prob[i].first) != guards.end()) {
			    swap(prob[offset], prob[i]);
			    cumprob += prob[offset++].second;
			}
		    }
		}

		stable_sort(prob.begin() + offset, prob.end(), 
			    greater_value<lemur::api::TERMID_T, double>());

		if (step_transform) {
		    double step = 1.0 / prob.size();
		    double sum = 0.5 * (prob.size() + 1);
		    for (size_t i = 0; i < prob.size(); ++i)
			prob[i].second = (1.0 - i * step) / sum;
		    cumprob = 0.0;
		    offset = 0;
		} 
		else if (uniform_transform) {
		    double step = 1.0 / prob.size();
		    for (size_t i = 0; i < prob.size(); ++i)
			prob[i].second = step;
		    cumprob = 0.0;
		    offset = 0;
		}

		for (size_t i = offset; i < prob.size(); ++i) {
		    cumprob += prob[i].second;
		    // prob[i].second = float(1.0) / cumprob;
		    prob[i].second = power_function(cumprob);
		}

		if (sampler)
		    for (size_t i = 0; i < prob.size(); ++i)  // FIXME: better go from offset instead of 0?
			sampler->add(prob[i].second);

		if (filter) {
		    for (size_t i = 0; i < offset; ++i) {
			lemur::api::TERMID_T seqID = prob[i].first;
			filter->operator[](seqID - 1).push_back(false);
		    }

		    for (size_t i = offset; i < prob.size(); ++i) {
			lemur::api::TERMID_T seqID = prob[i].first;
			double weighted_gain = prob[i].second;

			if (seqID)
			    filter->operator[](seqID - 1).push_back(weighted_gain < epsilon);
		    }
		}
	    }
	}
    }
}

// Document-centric pruning
//
void document_centric_pruning(indri::index::Index& index, 
			      indri::api::Parameters& parameters, 
			      Weight& weight, 
			      indri_contrib::QuantileSampler* sampler, 
			      indri_contrib::BinnedFilter* filter)
{
    string rule = parameters.get("rule", "method:tfidf");
    double lambda = parameters.get("epsilon", 0.0); // NOTE: lambda <=> epsilon

    // There are totally N term occurrences in D documents
    long N = index.termCount();
    int D = index.documentCount();

    vector<int> CF_cache;
    CF_cache.push_back(0);
    vector<int> DF_cache;
    DF_cache.push_back(0);

    int T = index.uniqueTermCount();
    vector<int> term2seq(T + 1); // map termID to seq # in inverted list
    term2seq[0] = 0;
    // unordered_map<string, int> pos;

    // Pass 1: Term-based scan, for setting up (empty) bitmaps
    {
	unique_ptr<indri::index::DocListFileIterator>
	    docListFileIterator(index.docListFileIterator());
	int seqID = 0;

	for (docListFileIterator->startIteration();
	     !docListFileIterator->finished();
	     docListFileIterator->nextEntry())
	{
	    indri::index::DocListFileIterator::DocListData* docListData = docListFileIterator->currentEntry();
	    assert(docListData);

	    indri::index::TermData* termData = docListData->termData;
	    assert(termData);

	    indri::index::DocListIterator* docListIterator = docListData->iterator;
	    assert(docListIterator);

	    ++seqID;
	    term2seq[index.term(termData->term)] = seqID;

	    int CF = termData->corpus.totalCount; // collection term frequency
	    CF_cache.push_back(CF);
	    
	    int DF = termData->corpus.documentCount; // document frequency
	    DF_cache.push_back(DF);

	    while (!docListIterator->finished())
		docListIterator->nextEntry(D + 1);

	    if (filter) {
		boost::dynamic_bitset<unsigned long> bitmap; 
		filter->add(bitmap);
	    }
	}
    }

    // Pass 2: Document-based scan, to compute scores
    {
	unique_ptr<indri::index::TermListFileIterator>
	    termListFileIterator(index.termListFileIterator());

	unordered_map<lemur::api::TERMID_T, int> bag;
	vector<pair<lemur::api::TERMID_T, double>> prob;
	lemur::api::DOCID_T docID = 0;
	
	for (termListFileIterator->startIteration();
	     !termListFileIterator->finished();
	     termListFileIterator->nextEntry())
	{
	    ++docID;

	    indri::index::TermList* termList = termListFileIterator->currentEntry();
	    assert(termList);
	
	    bag.clear();
	    const indri::utility::greedy_vector<lemur::api::TERMID_T>& vec = termList->terms();
	    for (size_t i = 0; i < vec.size(); ++i) ++bag[vec[i]];

	    int DL = vec.size();
	    unique_ptr<indri_contrib::DocumentCentricTermScoreFunction> scoreFunction(
		getDocumentCentricTermScoreFunction(rule, DL, N, D));
	
	    prob.clear();
	    for (unordered_map<lemur::api::TERMID_T, int>::const_iterator iter = 
		 bag.begin(); iter != bag.end(); ++iter)
	    {
		lemur::api::TERMID_T seqID = term2seq[iter->first];
		if (seqID == 0) continue; // FIXME: a quick hack for removing OOVs

		int TF = iter->second;
		double probScore = scoreFunction->compute(TF, CF_cache[seqID], DF_cache[seqID]);
		prob.push_back(make_pair(seqID, probScore));
	    }

	    if (!prob.empty()) {
		stable_sort(prob.begin(), prob.end(), 
			    greater_value<lemur::api::TERMID_T, double>());

		double step = 1.0 / prob.size();
		for (size_t i = 0; i < prob.size(); ++i)
		    prob[i].second = 1.0 - i * step; // comparable to lambda

		if (sampler) {
		    for (size_t i = 0; i < prob.size(); ++i) 
			sampler->add(prob[i].second);
		}

		if (filter) {
		    for (size_t i = 0; i < prob.size(); ++i) {
			lemur::api::TERMID_T seqID = prob[i].first;
			double l = prob[i].second;
			filter->operator[](seqID - 1).push_back(l <= lambda); // note the ``equality''
		    }
		}
	    }
	}
    }
}

template<typename T> void _loadVectorFile(std::vector<T>& v, std::string path)
{
    std::ifstream in(path);
    T val;
    while (in >> val) 
	v.push_back(val);
}

// Popularity-based pruning
//
void popularity_based_pruning(indri::index::Index& index,
			      indri::api::Parameters& parameters,
			      Weight& weight,
			      indri_contrib::QuantileSampler* sampler,
			      indri_contrib::BinnedFilter* filter)
{
    double epsilon = parameters.get("epsilon", 0.0);

    if (!parameters.exists("termGain"))
	die("Need to specify the path to term gain file (-termGain=PATH)");

    string termGainPath = parameters["termGain"];
    vector<double> termGain;
    termGain.push_back(0); // Orz...

    _loadVectorFile<double>(termGain, termGainPath);

    // long N = index.termCount();
    // int D = index.documentCount();

    unique_ptr<indri::index::DocListFileIterator>
	docListFileIterator(index.docListFileIterator());

    int seqID = 0;

    for (docListFileIterator->startIteration();
	 !docListFileIterator->finished();
	 docListFileIterator->nextEntry())
    {
	indri::index::DocListFileIterator::DocListData* docListData = docListFileIterator->currentEntry();
	assert(docListData);

	indri::index::TermData* termData = docListData->termData;
	assert(termData);

	indri::index::DocListIterator* docListIterator = docListData->iterator;
	assert(docListIterator);

	++seqID;

	// int CF = termData->corpus.totalCount; // collection term frequency
	int DF = termData->corpus.documentCount; // document frequency

	vector<double> scores;

	for (; !docListIterator->finished(); docListIterator->nextEntry()) {
	    indri::index::DocListIterator::DocumentData* documentData = 
		docListIterator->currentEntry();

	    int docID = documentData->document;
	    scores.push_back(termGain[seqID]);
	}

	if (sampler) 
	    sampler->add(scores.begin(), scores.end());

	if (filter) {
	    int offset = 0;
	    boost::dynamic_bitset<unsigned long> bitmap(DF); // allocate enough space first

	    BOOST_FOREACH (double score, scores) {
		if (score < epsilon)
		    bitmap.set(offset);
		++offset;
	    }

	    bitmap.resize(offset);
	    filter->add(bitmap);
	}
    }
}


// Two-sample two-proportion test (2N2P)
//
void two_proportion_test(indri::index::Index& index, 
			 indri::api::Parameters& parameters, 
			 Weight& weight, 
			 indri_contrib::QuantileSampler* sampler, 
			 indri_contrib::BinnedFilter* filter) 
{
    if (parameters.exists("rule"))
	die("Argument -rule has no effect in prp_based_pruning");

    double epsilon = parameters.get("epsilon", 0.0);
    double h = parameters.get("h", 0.0);

    // There are totally N term occurrences in D documents
    long N = index.termCount();
    int D = index.documentCount();

    unique_ptr<indri::index::DocListFileIterator> 
	docListFileIterator(index.docListFileIterator());

    for (docListFileIterator->startIteration(); 
	 !docListFileIterator->finished(); 
	 docListFileIterator->nextEntry()) 
    {
	indri::index::DocListFileIterator::DocListData* docListData = docListFileIterator->currentEntry();
	assert(docListData);

	indri::index::TermData* termData = docListData->termData;
	assert(termData);

	indri::index::DocListIterator* docListIterator = docListData->iterator;
	assert(docListIterator);

	int CF = termData->corpus.totalCount; // collection term frequency
	int DF = termData->corpus.documentCount; // document frequency

	double p_c = double(CF) / N;

	vector<double> scores;

	for (; !docListIterator->finished(); docListIterator->nextEntry()) {
	    indri::index::DocListIterator::DocumentData* documentData = 
		docListIterator->currentEntry();

	    int docID = documentData->document;
	    int TF = documentData->positions.size(); // term frequency (in this document)
	    int DL = index.documentLength(docID);

	    double p_d = double(TF) / DL;
	    double P = double(TF + CF) / (DL + N);
	    double scale = (1.0 / DL) + (1.0 / N);

	    double score = (p_d + p_c) / 
		std::sqrt(P * (1 - P) * scale);
	    if (h > 0.0) 
		score -= h * std::sqrt(scale);

	    scores.push_back(score);
	}

	if (sampler) 
	    sampler->add(scores.begin(), scores.end());

	if (filter) {
	    int offset = 0;
	    boost::dynamic_bitset<unsigned long> bitmap(DF); // allocate enough space first

	    BOOST_FOREACH (double score, scores) {
		if (score < epsilon)
		    bitmap.set(offset);
		++offset;
	    }

	    bitmap.resize(offset);
	    filter->add(bitmap);
	}
    }
}


// Probability-ranking principle 
//
void prp_based_pruning(indri::index::Index& index, 
		       indri::api::Parameters& parameters, 
		       Weight& weight, 
		       indri_contrib::QuantileSampler* sampler, 
		       indri_contrib::BinnedFilter* filter) 
{
    if (parameters.exists("rule"))
	die("Argument -rule has no effect in prp_based_pruning");

    double epsilon = parameters.get("epsilon", 0.0);

    // There are totally N term occurrences in D documents
    long N = index.termCount();
    int D = index.documentCount();

    RunningStat stat;
    for (int i = 1; i <= D; ++i)
	stat.push(index.documentLength(i));

    double dlMu = stat.mean();
    double dlSigma = stat.standardDeviation();

    unique_ptr<indri::index::DocListFileIterator> 
	docListFileIterator(index.docListFileIterator());

    for (docListFileIterator->startIteration(); 
	 !docListFileIterator->finished(); 
	 docListFileIterator->nextEntry()) 
    {
	indri::index::DocListFileIterator::DocListData* docListData = docListFileIterator->currentEntry();
	assert(docListData);

	indri::index::TermData* termData = docListData->termData;
	assert(termData);

	indri::index::DocListIterator* docListIterator = docListData->iterator;
	assert(docListIterator);

	int CF = termData->corpus.totalCount; // collection term frequency
	int DF = termData->corpus.documentCount; // document frequency

	// probability of term given nonrelevance p_q = p(q|C)
	double p_q = static_cast<double>(CF) / N;

	vector<double> scores;

	for (; !docListIterator->finished(); docListIterator->nextEntry()) {
	    indri::index::DocListIterator::DocumentData* documentData = 
		docListIterator->currentEntry();

	    int docID = documentData->document;
	    int TF = documentData->positions.size(); // term frequency (in this document)
	    int DL = index.documentLength(docID);

	    // document prior p_d = p(r|D)
	    double p_d = 0.5 + 0.1 * std::tanh(double(DL - dlMu) / dlSigma);

	    // query likelihood p_r = p(q|D,r)
	    double p_r = 0.4 * double(TF) / DL + 0.6 * p_q;

	    double score = (p_r / p_q) * (p_d / (1.0 - p_d));
	    scores.push_back(score);
	}

	if (sampler) 
	    sampler->add(scores.begin(), scores.end());

	if (filter) {
	    int offset = 0;
	    boost::dynamic_bitset<unsigned long> bitmap(DF); // allocate enough space first

	    BOOST_FOREACH (double score, scores) {
		if (score < epsilon)
		    bitmap.set(offset);
		++offset;
	    }

	    bitmap.resize(offset);
	    filter->add(bitmap);
	}
    }
}


// Uniform pruning
//
void uniform_pruning(indri::index::Index& index, 
		     indri::api::Parameters& parameters, 
		     Weight& weight,
		     indri_contrib::QuantileSampler* sampler,
		     indri_contrib::BinnedFilter* filter) 
{
    string rule = parameters.get("rule", "method:tfidf");
    double epsilon = parameters.get("epsilon", 0.0);

    // There are totally N term occurrences in D documents
    long N = index.termCount();
    int D = index.documentCount();

    unique_ptr<indri::index::DocListFileIterator> 
	docListFileIterator(index.docListFileIterator());

    for (docListFileIterator->startIteration(); 
	 !docListFileIterator->finished(); 
	 docListFileIterator->nextEntry()) 
    {
	indri::index::DocListFileIterator::DocListData* docListData = docListFileIterator->currentEntry();
	assert(docListData);

	indri::index::TermData* termData = docListData->termData;
	assert(termData);

	indri::index::DocListIterator* docListIterator = docListData->iterator;
	assert(docListIterator);

	int CF = termData->corpus.totalCount; // collection term frequency
	int DF = termData->corpus.documentCount; // document frequency

	unique_ptr<indri::query::TermScoreFunction> scoreFunction(
	    indri::query::TermScoreFunctionFactory::get(rule, CF, N, DF, D));
	vector<double> scores;

	for (; !docListIterator->finished(); docListIterator->nextEntry()) {
	    indri::index::DocListIterator::DocumentData* documentData = 
		docListIterator->currentEntry();

	    int docID = documentData->document;
	    int TF = documentData->positions.size(); // term frequency (in this document)
	    int DL = index.documentLength(docID);

	    double logWeight = std::log(weight[docID]); // always a catch!
	    scores.push_back(logWeight + scoreFunction->scoreOccurrence(TF, DL));
	}

	if (sampler) 
	    sampler->add(scores.begin(), scores.end());

	if (filter) {
	    int offset = 0;
	    boost::dynamic_bitset<unsigned long> bitmap(DF); // allocate enough space first

	    BOOST_FOREACH (double score, scores) {
		if (score < epsilon)
		    bitmap.set(offset);
		++offset;
	    }

	    bitmap.resize(offset);
	    filter->add(bitmap);
	}
    }
}


// Term-based pruning
//
void term_based_pruning(indri::index::Index& index, 
			indri::api::Parameters& parameters, 
			Weight& weight, 
			indri_contrib::QuantileSampler* sampler, 
			indri_contrib::BinnedFilter* filter)
{
    string rule = parameters.get("rule", "method:tfidf");
    int k = parameters.get("k", 10);
    double epsilon = parameters.get("epsilon", 0.0);

    // There are totally N term occurrences in D documents
    long N = index.termCount();
    int D = index.documentCount();

    unique_ptr<indri::index::DocListFileIterator> 
	docListFileIterator(index.docListFileIterator());

    for (docListFileIterator->startIteration(); 
	 !docListFileIterator->finished(); 
	 docListFileIterator->nextEntry())
    {
	indri::index::DocListFileIterator::DocListData* docListData = 
	    docListFileIterator->currentEntry();
	indri::index::TermData* termData = docListData->termData;
	indri::index::DocListIterator* docListIterator = docListData->iterator;

	int CF = termData->corpus.totalCount; // collection term frequency
	int DF = termData->corpus.documentCount; // document frequency

	unique_ptr<indri::query::TermScoreFunction> scoreFunction(
	    indri::query::TermScoreFunctionFactory::get(rule, CF, N, DF, D));
	vector<double> scores;

	for (; !docListIterator->finished(); docListIterator->nextEntry()) {
	    indri::index::DocListIterator::DocumentData* documentData = 
		docListIterator->currentEntry();

	    int TF = documentData->positions.size(); // term frequency (in this document)
	    int DL = index.documentLength(documentData->document);
	    scores.push_back(scoreFunction->scoreOccurrence(TF, DL));
	}

	if (sampler) {
	    if (scores.size() > k) {
		vector<double> sorter(scores.begin(), scores.end());
		nth_element(sorter.begin(), sorter.begin() + k - 1, sorter.end(), greater<double>());

		// NOTE: This is just a hack for skipping negative scores
		//       (only works for okapi)
		double denom = sorter[k - 1];
		if (denom > 0) { 
		    transform(sorter.begin(), sorter.end(), sorter.begin(), 
			      boost::bind(divides<double>(), _1, denom));

		    sampler->add(sorter.begin(), sorter.end());
		}
	    }
	    else {
		double fake_score = 1.0 + 1.0 / k;
		sampler->add(scores.size(), fake_score);
	    }
	}
	
	if (filter) {
	    int offset = 0;
	    boost::dynamic_bitset<unsigned long> bitmap(DF); // allocate enough space first

	    if (scores.size() > k) {
		vector<double> sorter(scores.begin(), scores.end());
		nth_element(sorter.begin(), sorter.begin() + k - 1, sorter.end(), greater<double>());

		double z = sorter[k - 1] * epsilon;
		BOOST_FOREACH (double score, scores) {
		    if (score < z) 
			bitmap.set(offset);
		    ++offset;
		}
	    }
	    else 
		offset += scores.size();

	    bitmap.resize(offset);
	    filter->add(bitmap);
	}
    }
}


//--------------------------------------------------
// Handlers
//-------------------------------------------------- 
void compile_guard(indri::api::Parameters& parameters)
{
    using indri::file::Path;

    std::vector<std::string> requiredKeys = { "index", "output", "rule", "k" };
    requires(parameters, requiredKeys);

    std::string sourcePath = parameters["index"];
    std::string outputPath = parameters["output"];
    std::string rule = parameters.get("rule", "method:tfidf");
    int k = parameters.get("k", 10);

    indri::collection::Repository repo;
    repo.openRead(sourcePath);
    indri::index::Index& index = *((*repo.indexes())[0]);

    // There are totally N term occurrences in D documents
    long N = index.termCount();
    int D = index.documentCount();

    std::vector<std::vector<int>> inv(D + 1, std::vector<int>());
    int seqID = 0;

    unique_ptr<indri::index::DocListFileIterator> 
	docListFileIterator(index.docListFileIterator());

    for (docListFileIterator->startIteration(); 
	 !docListFileIterator->finished(); 
	 docListFileIterator->nextEntry())
    {
	++seqID;

	indri::index::DocListFileIterator::DocListData* docListData = 
	    docListFileIterator->currentEntry();
	indri::index::TermData* termData = docListData->termData;
	indri::index::DocListIterator* docListIterator = docListData->iterator;

	int CF = termData->corpus.totalCount; // collection term frequency
	int DF = termData->corpus.documentCount; // document frequency

	unique_ptr<indri::query::TermScoreFunction> scoreFunction(
	    indri::query::TermScoreFunctionFactory::get(rule, CF, N, DF, D));
	vector<double> scores;
	vector<lemur::api::DOCID_T> docids;

	for (; !docListIterator->finished(); docListIterator->nextEntry()) {
	    indri::index::DocListIterator::DocumentData* documentData = 
		docListIterator->currentEntry();

	    int TF = documentData->positions.size(); // term frequency (in this document)
	    int DL = index.documentLength(documentData->document);
	    scores.push_back(scoreFunction->scoreOccurrence(TF, DL));
	    docids.push_back(documentData->document);
	}

	if (scores.size() > k) {
	    vector<double> sorter(scores.begin(), scores.end());
	    nth_element(sorter.begin(), sorter.begin() + k - 1, sorter.end(), greater<double>());

	    double z = sorter[k - 1];
	    for (size_t i = 0; i < scores.size(); ++i)
		if (scores[i] >= z)
		    inv[docids[i]].push_back(seqID);
	}
	else
	    for (size_t i = 0; i < docids.size(); ++i)
		inv[docids[i]].push_back(seqID);
    }

    std::ofstream out(outputPath.c_str());

    for (size_t i = 0; i < inv.size(); ++i) {
	for (size_t j = 0; j < inv[i].size(); ++j)
	    out << inv[i][j] << ' ';
	out << '\n';
    }
}

void compile_access(indri::api::Parameters& parameters) {
    using indri::file::Path;

    std::vector<std::string> requiredKeys = { "index", "viewFile", "output" };
    requires(parameters, requiredKeys);

    std::string sourcePath = parameters["index"];
    std::string viewPath = parameters["viewFile"];
    std::string outputPath = parameters["output"];

    indri::collection::Repository repo;
    repo.openRead(sourcePath);
    indri::collection::CompressedCollection& collection = *repo.collection();

    std::unordered_map<std::string, int> access;

    {
	std::ifstream vin(viewPath.c_str());
	std::string docno;

	while (vin >> docno) ++access[docno];
    }

    {
	std::ofstream out(outputPath.c_str());

	int D = (*repo.indexes())[0]->documentCount();
	for (size_t docID = 1; docID <= D; ++docID) {
	    std::string docno = collection.retrieveMetadatum(docID, "docno");

	    if (access.find(docno) != access.end())
		out << access[docno] << '\n';
	    else
		out << 0 << '\n';
	}
    }
}


void compile_term_gain(indri::api::Parameters& parameters) {
    using indri::file::Path;

    std::vector<std::string> requiredKeys = { "index", "queryFile", "hasCounts", "output" };
    requires(parameters, requiredKeys);

    std::string sourcePath = parameters["index"];
    std::string queryPath = parameters["queryFile"];
    std::string outputPath = parameters["output"];
    bool hasCounts = parameters.get("hasCounts", false);
    bool debug = parameters.get("debug", false);

    indri::collection::Repository repo;
    repo.openRead(sourcePath);
    indri::index::Index& index = *((*repo.indexes())[0]);

    // Step 1: accumulate stem counts
    std::unordered_map<std::string, int> stemPop;

    {
	std::unordered_map<std::string, int> bag;
	std::string term;

	std::ifstream in(queryPath.c_str());

	if (hasCounts) {
	    int queryCount;
	    std::string line;
	    std::istringstream line_in;

	    while (std::getline(in, line)) {
		line_in.str(line);
		line_in.clear();

		line_in >> queryCount;
		while (line_in >> term)
		    bag[term] += queryCount;
	    }
	}
	else
	    while (in >> term) ++bag[term];

	for (std::unordered_map<std::string, int>::iterator iter = bag.begin();
	     iter != bag.end(); ++iter)
	{
	    std::string stem = repo.processTerm(iter->first);
	    stemPop[stem] += iter->second;
	}
    }

    std::vector<double> termGain;
    std::vector<std::string> stems;
    std::vector<int> stemCounts;

    termGain.push_back(0);
    stems.push_back("");
    stemCounts.push_back(0);

    int D = index.documentCount();

    // Step 2: Term-based scan, for setting up (empty) bitmaps
    {
	unique_ptr<indri::index::DocListFileIterator>
	    docListFileIterator(index.docListFileIterator());
    
	for (docListFileIterator->startIteration();
	     !docListFileIterator->finished();
	     docListFileIterator->nextEntry())
	{
	    indri::index::DocListFileIterator::DocListData* docListData = docListFileIterator->currentEntry();
	    assert(docListData);
    
	    indri::index::TermData* termData = docListData->termData;
	    assert(termData);

	    int DF = termData->corpus.documentCount; // document frequency
    
	    std::unordered_map<std::string, int>::iterator match = 
		stemPop.find(termData->term);

	    stems.push_back(termData->term);

	    if (match == stemPop.end()) {
		termGain.push_back(0);
		stemCounts.push_back(0);
	    }
	    else {
		termGain.push_back(static_cast<double>(match->second) / DF);
		stemCounts.push_back(match->second);
	    }

	    // Still needs to traverse to the end
	    indri::index::DocListIterator* docListIterator = docListData->iterator;
	    assert(docListIterator);
    
	    while (!docListIterator->finished())
		docListIterator->nextEntry(D + 1);
	}
    }

    // Step 3: Output
    {
	std::ofstream out(outputPath.c_str());

	if (debug)
	    for (size_t i = 1; i < termGain.size(); ++i)
		out << boost::format("%3.17f\t%d\t%s") % termGain[i] % stemCounts.at(i) % stems.at(i) << '\n';
	else
	    for (size_t i = 1; i < termGain.size(); ++i) 
		out << boost::format("%3.17f") % termGain[i] << '\n';
    }
}

void compile_query_view(indri::api::Parameters& parameters) {
    using indri::file::Path;

    std::vector<std::string> requiredKeys = { "index", "queryFile", "viewFile", "output" };
    requires(parameters, requiredKeys);

    std::string sourcePath = parameters["index"];
    std::string queryPath = parameters["queryFile"];
    std::string viewPath = parameters["viewFile"];
    std::string outputPath = parameters["output"];

    indri::collection::Repository repo;
    repo.openRead(sourcePath);
    indri::index::Index& index = *((*repo.indexes())[0]);
    indri::collection::CompressedCollection& collection = *repo.collection();

    std::unordered_map<std::string, std::vector<lemur::api::DOCID_T>> post;

    {
	std::unordered_map<std::string, std::string> term2stem;

	std::ifstream qin(queryPath.c_str()), vin(viewPath.c_str());
	std::string line;
	std::istringstream iss;
	std::vector<std::string> stems;
	std::vector<std::string> docnos;
	std::vector<lemur::api::DOCID_T> docids;

	while (true) {

	    // first let us retrieve the stems
	    if (!getline(qin, line)) break;

	    boost::trim(line);

	    stems.clear();
	    boost::split(stems, line, boost::is_any_of(" "));
	    for (size_t i = 0; i < stems.size(); ++i) {
		if (term2stem.find(stems[i]) == term2stem.end())
		    term2stem[stems[i]] = repo.processTerm(stems[i]);
		stems[i] = term2stem[stems[i]];
	    }

	    // then we look at the document views
	    if (!getline(vin, line)) break;

	    boost::trim(line);

	    docnos.clear();
	    docids.clear();
	    boost::split(docnos, line, boost::is_any_of(" "));
	    for (size_t i = 0; i < docnos.size(); ++i) {
		std::vector<lemur::api::DOCID_T> ids = 
		    collection.retrieveIDByMetadatum("docno", docnos[i]);
		docids.push_back(ids[0]);
	    }

	    // save this to the hashtable
	    for (size_t i = 0; i < stems.size(); ++i) 
		for (size_t j = 0; j < docids.size(); ++j) 
		    post[stems[i]].push_back(docids[j]);
	}
    }


    std::ofstream out(outputPath.c_str());

    for (std::unordered_map<std::string, std::vector<lemur::api::DOCID_T>>::iterator iter = post.begin();
	 iter != post.end(); ++iter)
    {
	std::vector<lemur::api::DOCID_T>& l = iter->second;
	std::sort(l.begin(), l.end());
	std::vector<lemur::api::DOCID_T>::iterator last = std::unique(l.begin(), l.end());
	l.erase(last, l.end());
    }

    int D = index.documentCount();
    
    {
	unique_ptr<indri::index::DocListFileIterator>
	    docListFileIterator(index.docListFileIterator());
    
	for (docListFileIterator->startIteration();
	     !docListFileIterator->finished();
	     docListFileIterator->nextEntry())
	{
	    indri::index::DocListFileIterator::DocListData* docListData = docListFileIterator->currentEntry();
	    assert(docListData);
    
	    indri::index::TermData* termData = docListData->termData;
	    assert(termData);
    
	    int DF = termData->corpus.documentCount; // document frequency

	    indri::index::DocListIterator* docListIterator = docListData->iterator;
	    assert(docListIterator);
    
	    if (post.find(termData->term) != post.end()) {
		std::vector<lemur::api::DOCID_T>& l = post[termData->term];
		std::vector<lemur::api::DOCID_T>::iterator curr = l.begin();

		size_t offset = 0;
		while (!docListIterator->finished() && curr != l.end()) {
		    ++offset;

		    indri::index::DocListIterator::DocumentData* docData = docListIterator->currentEntry();
		    lemur::api::DOCID_T docID = docData->document;

		    if (*curr == docID) {
			out << ' ' << offset;
			++curr;
		    }

		    docListIterator->nextEntry();
		}
	    }

	    out << '\n';
    
	    while (!docListIterator->finished())
		docListIterator->nextEntry(D + 1);
	}
    }
}

void prune_index(indri::api::Parameters& parameters) {
    using indri::file::Path;

    std::vector<std::string> requiredKeys = { "index", "method" };
    requires(parameters, requiredKeys);

    std::string sourcePath = parameters["index"];
    std::string method = parameters["method"];

    struct PruneMethod {
	typedef void (*Handler)(indri::index::Index&, 
				indri::api::Parameters&, 
				Weight&, 
				indri_contrib::QuantileSampler*,
				indri_contrib::BinnedFilter*);

	std::string name;
	Handler handler;
    };

    std::vector<PruneMethod> listOfPruneMethods = {
	{ "gup_kl", generalized_uniform_pruning<LogarithmGain> },
	{ "gup_var", generalized_uniform_pruning<IdentityGain> },
	{ "gup_hellinger", generalized_uniform_pruning<SquareRootGain> },
	{ "gup_x2", generalized_uniform_pruning<MinusReciprocalGain> },
	{ "gup_renyi", generalized_uniform_pruning<MinusExponentialGain> },

	{ "guarded_x2", generalized_uniform_pruning<MinusReciprocalGain> },

	{ "mix_kl_x2", generalized_uniform_pruning<MixtureGain<LogarithmGain, MinusReciprocalGain>> },
	{ "mix_hellinger_x2", generalized_uniform_pruning<MixtureGain<SquareRootGain, MinusReciprocalGain>> },
	{ "mix_renyi_x2", generalized_uniform_pruning<MixtureGain<MinusExponentialGain, MinusReciprocalGain>> },

	{ "dcp", document_centric_pruning },
	{ "pp", popularity_based_pruning },
	{ "up", uniform_pruning },
	{ "tcp", term_based_pruning },
	{ "prp", prp_based_pruning },
	{ "2n2p", two_proportion_test },

	{ "infinite_renyi", renyi_divergence_infinite_order }
    };

    PruneMethod::Handler handler = 0;
    BOOST_FOREACH (const PruneMethod& pruneMethod, listOfPruneMethods) {
	if (method == pruneMethod.name) {
	    handler = pruneMethod.handler;
	    break;
	}
    }

    if (!handler) 
	die(boost::format("No such method '%s'") % method);

    // Load up the old index
    indri::index::DiskIndex diskIndex;

    {
	indri::api::Parameters manifest;
	manifest.loadFile(Path::combine(sourcePath, "manifest"));
	diskIndex.open(Path::combine(sourcePath, "index"), 
		       i64_to_string((INT64)manifest["indexes.index"]));
    }

    indri_contrib::BinnedFilter filter;

    int documentCount = diskIndex.documentCount();
    Weight weight(documentCount + 1, 1.0); // uniform weight

    if (parameters.exists("prior")) {
	std::string priorPath = parameters["prior"];
	double beta = parameters.get("beta", 0.0);

	indri::collection::Repository repo;
	repo.openRead(sourcePath);
	indri::collection::CompressedCollection* collection = repo.collection();
	weight.assign(documentCount + 1, 0.0); 

	{
	    std::ifstream in(priorPath.c_str());
	    std::vector<lemur::api::DOCID_T> result;
	    double cumsum = 0.0;

	    std::string docno;
	    double prob;
	    while (in >> docno >> prob) {
		result = collection->retrieveIDByMetadatum("docno", docno);
		assert(result.size() == 1);
		weight[result.front()] = prob;
		cumsum += prob;
	    }

	    for (int i = 1; i <= documentCount; ++i) {
		if (weight[i] == 0.0) weight[i] = 1.0;
		else weight[i] = 1.0 + weight[i] / (1.0 - cumsum + beta) * documentCount;
	    }
	}
    }

    if (parameters.exists("estimate")) {
	indri::api::Parameters estimateSpec;
	parseSpec(estimateSpec, parameters.get("estimate"));

	std::string method = estimateSpec.get("method", "p2"); // default using p-square algorithm
	std::vector<double> poi = {
	    0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
	    0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99 };

	// for (double f = 0.1; f <= 0.9; f += 0.1) poi.push_back(f);
	// for (double f = 0.91; f <= 0.99; f += 0.01) poi.push_back(f);

	// Method 1: Bisection
	if (method == "bisection") {
	    double lb = estimateSpec.get("lb", 0.0);
	    double ub = estimateSpec.get("ub", 1.0);

	    std::unordered_map<double, double> estimates;

	    struct Bracket {
		double a, b;
		double f_a, f_b;
		std::vector<double> val;
	    };

	    Bracket root;
	    root.a = lb;
	    root.b = ub;
	    root.f_a = 0;
	    root.f_b = 1;
	    root.val = poi;

	    std::queue<Bracket> q;
	    q.push(root);

	    while (!q.empty()) {
		Bracket r = q.front();
		q.pop();

		std::cerr << boost::format("(%3.17f, %3.17f) -> (%3.17f, %3.17f): ") % r.a % r.b % r.f_a % r.f_b;
		BOOST_FOREACH (double v, r.val) {
		    std::cerr << v << ' ';
		}
		std::cerr << '\n';

		double p = (r.a + r.b) / 2;
		parameters.set("epsilon", p);
		handler(diskIndex, parameters, weight, 0, &filter);

		double f_p = filter.ratio();
		filter.clear();

		std::cerr << boost::format("  epsilon = %3.17f, rho = %3.17f") % p % f_p << '\n';

		const double SIGMA = 0.000000001; // change in x
		const double DELTA = 0.000001; // change in f(x)

		if (r.b - r.a <= SIGMA) {
		    estimates[p] = f_p; 
		    BOOST_FOREACH (double v, r.val) {
			std::cerr << boost::format("  DONE %f") % v << '\n';
		    }
		}
		else {
		    Bracket r_ap, r_pb;
		    r_ap.a = r.a; r_ap.f_a = r.f_a;
		    r_ap.b = p; r_ap.f_b = f_p;
		    r_pb.a = p; r_pb.f_a = f_p;
		    r_pb.b = r.b; r_pb.f_b = r.f_b;

		    BOOST_FOREACH (double v, r.val) {
			if (v >= f_p - DELTA && v <= f_p + DELTA) {
			    estimates[p] = f_p; 
			    std::cerr << boost::format("  DONE %f") % v << '\n';
			}
			else if (v < f_p - DELTA)
			    r_ap.val.push_back(v);
			else if (v > f_p + DELTA)
			    r_pb.val.push_back(v);
		    }

		    if (!r_ap.val.empty()) q.push(r_ap);
		    if (!r_pb.val.empty()) q.push(r_pb);
		}
	    }

	    std::vector<double> keys;
	    for (std::unordered_map<double, double>::iterator iter = 
		 estimates.begin(); iter != estimates.end(); ++iter) 
		keys.push_back(iter->first);

	    std::sort(keys.begin(), keys.end());
	    BOOST_FOREACH (double key, keys) {
		std::cout << boost::format("%3.17f\t%3.17f") % key % estimates[key] << '\n';
	    }
	}
	// Case 2: P-square algorithm
	else if (method == "p_square" || method == "p2") {
	    double step = estimateSpec.get("step", 0.1);
	    double window = estimateSpec.get("window", 0.0);

	    std::vector<double> points;

	    if (step > 0)
		for (double q = step; q < 1.0; q += step) {
		    if (window > 0) points.push_back(q - window);

		    points.push_back(q);

		    if (window > 0) points.push_back(q + window);
		}

	    std::stable_sort(points.begin(), points.end());

	    indri_contrib::PSquareQuantileSampler sampler(points);
	    handler(diskIndex, parameters, weight, &sampler, 0); 

	    std::vector<double> estimates = sampler.query(poi);
	    for (size_t i = 0; i < poi.size(); ++i)
		std::cout << boost::format("%3.17f\t%3.17f") % estimates[i] % poi[i] << '\n';
	}
	// Case 3: Reservoir sampling
	else if (method == "reservoir") {
	    size_t reservoirSize = estimateSpec.get("size", 1000000);

	    indri_contrib::ReservoirQuantileSampler sampler(reservoirSize);
	    handler(diskIndex, parameters, weight, &sampler, 0); 

	    std::vector<double> estimates = sampler.query(poi);
	    for (size_t i = 0; i < poi.size(); ++i)
		std::cout << boost::format("%3.17f\t%3.17f") % estimates[i] % poi[i] << '\n';
	}
	else 
	    die(boost::format("No such estimation method '%s'") % method);
    }
    else {
	// The pruning method returns a prune (filter)
	handler(diskIndex, parameters, weight, 0, &filter);
	cout << filter.ratio() << '\n';

	indri_contrib::FilteredIndex<indri_contrib::BinnedFilter> 
	    filteredIndex(diskIndex, filter);

	// Apply and produce a new index ONLY when the name is given
	if (parameters.exists("newindex")) {
	    string path = parameters["newindex"];
	    if (!Path::exists(path))
		Path::create(path);

	    // Load up the manifest 
	    indri::api::Parameters manifest;
	    manifest.loadFile(Path::combine(sourcePath, "manifest"));

	    if (manifest["indexes.index"].size() > 1)
		die("Source contains more than one index (run 'dump-index compact' first)");

	    // Load up the deletedList
	    indri::index::DeletedDocumentList deletedList;
	    deletedList.read(Path::combine(sourcePath, "deleted"));
	    
	    if (deletedList.deletedCount() > 0)
		die("Source has non-empty deleted list (run 'dump-index compact' first)");

	    // Write out the index
	    try {
		indri::index::IndexWriter writer;
		string indexPath = Path::combine(path, "index");
		if (!Path::exists(indexPath)) Path::create(indexPath);

		vector<indri::index::Index::FieldDescription> fieldDescriptions;
		writer.write(filteredIndex, fieldDescriptions, deletedList, Path::combine(indexPath, "0"));
	    }
	    catch (lemur::api::Exception& e) {
		// LEMUR_RETHROW(e, "sip");
		std::cout << e.what() << std::endl;
	    }

	    // Copy the compressed collection
	    {
		indri::collection::CompressedCollection sourceCollection, collection; 
		sourceCollection.openRead(Path::combine(sourcePath, "collection"));
		
		string collectionPath = Path::combine(path, "collection");
		Path::create(collectionPath);
		collection.create(collectionPath,
				  sourceCollection.forwardFields(),
				  sourceCollection.reverseFields());
		collection.append(sourceCollection, deletedList, 0);
	    }
	    
	    // Write out the deleted list
	    deletedList.write(Path::combine(path, "deleted"));
	    
	    // Write out the manifest
	    manifest["indexes"].set("index", 0);
	    manifest.writeFile(Path::combine(path, "manifest"));
	}
    }
}

//--------------------------------------------------
// Main program
//-------------------------------------------------- 
int main(int argc, char** argv) {
    indri::api::Parameters parameters;

    try {
	parameters.loadCommandLine(argc, argv);
    }
    catch (lemur::api::Exception& e) {
	LEMUR_ABORT(e);
    }

    struct Command {
	typedef void (*Handler)(indri::api::Parameters&);

	std::string name;
	Handler handler;
    };

    std::vector<Command> listOfCommands = {
	{ "pruneIndex", prune_index },
	{ "printInvertedList", print_inverted_list },
	{ "printDirectFile", print_direct_file },
	{ "copyRepository", copy_repository },
	{ "checkIndexHealth", check_index_health },
	{ "checkDocumentPrior", check_document_prior },
	{ "compileGuard", compile_guard },
	{ "compileAccess", compile_access },
	{ "compileTermGain", compile_term_gain },
	{ "compileQueryView", compile_query_view }
    };

    Command::Handler handler = 0;
    BOOST_FOREACH (const Command& cmd, listOfCommands) {
	if (parameters.exists(cmd.name)) {
	    handler = cmd.handler;
	    break;
	}
    }

    if (handler)
	handler(parameters);
    else {
	std::string padding = "  ";
	std::cout << "Available commands:" << '\n';
	BOOST_FOREACH (const Command& cmd, listOfCommands) {
	    std::cout << padding << '-' << cmd.name << '\n';
	}
    }

    return 0;
}
