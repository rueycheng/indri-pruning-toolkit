#include "IndriPruneIndex.h"

#include <unordered_map>
#include <boost/foreach.hpp>

#include "indri/CompressedCollection.hpp"
#include "indri/DiskIndex.hpp"
#include "indri/IndexWriter.hpp"
#include "indri/Path.hpp"

#define __atomic_add __gnu_cxx::__atomic_add
#include "indri/Repository.hpp"

void copy_repository(indri::api::Parameters& parameters) {
    using indri::file::Path;

    vector<string> requiredKeys = { "from", "to" };
    requires(parameters, requiredKeys);

    string sourcePath = parameters["from"];
    string path = parameters["to"];

    // Load up the manifest 
    indri::api::Parameters manifest;
    manifest.loadFile(Path::combine(sourcePath, "manifest"));

    if (manifest["indexes.index"].size() > 1)
	die("Source contains more than one index (run 'dumpindex compact' first)");

    // Load up the deletedList
    indri::index::DeletedDocumentList deletedList;
    deletedList.read(Path::combine(sourcePath, "deleted"));
    
    if (deletedList.deletedCount() > 0)
    	die("Source has non-empty deleted list (run 'dumpindex compact' first)");

    // Load up the index in the ``pruned mode''
    indri::index::DiskIndex diskIndex;
	
    diskIndex.open(Path::combine(sourcePath, "index"),
		   i64_to_string((INT64)manifest["indexes.index"]));
    
    if (diskIndex.documentMaximum() == 0)
	die("The index is empty");

    // Create the destination path as necessary
    if (!Path::exists(path)) Path::create(path);

    // Write out the new index
    {
	indri::index::IndexWriter writer;
	string indexPath = Path::combine(path, "index");
	Path::create(indexPath);

	vector<indri::index::Index::FieldDescription> fieldDescriptions;
	writer.write(diskIndex, fieldDescriptions, deletedList, Path::combine(indexPath, "0"));
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

void check_index_health(indri::index::Index& index, 
			bool dumpPostingList) {
    unique_ptr<indri::index::DocListFileIterator> 
	docListFileIterator(index.docListFileIterator());

    if (dumpPostingList) {
	for (docListFileIterator->startIteration(); 
	     !docListFileIterator->finished(); 
	     docListFileIterator->nextEntry()) 
	{
	    indri::index::DocListFileIterator::DocListData* docListData = docListFileIterator->currentEntry();
	    indri::index::TermData* termData = docListData->termData;
	    indri::index::DocListIterator* docListIterator = docListData->iterator;

	    assert(termData);
	    int documentCount = 0;
	    for (; !docListIterator->finished(); docListIterator->nextEntry()) 
	    {
		indri::index::DocListIterator::DocumentData* documentData = 
		    docListIterator->currentEntry();
		assert(documentData);
		++documentCount;
	    }

	    cout << documentCount << '\n';
	}
    }
    else {
	int termCount = 0;
	long postingCount = 0;

	for (docListFileIterator->startIteration(); 
	     !docListFileIterator->finished(); 
	     docListFileIterator->nextEntry()) 
	{
	    indri::index::DocListFileIterator::DocListData* docListData = docListFileIterator->currentEntry();
	    indri::index::TermData* termData = docListData->termData;
	    indri::index::DocListIterator* docListIterator = docListData->iterator;

	    assert(termData);
	    ++termCount;

	    for (; !docListIterator->finished(); docListIterator->nextEntry()) 
	    {
		indri::index::DocListIterator::DocumentData* documentData = 
		    docListIterator->currentEntry();
		assert(documentData);
		++postingCount;
	    }
	}

	cout << "term count: " << termCount << '\n';
	cout << "posting count: " << postingCount << '\n';
    }
}


void print_direct_file(indri::api::Parameters& parameters) {
    using indri::file::Path;

    vector<string> requiredKeys = { "index", "docID" };
    requires(parameters, requiredKeys);

    string sourcePath = parameters["index"];
    int docID = static_cast<int>(parameters["docID"]);
    bool outputText = parameters.exists("outputText");

    indri::index::DiskIndex diskIndex;
    diskIndex.open(Path::combine(sourcePath, "index"), "0");

    indri::index::Index& original = diskIndex;
    unique_ptr<const indri::index::TermList> termList(original.termList(docID));

    if (outputText) {
	unordered_map<int, string> cache;
	BOOST_FOREACH (int termID, termList->terms()) {
	    unordered_map<int, string>::iterator match = cache.find(termID);
	    if (match != cache.end())
		cout << match->second << '(' << match->first << ')' << ' ';
	    else {
		string term = original.term(termID);
		cache[termID] = term;
		cout << term << '(' << termID << ')' << ' ';
	    }
	}
    }
    else{
	BOOST_FOREACH (int termID, termList->terms()) {
	    cout << termID << ' ';
	}
    }
    cout << endl;
}


// DocListFileIterator manages its own DocListIterators.  Manually deleting the
// iterators (e.g., `delete') can cause segfaults.  Each of these
// DocListIterator has been started (via startIteration()) once they are ready
// in the host object.  So additional calls are also harmful.
//
// Be sure to let these DocListIterators run to the end once initialized.
// Expect segfaults otherwise.
//
void print_inverted_list(indri::api::Parameters& parameters) { 
    vector<string> requiredKeys = { "index", "term" };
    requires(parameters, requiredKeys);

    string sourcePath = parameters["index"];
    string term = parameters["term"];

    indri::index::DiskIndex diskIndex;
    diskIndex.open(indri::file::Path::combine(sourcePath, "index"), "0");

    indri::index::Index& original = diskIndex;
    int documentCount = original.documentCount();

    unique_ptr<indri::index::DocListFileIterator> 
	docListFileIterator(original.docListFileIterator());

    docListFileIterator->startIteration();
    while (!docListFileIterator->finished()) {
	indri::index::DocListFileIterator::DocListData* docListData = docListFileIterator->currentEntry();
	indri::index::TermData* termData = docListData->termData;
	indri::index::DocListIterator* docListIterator = docListData->iterator;

	if (term == termData->term) {
	    while (!docListIterator->finished()) {
		indri::index::DocListIterator::DocumentData* documentData = 
		    docListIterator->currentEntry();
		cout << documentData->document << ' ';
		docListIterator->nextEntry();
	    }
	    cout << endl; 
	}
	else
	    while (!docListIterator->finished()) 
		docListIterator->nextEntry(documentCount + 1);
    
	docListFileIterator->nextEntry();
    }
}


void check_index_health(indri::api::Parameters& parameters) {
    using indri::file::Path;

    std::vector<std::string> requiredKeys = { "index" };
    requires(parameters, requiredKeys);

    std::string sourcePath = parameters["index"];
    bool dumpPostingList = parameters.exists("dump");

    indri::index::DiskIndex diskIndex;
    indri::api::Parameters manifest;
    manifest.loadFile(Path::combine(sourcePath, "manifest"));
    diskIndex.open(Path::combine(sourcePath, "index"), 
		   i64_to_string((INT64)manifest["indexes.index"]));

    check_index_health(diskIndex, dumpPostingList);
}


void check_document_prior(indri::api::Parameters& parameters) {
    using indri::file::Path;

    std::vector<std::string> requiredKeys = { "index", "priorFile", "beta" };
    requires(parameters, requiredKeys);

    std::string sourcePath = parameters["index"];
    std::string priorPath = parameters["priorFile"];
    double beta = parameters["beta"];

    indri::collection::Repository repo;
    repo.openRead(sourcePath);

    indri::index::Index* index = (*repo.indexes())[0];
    indri::collection::CompressedCollection* collection = repo.collection();

    int documentCount = index->documentCount();

    std::vector<double> prior(documentCount + 1, 0.0);
    double cumsum = 0.0;

    {
	std::ifstream in(priorPath.c_str());

	std::vector<lemur::api::DOCID_T> result;
	std::string docno;
	double prob;

	while (in >> docno >> prob) {
	    result = collection->retrieveIDByMetadatum("docno", docno);
	    assert(result.size() == 1);
	    prior[result.front()] = prob;
	    cumsum += prob;
	}
    }

    for (int i = 1; i <= documentCount; ++i) {
	if (prior[i] == 0.0)
	    cout << i << ' ' << '*' << 1.0 << '\n';
	else
	    cout << i << ' ' << ' ' << 1.0 + prior[i] / (1 - cumsum + beta) * documentCount << '\n'; 
    }
}

