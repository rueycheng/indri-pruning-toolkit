#ifndef DOCUMENT_CENTRIC_TERM_SCORE_FUNCTION_H_
#define DOCUMENT_CENTRIC_TERM_SCORE_FUNCTION_H_

namespace indri_contrib {

// DocumentCentricTermScoreFunction
//
class DocumentCentricTermScoreFunction {
public:
    virtual double compute(double termCount, double collectionTermCount, int documentFrequency) = 0;
    virtual double computeLogProb(double termCount, double collectionTermCount, int documentFrequency) = 0;
    virtual bool isProbability() const = 0; // does it output probability?
};


// DirichletDocumentCentricTermScoreFunction
//
class DirichletDocumentCentricTermScoreFunction: public DocumentCentricTermScoreFunction {
    double _mu;
    int _documentLength; // DL
    int _collectionLength; // N
    double _documentLengthPlusMu;

public:
    DirichletDocumentCentricTermScoreFunction(double mu, int documentLength, int collectionLength): 
	_mu(mu), _documentLength(documentLength), _collectionLength(collectionLength) 
    { 
	_documentLengthPlusMu = _documentLength + _mu;
    }

    double compute(double termCount, double collectionTermCount, int documentFrequency) {
	return (termCount + _mu * collectionTermCount / _collectionLength) / _documentLengthPlusMu;
    }

    double computeLogProb(double termCount, double collectionTermCount, int documentFrequency) {
	return log(compute(termCount, collectionTermCount, documentFrequency));
    }

    bool isProbability() const { return true; }
};


// JelinekMercerDocumentCentricTermScoreFunction
//
class JelinekMercerDocumentCentricTermScoreFunction: public DocumentCentricTermScoreFunction {
    double _collectionLambda;
    double _foregroundLambda;
    int _documentLength; // DL
    int _collectionLength; // N

public:
    JelinekMercerDocumentCentricTermScoreFunction(
	double collectionLambda, int documentLength, int collectionLength):
	_collectionLambda(collectionLambda), _foregroundLambda(1 - collectionLambda),
	_documentLength(documentLength), _collectionLength(collectionLength) {}

    double compute(double termCount, double collectionTermCount, int documentFrequency) {
	return _foregroundLambda * termCount / _documentLength + 
	    _collectionLambda * collectionTermCount / _collectionLength;
    }

    double computeLogProb(double termCount, double collectionTermCount, int documentFrequency) {
	return log(compute(termCount, collectionTermCount, documentFrequency));
    }

    bool isProbability() const { return true; }
};


// TFIDFDocumentCentricTermScoreFunction
//
class TFIDFDocumentCentricTermScoreFunction: public DocumentCentricTermScoreFunction {
    int _qTF;
    double _k1, _b, _k3;
    bool _okapi;
    double _documentLength, _collectionLength;
    int _documentCount;

    double _averageDocumentLength;
    double _k1TimesSomething;
    double _weight;
      
public:
    TFIDFDocumentCentricTermScoreFunction(
	int qTF, double k1, double b, double k3, bool okapi, double documentLength, 
	double collectionLength, int documentCount): 
	_qTF(qTF), _k1(k1), _b(b), _k3(k3), _okapi(okapi), _documentLength(documentLength), 
	_collectionLength(collectionLength), _documentCount(documentCount) 
    {
	_averageDocumentLength = _collectionLength / _documentCount;
	_k1TimesSomething = _k1 * ((1 - _b) + _b * _documentLength / _averageDocumentLength);
	_weight = _okapi? ((_k3 + 1) / (_k3 + _qTF)): (1000.0 / (1000.0 + _qTF));
    }

    double compute(double termCount, double collectionTermCount, int documentFrequency) {
	double idf = _okapi? 
	    log((_collectionLength - documentFrequency + 0.5) / (documentFrequency + 0.5)):
	    log((_collectionLength + 1.0) / (documentFrequency + 0.5));

	if (_okapi)
	    return _weight * idf * (_k1 + 1) * termCount / (termCount + _k1TimesSomething);
	else
	    return _weight * idf * idf * _k1 * termCount / (termCount + _k1TimesSomething);
    }

    double computeLogProb(double termCount, double collectionTermCount, int documentFrequency) {
	throw 0;
    }

    bool isProbability() const { return false; }
};

// KLDivergenceDocumentCentricTermScoreFunction
//
class KLDivergenceDocumentCentricTermScoreFunction: public DocumentCentricTermScoreFunction {
    double _collectionLambda;
    double _foregroundLambda;
    int _documentLength; // DL
    int _collectionLength; // N

    double _foregroundLambdaOverDocumentLength;
    double _collectionLambdaOverCollectionLength;
    double _foregroundLambdaTimesRatio;

public:
    KLDivergenceDocumentCentricTermScoreFunction(
	double collectionLambda, int documentLength, int collectionLength):
	_collectionLambda(collectionLambda), _foregroundLambda(1 - collectionLambda),
	_documentLength(documentLength), _collectionLength(collectionLength) 
    {
	_foregroundLambdaOverDocumentLength = _foregroundLambda / _documentLength;
	_collectionLambdaOverCollectionLength = _collectionLambda / _collectionLength;
	_foregroundLambdaTimesRatio = _foregroundLambda * (double(_collectionLength) / _documentLength);
    }

    double compute(double termCount, double collectionTermCount, int documentFrequency) {
	return 
	    (_foregroundLambdaOverDocumentLength * termCount +
		_collectionLambdaOverCollectionLength * collectionTermCount) * 
	    std::log(_foregroundLambdaTimesRatio * 
		     (collectionTermCount / termCount) + _collectionLambda);
    }

    double computeLogProb(double termCount, double collectionTermCount, int documentFrequency) {
	throw 0;
    }

    bool isProbability() const { return false; }
};


}

#endif
