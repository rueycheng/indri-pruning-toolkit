#ifndef FILTERED_INDEX_H_
#define FILTERED_INDEX_H_

#include <boost/dynamic_bitset.hpp>

#include "indri/Index.hpp"

namespace indri_contrib {

// FilteredDocListIterator
//
class FilteredDocListIterator: public indri::index::DocListIterator {
    indri::index::DocListIterator* _iter;
    const boost::dynamic_bitset<unsigned long>* _bitmap;
    size_t _pos;
    bool _ownsIterator;

public:
    FilteredDocListIterator(indri::index::DocListIterator* iter,
				 const boost::dynamic_bitset<unsigned long>* bitmap)
	: _iter(iter), _bitmap(bitmap), _pos(-1), _ownsIterator(true) {}

    FilteredDocListIterator() 
	:_iter(0), _bitmap(0), _pos(-1), _ownsIterator(false) { }
    
    ~FilteredDocListIterator() { 
	if (_ownsIterator) delete _iter; 
    }
    
    void reassign(indri::index::DocListIterator* iter, 
		  const boost::dynamic_bitset<unsigned long>* bitmap)
    {
	assert(!_ownsIterator);
	_iter = iter;
	_bitmap = bitmap;
	_pos = 0;

	// Update document frequencies and term frequencies (approximately)
	int filteredDocumentCount = bitmap->count();
	indri::index::TermData* termData = iter->termData();
	termData->corpus.totalCount -= filteredDocumentCount;
	termData->corpus.documentCount -= filteredDocumentCount;

	// TODO: Add a quick skip to the end if all bits are set
	if (_bitmap->test(_pos)) nextEntry();
    }
    
    void startIteration() {
	assert(_ownsIterator);
	assert(_iter);
	_iter->startIteration();
	_pos = 0;
    }
    
    indri::index::TermData* termData() {
	assert(_iter);
	return _iter->termData();
    }
    
    const indri::utility::greedy_vector<TopDocument>& topDocuments() {
	cerr << "No support\n";
	throw 0;
    }
    
    indri::index::DocListIterator::DocumentData* currentEntry() {
	assert(_iter);
	assert(!_bitmap->test(_pos));
	return _iter->currentEntry();
    }
    
    bool nextEntry() {
	assert(_iter);
	while (_iter->nextEntry() && _bitmap->test(++_pos)) ;
	return !finished();
    }
    
    bool nextEntry(lemur::api::DOCID_T documentID) {
	cerr << "No support\n";
	throw 0;
    }
    
    bool finished() {
	assert(_iter);
	return _iter->finished();
    }

    bool isFrequent() const {
	return _iter->isFrequent();
    }
};


// FilteredDocListFileIterator
//
template<typename Filter>
class FilteredDocListFileIterator: public indri::index::DocListFileIterator {
    indri::index::DocListFileIterator* _iter;
    const Filter& _filter;

    size_t _pos;
    FilteredDocListIterator _docListIterator;

    void _fixCurrentEntry() {
	indri::index::DocListFileIterator::DocListData* docListData = _iter->currentEntry();
	_docListIterator.reassign(docListData->iterator, &_filter[++_pos]);
	docListData->iterator = &_docListIterator;
    }

public:
    FilteredDocListFileIterator(indri::index::DocListFileIterator* iter, 
				     const Filter& filter)
	: _iter(iter), _filter(filter), _pos(-1) { }

    ~FilteredDocListFileIterator() {
	delete _iter;
    }

    bool finished() const {
	return _iter->finished();
    }

    void startIteration() {
	_iter->startIteration();
	if (!finished()) _fixCurrentEntry();
    }

    bool nextEntry() {
	if (_iter->nextEntry()) _fixCurrentEntry();
	return !finished();
    }

    indri::index::DocListFileIterator::DocListData* currentEntry() {
	return _iter->currentEntry();
    }

    const indri::index::DocListFileIterator::DocListData* currentEntry() const {
	return _iter->currentEntry();
    }
};


// FilteredIndex
//
template<typename Filter>
class FilteredIndex: public indri::index::Index {
    const Filter& _filter;
    indri::index::Index& _index;
    
public:
    FilteredIndex(indri::index::Index& index, const Filter& filter)
	: _index(index), _filter(filter) {}

    void close() { _index.close(); }

    lemur::api::DOCID_T documentBase() { return _index.documentBase(); }
    lemur::api::DOCID_T documentMaximum() { return _index.documentMaximum(); }
    lemur::api::TERMID_T term(const char* term) { return _index.term(term); }
    lemur::api::TERMID_T term(const std::string& term) { return _index.term(term); }
    std::string term(lemur::api::TERMID_T termID) { return _index.term(termID); }

    int field(const char* fieldName) { return _index.field(fieldName); }
    int field(const std::string& fieldName) { return _index.field(fieldName); }
    std::string field(int fieldID) { return _index.field(fieldID); }

    int documentLength(lemur::api::DOCID_T documentID) { return _index.documentLength(documentID); }
    UINT64 documentCount() { return _index.documentCount(); }
    UINT64 documentCount(const std::string& term) { return _index.documentCount(term); }

    UINT64 uniqueTermCount() { return _index.uniqueTermCount(); }

    UINT64 termCount(const std::string& term) { return _index.termCount(term); }
    UINT64 termCount() { return _index.termCount(); }

    UINT64 fieldTermCount(const std::string& field) { return _index.fieldTermCount(field); }
    UINT64 fieldTermCount(const std::string& field, const std::string& term) { return _index.fieldTermCount(field, term); }

    UINT64 fieldDocumentCount(const std::string& field) { return _index.fieldDocumentCount(field); }
    UINT64 fieldDocumentCount(const std::string& field, const std::string& term) { return _index.fieldDocumentCount(field, term); }
    
    indri::index::DocListIterator* docListIterator(lemur::api::TERMID_T termID) { 
	return _index.docListIterator(termID); 
    }

    indri::index::DocListIterator* docListIterator(const std::string& term) { 
	return _index.docListIterator(term); 
    }

    indri::index::DocListFileIterator* docListFileIterator() { 
	return new FilteredDocListFileIterator<Filter>(_index.docListFileIterator(), _filter); 
    }
                        
    indri::index::DocExtentListIterator* fieldListIterator(int fieldID) { return _index.fieldListIterator(fieldID); }
    indri::index::DocExtentListIterator* fieldListIterator(const std::string& field) { return _index.fieldListIterator(field); }
    const indri::index::TermList* termList(lemur::api::DOCID_T documentID) { return _index.termList(documentID); }
    indri::index::TermListFileIterator* termListFileIterator() { return _index.termListFileIterator(); }
    indri::index::DocumentDataIterator* documentDataIterator() { return _index.documentDataIterator(); }

    indri::index::VocabularyIterator* frequentVocabularyIterator() { return _index.frequentVocabularyIterator(); }
    indri::index::VocabularyIterator* infrequentVocabularyIterator() { return _index.infrequentVocabularyIterator(); }
    indri::index::VocabularyIterator* vocabularyIterator() { return _index.vocabularyIterator(); }

    indri::thread::Lockable* iteratorLock() { return _index.iteratorLock(); }
    indri::thread::Lockable* statisticsLock() { return _index.statisticsLock(); }
};

}

#endif
