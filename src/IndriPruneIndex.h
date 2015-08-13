#ifndef SIP_H_
#define SIP_H_

#include <vector>

#include <boost/format.hpp>

#include "indri/Index.hpp"
#include "indri/Parameters.hpp"

template<typename Key, typename Value> struct less_key {
    bool operator()(const pair<Key, Value>& lhs, const pair<Key, Value>& rhs) {
	return lhs.first < rhs.first;
    }
};

template<typename Key, typename Value> struct greater_key {
    bool operator()(const pair<Key, Value>& lhs, const pair<Key, Value>& rhs) {
	return lhs.first > rhs.first;
    }
};

template<typename Key, typename Value> struct less_value {
    bool operator()(const pair<Key, Value>& lhs, const pair<Key, Value>& rhs) {
	return lhs.second < rhs.second || lhs.second == rhs.second && lhs.first < rhs.first;
    }
};

template<typename Key, typename Value> struct greater_value {
    bool operator()(const pair<Key, Value>& lhs, const pair<Key, Value>& rhs) {
	return lhs.second > rhs.second || lhs.second == rhs.second && lhs.first < rhs.first;
    }
};


extern void warn(const char* message);
extern void warn(boost::format& fmt);
extern void die(const char* message);
extern void die(boost::format& fmt);

extern void requires(indri::api::Parameters& parameters, 
		     const std::vector<std::string>& keys);
extern void requires(indri::api::Parameters& parameters, 
		     const char* key);
extern void parseSpec(indri::api::Parameters& converted, 
		      const std::string& spec);

extern void copy_repository(indri::api::Parameters& parameters);
extern void check_index_health(indri::index::Index& index, 
			       bool dumpPostingList);
extern void print_direct_file(indri::api::Parameters& parameters);
extern void print_inverted_list(indri::api::Parameters& parameters);
extern void check_index_health(indri::api::Parameters& parameters);
extern void check_document_prior(indri::api::Parameters& parameters);

#endif
