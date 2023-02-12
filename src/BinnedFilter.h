#ifndef BINNED_FILTER_H_
#define BINNED_FILTER_H_

#include <vector>

#include <boost/dynamic_bitset.hpp>
#include <boost/foreach.hpp>

namespace indri_contrib {

class BinnedFilter {
    std::vector<boost::dynamic_bitset<unsigned long>> _bitmaps;

public:
    void clear() {
	_bitmaps.clear();
    }

    size_t size() const {
	return _bitmaps.size();
    }

    void add(const boost::dynamic_bitset<unsigned long>& bitmap) {
	_bitmaps.push_back(bitmap);
    }

    double ratio() const {
	size_t totalCount = 0, totalSize = 0;

	BOOST_FOREACH (const boost::dynamic_bitset<unsigned long>& bitmap, _bitmaps) {
	    totalCount += bitmap.count();
	    totalSize += bitmap.size();
	}

	return double(totalCount) / totalSize;
    }

    const boost::dynamic_bitset<unsigned long>& operator[](size_t pos) const {
	return _bitmaps[pos];
    }

    boost::dynamic_bitset<unsigned long>& operator[](size_t pos) {
	return _bitmaps[pos];
    }
};

}

#endif
