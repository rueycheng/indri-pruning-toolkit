#ifndef QUANTILE_SAMPLER_H_
#define QUANTILE_SAMPLER_H_

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/extended_p_square_quantile.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>

namespace indri_contrib {

// QuantileSampler
//
class QuantileSampler {
public:
    virtual void add(double v) = 0;
    virtual double query(double q) = 0;
    virtual std::vector<double> query(const std::vector<double>& qvec) = 0;

    void add(size_t n, const double v) {
	for (size_t i = 0; i < n; ++i) add(v);
    }

    template<typename Iterator> void add(Iterator first, Iterator last) {
	while (first != last) add(*first++);
    }
};


// FIXME: This would pollute the global namespace
using namespace boost::accumulators;


// PSquareQuantileSampler
//
class PSquareQuantileSampler: public QuantileSampler {
    std::vector<double> _probs;

    typedef accumulator_set<double, stats<tag::extended_p_square_quantile>> accumulator_t;
    accumulator_t _acc;

public:
    PSquareQuantileSampler(const std::vector<double>& probs): 
	 _acc(extended_p_square_probabilities = probs) { }

    void add(double v) { _acc(v); }

    double query(double q) {
	assert(q >= 0.0 && q <= 1.0);
	return quantile(_acc, quantile_probability = q);
    }

    std::vector<double> query(const std::vector<double>& qvec) {
	std::vector<double> result;
	for (size_t i = 0; i < qvec.size(); ++i)
	    result.push_back(query(qvec[i]));
	return result;
    }
};


// ReservoirQuantileSampler
//
class ReservoirQuantileSampler: public QuantileSampler {
    boost::mt19937 _rng;
    boost::uniform_01<> _prob;

    std::vector<double> _reservoir;
    size_t _count;

public:
    ReservoirQuantileSampler(size_t reservoirSize):
	_reservoir(reservoirSize), _count(0) {}

    void add(double v) {
	if (_count < _reservoir.size())
	    _reservoir[_count] = v;
	else {
	    size_t pos = _count *_prob(_rng);
	    if (pos < _reservoir.size())
		_reservoir[pos] = v;
	}

	++_count;
    }

    double _estimateQuantile(const std::vector<double>& x, double q) {
	assert(q >= 0.0 && q <= 1.0);

	size_t N = x.size();
	double h = (N + 1.0 / 3.0) * q + 1.0 / 3.0 - 1; // substract one because array index starts from 0
	size_t h_floor = std::floor(h);

	return x[h_floor] + (h - h_floor) * (x[h_floor + 1] - x[h_floor]);
    }

    double query(double q) {
	std::vector<double> x(_reservoir.begin(), _reservoir.end());
	std::stable_sort(x.begin(), x.end());

	return _estimateQuantile(x, q);
    }

    std::vector<double> query(const std::vector<double>& qvec) {
	std::vector<double> x(_reservoir.begin(), _reservoir.end());
	std::stable_sort(x.begin(), x.end());

	std::vector<double> result;
	for (size_t i = 0; i < qvec.size(); ++i)
	    result.push_back(_estimateQuantile(x, qvec[i]));
	return result;
    }
};

}

#endif
