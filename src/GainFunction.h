#ifndef GAIN_FUNCTION_H_
#define GAIN_FUNCTION_H_

#include <cmath>

#include "indri/Parameters.hpp"

// Gain functions
struct LogarithmGain {
    void loadParameters(indri::api::Parameters& parameters) {}
    double operator()(double p) { return std::log(p); }
};

struct IdentityGain {
    void loadParameters(indri::api::Parameters& parameters) {}
    double operator()(double p) { return p; }
};

struct SquareRootGain {
    void loadParameters(indri::api::Parameters& parameters) {}
    double operator()(double p) { return std::sqrt(p); }
};

struct MinusReciprocalGain {
    void loadParameters(indri::api::Parameters& parameters) {}
    double operator()(double p) { return -1.0 / p; }
};

struct MinusExponentialGain {
    double _oneMinusAlpha;

    MinusExponentialGain(): _oneMinusAlpha(-1.0) {}

    void loadParameters(indri::api::Parameters& parameters) { 
	double alpha = parameters.get("alpha", 3.0);
	_oneMinusAlpha = 1.0 - alpha;
    }

    double operator()(double p) { return -std::pow(p, _oneMinusAlpha); }
};

template<class Gain1, class Gain2> struct MixtureGain {
    Gain1 _gain1;
    Gain2 _gain2;

    double _lambda;
    double _oneMinusLambda;

    MixtureGain(): _lambda(0.5), _oneMinusLambda(0.5) {}

    void loadParameters(indri::api::Parameters& parameters) {
	// FIXME: sharing namespace is not safe
	_gain1.loadParameters(parameters);
	_gain2.loadParameters(parameters);
	_lambda = parameters.get("mixture.lambda", 0.5);
	_oneMinusLambda = 1.0 - _lambda;
    }

    double operator()(double p) { return _lambda * _gain1(p) + _oneMinusLambda * _gain2(p); }
};

#endif
