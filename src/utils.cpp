#include "IndriPruneIndex.h"

#include <boost/foreach.hpp>

void warn(const char* message) {
    std::cerr << message << "\n";
}

void warn(boost::format& fmt) {
    std::cerr << fmt << "\n";
}

void die(const char* message) {
    warn(message);
    exit(1);
}

void die(boost::format& fmt) {
    warn(fmt);
    exit(1);
}

void requires(indri::api::Parameters& parameters, 
	      const std::vector<std::string>& keys) {
    bool found = false;

    BOOST_FOREACH (const std::string& key, keys) {
	if (!parameters.exists(key)) {
	    found = true;
	    warn(boost::format("Parameter '%s' is required") % key);
	}
    }

    if (found) exit(1);
}

void requires(indri::api::Parameters& parameters, const char* key) {
    std::vector<std::string> keys = { key };
    requires(parameters, keys);
}

void parseSpec(indri::api::Parameters& converted, const std::string& spec) {
    int nextComma = 0, nextColon = 0;

    for (int location = 0; location < spec.length();) {
	nextComma = spec.find(',', location);
	nextColon = spec.find(':', location);

	std::string key = spec.substr(location, nextColon-location);
	std::string value = spec.substr(nextColon+1, nextComma-nextColon-1);
	converted.set(key, value);

	if (nextComma > 0) location = nextComma+1;
	else location = spec.size();
    }
}

