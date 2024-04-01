#ifndef HELPERS_HH
#define HELPERS_HH

#include <algorithm>
#include <cmath> 
#include <memory>
#include <random>
#include <vector>
#include <boost/format.hpp>
#include <boost/any.hpp>

#include "DataFormats/PatCandidates/interface/Jet.h"

class Helper {
    public:
        Helper();
        ~Helper();
        bool JetID(const pat::Jet &jet, std::string year);
};

#endif