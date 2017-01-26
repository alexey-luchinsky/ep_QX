//
// Created by Poslavsky on 24/03/16.
//

#ifndef GGGPSIPSI_RANDOM_H
#define GGGPSIPSI_RANDOM_H

#include <stdexcept>

#include "algebra.h"
#include <boost/random.hpp>
#include <boost/random/mersenne_twister.hpp>


class Random {
private:
    long seed;
    boost::mt11213b inner_random;
    boost::uniform_real<dbl_type> flat;

public:
    Random() { }

    Random(long _seed_) {
        seed = _seed_;
        inner_random.seed(seed);
    }

    unsigned int get_seed() { return seed; };

    dbl_type rand(dbl_type min, dbl_type max) {
        if (min > max)
            throw runtime_error("max < min");
        return min + (max - min) * flat(inner_random);
    }
};


#endif //GGGPSIPSI_RANDOM_H
