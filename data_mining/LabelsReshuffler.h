//
// Created by Filippo Valle on 2019-05-22.
//

#ifndef THESIS_DATA_MINING_LABELSRESHUFFLER_H
#define THESIS_DATA_MINING_LABELSRESHUFFLER_H


#include <fstream>
#include <sstream>
#include <vector>
#include <iostream>
#include <algorithm>

#include <omp.h>

#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>
#include <boost/random/uniform_int_distribution.hpp>

#include "RandomGen.h"

using namespace std;

class LabelsReshuffler {
public:
    static void Shuffle();

    static std::vector<std::string> tokenize(const string&);

    static void shuffleVector(std::vector<string> &vec, RandomGen &rng,
                              boost::random::uniform_int_distribution<uint64_t> &distr);
};


#endif //THESIS_DATA_MINING_LABELSRESHUFFLER_H
