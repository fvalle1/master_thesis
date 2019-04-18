//
// Created by Filippo Valle on 2019-01-16.
//

#ifndef THESIS_DATA_MINING_NULLMODEL_H
#define THESIS_DATA_MINING_NULLMODEL_H


#include "MainTable.h"
#include <boost/random/uniform_01.hpp>
#include <boost/random/discrete_distribution.hpp>
#include "pcg_random.hpp"

class SamplingModel : public MainTable{
public:
    SamplingModel():MainTable(){};
    ~SamplingModel() = default;

    void GenerateNullData(int maxStatistic = 10);
    void GenerateNullBinary() const;


};


#endif //THESIS_DATA_MINING_NULLMODEL_H
