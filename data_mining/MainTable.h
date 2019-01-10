//
// Created by Filippo Valle on 2019-01-08.
//

#ifndef THESIS_DATA_MINING_MAINTABLE_H
#define THESIS_DATA_MINING_MAINTABLE_H


#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <boost/tokenizer.hpp>
#include <boost/random/random_device.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <pcg_random.hpp>

using namespace std;


typedef bool correlated;
typedef bool* matrix;


class MainTable {
public:
    MainTable();

    ~MainTable();

    void read(bool saveAbundancesOccurrences = false);

    void readBinary();

    void SaveBinary();

    void ExtimateCorrelations();

    uint64_t get(uint64_t component, uint64_t realization) {
        if (component < fNComponents && realization < fNRealizations)
            return fData[fNRealizations * component + realization] ? 1 : 0;
        else return 0;
    };

private:
    uint64_t fNComponents;
    uint64_t fNRealizations;
    matrix fData;

    std::vector<std::string> tokenize(const std::string &);


    void SaveTotalArray(const char *filename, uint64_t length, double *X);
    void SaveHeapData(double *VocabularySize);
    double GetEntropy(double *X, uint64_t l);
    void ExtimateHXY();

    static const double constexpr fThreshold = 0.;

};

#endif //THESIS_DATA_MINING_MAINTABLE_H
