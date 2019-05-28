//
// Created by Filippo Valle on 2019-02-04.
//

#ifndef THESIS_DATA_MINING_RANDOMGEN_H
#define THESIS_DATA_MINING_RANDOMGEN_H

#include <pcg_random.hpp>

class RandomGen : public pcg32_fast {
public:
    static RandomGen Instance(uint64_t seed = 42);
    static void Destroy();

private:
    explicit RandomGen(uint64_t seed) : pcg32_fast(seed){ }
    RandomGen() = default;

    static RandomGen* fgRandomGen;
};


#endif //THESIS_DATA_MINING_RANDOMGEN_H
