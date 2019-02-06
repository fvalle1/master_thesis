//
// Created by Filippo Valle on 2019-02-04.
//

#include "RandomGen.h"

RandomGen* RandomGen::fgRandomGen = nullptr;

RandomGen RandomGen::Instance(uint64_t seed) {
    if(!fgRandomGen) fgRandomGen = new RandomGen(seed);
    return *fgRandomGen;
}
