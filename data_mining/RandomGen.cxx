//
// Created by Filippo Valle on 2019-02-04.
//

#include "RandomGen.h"

RandomGen* RandomGen::fgRandomGen = nullptr;

RandomGen RandomGen::Instance(uint64_t seed) {
    // if(!fgRandomGen) printf("NO existing rng!!!!!!!!!!");
    // else printf("already existing rng");
    if(!fgRandomGen) fgRandomGen = new RandomGen(seed);
    return *fgRandomGen;
}

void RandomGen::Destroy() {
    delete fgRandomGen;
}
