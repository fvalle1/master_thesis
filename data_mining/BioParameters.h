//
// Created by Filippo Valle on 2019-04-12.
//

#ifndef THESIS_DATA_MINING_BIOPARAMETERS_H
#define THESIS_DATA_MINING_BIOPARAMETERS_H


class BioParameters {
public:
    static constexpr int getENSLenght(){
        //return 15;//Home sapiens
        return 18; //mouse
    }

    static constexpr int getSampleIdLenght(){
        //return 36;//TCGA
        return 18; //mouse
    }
};


#endif //THESIS_DATA_MINING_BIOPARAMETERS_H
