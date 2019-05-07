//
// Created by Filippo Valle on 2019-04-12.
//

#ifndef THESIS_DATA_MINING_BIOPARAMETERS_H
#define THESIS_DATA_MINING_BIOPARAMETERS_H

enum species{
  kTCGA,
  kMouse,
  kGTEX
};

class BioParameters {
public:
  static constexpr int getENSLenght(){
    switch(fSpecie){
      case kTCGA:
      return 15; //Home sapiens
      case kMouse:
      return 18; //mouse
      case kGTEX:
      return 15;
    }
  }

  static constexpr int getSampleIdLenght(){
    switch(fSpecie){
      case kTCGA:
      return 36;//TCGA
      case kMouse:
      return 18; //mouse
      case kGTEX:
      return 23;
    }
  }

private:
  BioParameters() = default;
  ~BioParameters() = default;

  static const species fSpecie=kMouse;
};


#endif //THESIS_DATA_MINING_BIOPARAMETERS_H
