//
// Created by Filippo Valle on 2019-05-22.
//

#include "LabelsReshuffler.h"

void LabelsReshuffler::Shuffle(){
  fstream inputFile("files.dat", ios::in);

  if(!inputFile.is_open()){
    cerr<<"error opening files.dat"<<endl;
    return;
  }

  string line;
  getline(inputFile, line).good();//header
  auto header = tokenize(line);

  auto nClasses=header.size()-1;
  if(nClasses<1){
    cerr<<"not enough classes"<<endl;
    return;
  }else{
    cout<<"Reading "<<nClasses<<" classes"<<endl;
    BOOST_FOREACH(auto var , header) cout<<var<<"\t";
    cout<<endl;
  }

  vector<vector<string>> allFiles;
  allFiles.reserve(nClasses+1);
  for(int i=0; i<=nClasses;i++) allFiles.emplace_back(vector<string>());


  while (getline(inputFile, line).good()) {
    auto tokenizedLine = tokenize(line);
    int column = 0;
    BOOST_FOREACH(auto var, tokenizedLine) allFiles[column++].emplace_back(var);
  }

  inputFile.close();

  cout<<"shuffling.."<<endl;
  auto rng = RandomGen::Instance();
  static boost::random::uniform_int_distribution<uint64_t> distribution(0, allFiles[0].size()-1);
  for(uint8_t c=1; c<=nClasses; c++){
    cout<<"shuffling "<<header[c]<<endl;
    shuffleVector(allFiles[c], rng, distribution);
  }


  cout<<"saving output.."<<endl;
  fstream outputFile("files_shuffles.dat",ios::out);
  for(uint8_t c = 1; c<=nClasses;c++){
    outputFile<<","<<header[c];
  }
  outputFile<<endl;

  for(uint64_t iData=0; iData<allFiles[0].size();iData++){
    outputFile<<allFiles[0][iData];
    for(uint8_t c = 1; c<=nClasses;c++){
      outputFile<<","<<allFiles[c][iData];
    }
    outputFile<<endl;
  }

  outputFile.close();

}

void LabelsReshuffler::shuffleVector(std::vector<string> &vec, RandomGen &rng,
  boost::random::uniform_int_distribution<uint64_t> &distr, int nstat) {
    for(int stat=0; stat<nstat;stat++)
    for(uint64_t i =0; i< vec.size();i++)
    iter_swap(vec.begin()+i, vec.begin()+distr(rng));
  }

  std::vector<std::string> LabelsReshuffler::tokenize(const string& line){
    // escape char is \ , fields are seperated by , , some fields may be quoted with
    boost::escaped_list_separator<char> sep( '\\', ',', '\0' ) ;
    boost::tokenizer< boost::escaped_list_separator<char> > tokenizer( line, sep ) ;
    return std::vector<std::string>(tokenizer.begin(), tokenizer.end()) ;
  }
