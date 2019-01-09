//
// Created by Filippo Valle on 2019-01-08.
//

#include "MainTable.h"


MainTable::MainTable() {
    //init size of the table
    fNComponents = 60483;
    fNRealizations = 11571;

    fData = new correlated[fNRealizations*fNComponents];
}


MainTable::~MainTable() {
    delete[] fData;
    printf("\n\n");
}


void MainTable::read(bool saveAbundancesOccurrences) {
    cout<<"Reading maintable.csv"<<endl;
    fstream file("mainTable.csv", std::ios::in);
    string line;
    uint64_t idata = 0;

    double A[fNComponents];
    double O[fNComponents];
    double VocabularySize[fNRealizations];
    if(saveAbundancesOccurrences) {
        for (uint64_t i = 0; i < fNComponents; i++) A[i] = 0.;
        for (uint64_t i = 0; i < fNComponents; i++) O[i] = 0.;
        for (uint64_t i = 0; i < fNRealizations; i++) VocabularySize[i] = 0.;
    }

    //read header line
    getline(file, line).good();


    bool firstRead = false;
    while (getline(file, line).good()) {
        auto tokenizedLine = tokenize(line);

        if(!firstRead) {
            fNRealizations=tokenizedLine.size() - 1;
            firstRead = true;
        }
        if(idata%10000==0) printf("\r%llu/%llu",idata,fNRealizations*fNComponents);


        //+1 because first column is component-id
        for (auto token = tokenizedLine.begin() + 1; token != tokenizedLine.end(); token++) {
            //cout<<*token<<endl;
            double value = std::stod(*token);
            bool binaryValue = value > MainTable::fThreshold;
            fData[idata++] = binaryValue;

            if(saveAbundancesOccurrences) {
                A[idata / fNRealizations] += value;
                O[idata / fNRealizations] += (binaryValue ? 1. : 0.) / fNRealizations;
            }
            VocabularySize[idata%fNRealizations] += value;
        }
    }

    file.close();
    cout<<endl;

    if(saveAbundancesOccurrences) {
        SaveTotalArray("A.dat", fNComponents, A);
        SaveTotalArray("O.dat", fNComponents, O);
    }
    SaveTotalArray("vocabulary_size.dat", fNRealizations, VocabularySize);
    SaveHeapData(VocabularySize);
}


void MainTable::readBinary() {
    cout<<"Reading maintable.csv"<<endl;
    fstream file("binaryTable.csv", std::ios::out);
    string line;
    uint64_t idata = 0;

    bool firstRead = false;
    while (getline(file, line).good()) {
        auto tokenizedLine = tokenize(line);

        if(!firstRead) {
            fNRealizations=tokenizedLine.size() - 1;
            firstRead = true;
        }
        if(idata%10000==0) printf("\r%llu/%llu",idata,fNRealizations*fNComponents);

        //+1 because first column is component-id
        for (auto token = tokenizedLine.begin() + 1; token != tokenizedLine.end(); token++) {
            //cout<<*token<<endl;
            bool value = std::stoi(*token)==1;
            fData[idata++] = value;
        }

    }

    file.close();
    cout<<endl;
}

std::vector<std::string> MainTable::tokenize( const std::string& line )
{
    // escape char is \ , fields are seperated by , , some fields may be quoted with
    boost::escaped_list_separator<char> sep( '\\', ',', '\0' ) ;
    boost::tokenizer< boost::escaped_list_separator<char> > tokenizer( line, sep ) ;
    return std::vector<std::string>( tokenizer.begin(), tokenizer.end() ) ;
}

void MainTable::SaveBinary() {
    cout<<"Saving binary matrix"<<endl;
    fstream file("binaryTable.csv", std::ios::out);
    for(uint64_t idata = 0; idata < fNComponents*fNRealizations; idata++){
        file<<(fData[idata]?1:0);
        if((idata+1)%fNRealizations==0) file<<"\n";
        else file <<",";
    }
    file.close();
}

void MainTable::SaveTotalArray(const char *filename, uint64_t length, double *X) {
    printf("Saving total %s\n",filename);

    fstream file(filename, std::ios::out);
    file<<X[0];
    for(uint64_t i=1; i<length; i++){
        file<<"\n"<<X[i];
    }

    file.close();
}

void MainTable::SaveHeapData(double *VocabularySize) {
    cout << "Saving Heap data" << endl;
    fstream file("heap.dat", std::ios::out);
    uint64_t vocabulary[fNComponents];
    for (uint64_t i = 0; i < fNComponents; i++) vocabulary[i] = 0;

    uint64_t cNumberOfDifferentWords = 0;
    long double cVocabularySize = 0.;

    pcg_extras::seed_seq_from<boost::random::random_device> seed_seq;
    pcg32_fast rng(seed_seq);
    boost::random::uniform_int_distribution<uint64_t> distribution(0, fNRealizations);

    for (uint64_t r = 0; r < fNRealizations; r++) {
        printf("\r%llu/%llu",r,fNRealizations);

        uint64_t realization = distribution(rng);
        cVocabularySize += VocabularySize[realization];
        for (uint64_t component = 0; component < fNComponents; component++) {
            if (get(component, realization) != 0) vocabulary[component]++;
        }

        cNumberOfDifferentWords = 0;
        for (uint64_t i = 0; i < fNComponents; i++) cNumberOfDifferentWords += vocabulary[i] >= 1 ? 1 : 0;

        file<<cVocabularySize<<","<<cNumberOfDifferentWords<<endl;
        file.flush();
    }


    file.close();
    cout<<endl;
}