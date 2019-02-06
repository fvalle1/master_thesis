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


void MainTable::read(const char *tableFilename, bool saveAbundancesOccurrences) {
    cout << "Reading maintable.csv" << endl;
    fstream file(tableFilename, std::ios::in);
    if (!file.is_open()) {
        cerr << "file not found" << endl;
    } else {
        string line;
        uint64_t idata = 0;

        double A[fNComponents];
        double O[fNComponents];
        double VocabularySize[fNRealizations];
        if (saveAbundancesOccurrences) {
            for (uint64_t i = 0; i < fNComponents; i++) A[i] = 0.;
            for (uint64_t i = 0; i < fNComponents; i++) O[i] = 0.;
            for (uint64_t i = 0; i < fNRealizations; i++) VocabularySize[i] = 0.;
        }

        //read header line
        getline(file, line).good();

        bool firstRead = false;
        while (getline(file, line).good()) {
            auto tokenizedLine = tokenize(line);

            if (!firstRead) {
                fNRealizations = tokenizedLine.size() - 1;
                firstRead = true;
            }
            if (idata % 10000 == 0) printf("\r%llu/%llu", idata, fNRealizations * fNComponents);


            //+1 because first column is component-id
            for (auto token = tokenizedLine.begin() + 1; token != tokenizedLine.end(); token++) {
                //cout<<*token<<endl;
                double value = std::stod(*token);
                bool binaryValue = value > MainTable::fThreshold;
                fData[idata++] = binaryValue;

                if (saveAbundancesOccurrences) {
                    A[idata / fNRealizations] += value;
                    O[idata / fNRealizations] += (binaryValue ? 1. : 0.) / fNRealizations;
                }
                VocabularySize[idata % fNRealizations] += value;
            }
        }

        file.close();
        cout << endl;

#pragma omp parallel sections
        {
#pragma omp section
            {
                if (saveAbundancesOccurrences) {
                    SaveTotalArray("A.dat", A, fNComponents);
                    SaveTotalArray("O.dat", O, fNComponents);
                }
                cout<<endl;
            }
#pragma omp section
            {
                SaveTotalArray("vocabulary_size.dat", VocabularySize, fNRealizations);
                SaveHeapsData(VocabularySize);
                cout<<endl;
            }
        }
    }
}

void MainTable::readBinary() {
    cout << "Reading maintable.csv" << endl;
    fstream file("binaryTable.csv", std::ios::in);

    if (!file.is_open()) {
        cerr << "file not found" << endl;
    } else {

        string line;
        uint64_t idata = 0;

        bool firstRead = false;
        while (getline(file, line).good()) {
            auto tokenizedLine = tokenize(line);

            if (!firstRead) {
                fNRealizations = tokenizedLine.size() - 1;
                firstRead = true;
            }
            if (idata % 1000 == 0) printf("\r%llu/%llu", idata, fNRealizations * fNComponents);

            //+1 because first column is component-id
            for (auto token = tokenizedLine.begin() + 1; token != tokenizedLine.end(); token++) {
                //cout<<*token<<endl;
                bool value = std::stoi(*token) == 1;
                fData[idata++] = value;
            }
        }

        file.close();
        cout << endl;
    }
}

std::vector<std::string> MainTable::tokenize( const std::string& line )
{
    // escape char is \ , fields are seperated by , , some fields may be quoted with
    boost::escaped_list_separator<char> sep( '\\', ',', '\0' ) ;
    boost::tokenizer< boost::escaped_list_separator<char> > tokenizer( line, sep ) ;
    return std::vector<std::string>(tokenizer.begin(), tokenizer.end()) ;
}

void MainTable::SaveBinary(const char *filename) {
    cout<<"Saving binary matrix"<<endl;
    fstream file(filename, std::ios::out);
    for(uint64_t idata = 0; idata < fNComponents*fNRealizations; idata++){
        if(idata%100000==0) printf("\r%llu/%llu",idata,fNRealizations*fNComponents);
        file<<(fData[idata]?1:0);
        if((idata+1)%fNRealizations==0) file<<"\n";
        else file <<",";
    }
    file.close();
    cout<<endl;
}

void MainTable::SaveTotalArray(const char *filename, const double *X, uint64_t length) {
    printf("Saving total %s\n",filename);

    fstream file(filename, std::ios::out);
    file<<X[0];
    for(uint64_t i=1; i<length; i++){
        file<<"\n"<<X[i];
    }

    file.close();
}

void MainTable::SaveHeapsData(const double *VocabularySize) {
    cout << "Saving Heaps data" << endl;
    fstream file("heaps.dat", std::ios::out);

    for (uint64_t realization = 0; realization < fNRealizations; ++realization) {
        printf("\r%llu/%llu", realization, fNRealizations);

        uint64_t cNumberOfDifferentWords = 0;
        long double cVocabularySize = VocabularySize[realization];

#pragma omp parallel for reduction(+:cNumberOfDifferentWords)
        for (uint64_t component = 0; component < fNComponents; component++) {
            if (get(component, realization) != 0) cNumberOfDifferentWords++;
        }

        file << cVocabularySize << "," << cNumberOfDifferentWords << endl;
        file.flush();
    }

//#pragma omp master
    file.close();
    cout<<endl;
}


void MainTable::ExtimateCorrelations(const char *filename) {
    cout<<"Extimating correlations"<<endl;

    fNComponents = 1000;
    //fNRealizations = 100;
    ExtimateHXY(filename);
}

void MainTable::ExtimateHXY(const char *filename) {
    cout<<"Extimating H(X,Y)"<<endl;
    cout<<"words: "<<fNComponents<<"\t documents: "<<fNRealizations<<endl;

    fstream file(filename, ios_base::out);

    double norm = 1./fNRealizations;
    double H, h;
    double hx[2];
    double hy[2];

    for(uint64_t firstComponent = 0; firstComponent < fNComponents; firstComponent++){
        for(uint64_t secondComponent = firstComponent + 1; secondComponent < fNComponents; secondComponent++){
            printf("\r%llu/%llu",firstComponent,secondComponent);
            double P[4] = {0.};
            for (uint64_t realization = 0; realization < fNRealizations; ++realization) {
                auto x = get(firstComponent, realization);
                auto y = get(secondComponent, realization);
                if(x==y){
                    if(x==0){//00
                        P[0]+=norm;
                    }else{//11
                        P[3]+=norm;
                    }
                }else{
                    if(x==0){//01
                        P[1]+=norm;
                    }else{//10
                        P[2]+=norm;
                    }
                }
            }

            h = GetEntropy(4, P);
            hx[0] = P[0] + P[1]; //prob of Zeros in first
            hx[1] = 1. - hx[0]; //prob of Ones in first

            hy[0] = P[0] + P[2]; //prob of Zeros in second
            hy[1] = 1. - hy[0]; //prob of Ones in second


            H = GetEntropy(2, hx, firstComponent) + GetEntropy(2, hy, secondComponent) - h;
//            if(abs(H)>1) cerr<<endl<<"{"<<hx[0]<<","<<hx[1]<<"}"<<" {"<<hy[0]<<","<<hy[1]<<"}"<<"\t["<<P[0]<<"\t"<<P[1]<<"\t"<<P[2]<<"\t"<<P[3]<<"]\t"<<GetEntropy(2, hx, firstComponent)<<"\t"<<GetEntropy(2, hy, secondComponent)<<"\t"<<h<<"\t"<<H<<endl;

            file << firstComponent << "\t" << secondComponent << "\t" << H << endl;

        }

        file.flush();
    }

    file.close();
    cout<<endl;
}


double MainTable::GetEntropy(uint64_t l, double *X, const uint64_t component) {
    static std::map<uint64_t, double> cache;

    if(component<=fNComponents){
        auto it = cache.find(component);
        if (it != cache.end()){
            return it->second;
        }else{
            double H = SumEntropy(l, X);
            cache.insert(std::pair<uint64_t, double>(component, H));
            return H;
        }
    }else{
        return SumEntropy(l, X);
    }
}

double MainTable::SumEntropy(uint64_t l, double *X){
    double H = 0.;
    for(uint64_t i = 0; i < l; i++){
        double x = X[i];
        if(x>1e-5) H += x*log2(x);
    }
    return -H;
}
