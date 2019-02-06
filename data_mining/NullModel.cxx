//
// Created by Filippo Valle on 2019-01-16.
//

#include "NullModel.h"

void NullModel::GenerateNullData() {
    cout << "Generating null model" << endl;
    fstream occ("A.dat", ios_base::in);
    fstream voc("vocabulary_size.dat", std::ios::in);

    if(occ.is_open()&&voc.is_open()) {
        long double o;
        std::vector<long double> probabilities;
        probabilities.reserve(fNComponents);
        while(occ>>o){
            probabilities.push_back(o);
        }
        occ.close();
        cout<<"loaded occurences.."<<endl;

        // Make a random number engine
        auto rng = RandomGen::Instance();

        // Choose a random multinomial
        boost::random::discrete_distribution<uint16_t> distr(probabilities);

        probabilities.clear(); //a copy is stored in boost::ranodm::dicrete_distribution
        delete[] fData; //just to free some RAM I'll reallocate in the future

        auto nullData = new uint16_t[fNRealizations*fNComponents];

        double M;
        uint64_t cRealization = 0;
        auto counts = new uint16_t[fNComponents];
        for(uint64_t i = 0; i < fNComponents; i++) counts[i]=0;

        cout<<"generating data.."<<endl;
        //for realization
        while (voc >> M) {
            for(uint64_t word = 0; word < M; word++) counts[distr(rng)]++;
            printf("\r%llu/%llu",cRealization,fNRealizations);

            for(uint64_t component = 0; component < fNComponents; component++)
                nullData[fNRealizations*component + cRealization] = counts[component];


            for(uint64_t i = 0; i < fNComponents; i++) counts[i]=0;
            cRealization++;
        }
        delete[] counts;
        voc.close();

        cout<<endl<<"saving file.."<<endl;
        fstream file("nullTable.csv", std::ios::out);
        file<<endl; //mime header line
        for(uint64_t component = 0; component < fNComponents; component++) {
          file<<",";
            printf("\r%llu/%llu", component, fNComponents);
            for (uint64_t realization = 0; realization < cRealization; realization++) {
                file << nullData[fNRealizations * component + realization];
                if(realization < cRealization - 1) file << ",";
            }
            file<<"\n";
        }

        delete[] nullData;

        file.close();

        cout << endl;

        fData = new correlated[fNRealizations*fNComponents];
    }else{
        cerr<<"Error reading files"<<endl;
    }
}

void NullModel::GenerateNullBinary() const {
    cout << "Generating null model binary" << endl;
    fstream occ("O.dat", ios_base::in);

    if (!occ.is_open()) {
        cerr << "File not found" << endl;
    } else {

        double o;
        uint64_t idata = 0;

        // Make a random number engine
        auto rng = RandomGen::Instance();

        // Choose a random mean between 0 and 1
        boost::random::uniform_01<double> uniform_dist;

        while (occ >> o) {
            for (int realization = 0; realization < fNRealizations; ++realization) {
                if (idata % 100000 == 0) printf("\r%llu/%llu", idata, fNRealizations * fNComponents);
                bool point = uniform_dist(rng) < o;
                fData[idata++] = point;
            }
        }

        occ.close();
        cout << endl;
    }
}
