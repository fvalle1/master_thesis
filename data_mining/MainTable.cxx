//
// Created by Filippo Valle on 2019-01-08.
//

#include "MainTable.h"

const std::pair<double, double> MainTable::fThreshold = {0.1, 1e5}; //as stated by Matteo Cereda


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


void MainTable::read(const char *tableFilename, bool saveAbundancesOccurrences, bool saveTitles) {
    cout << "Reading " << tableFilename << endl;
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

        //read header line and generate titles.txt
        //NB begins()+1 because first row is just "gene" word
        std::string header;
        if(getline(file, header).good()){
            if(saveTitles) {
                fstream titles("titles.txt", ios::out);
                auto names = tokenize(header);
                std::for_each(names.begin() + 1, names.end(), [&](string name) {
                    titles << name << "\n";
                });
                titles.close();
            }
        }


        bool firstRead = false;
        uint64_t actualWords = 0;
        while (getline(file, line).good()) {
            auto tokenizedLine = tokenize(line);

            if (!firstRead) {
                fNRealizations = tokenizedLine.size() - 1;
                firstRead = true;
            }
            if (idata % 5000 == 0) printf("\r%llu/%llu", idata, fNRealizations * fNComponents);

            //+1 because first column is component-id
            for (auto token = tokenizedLine.begin() + 1; token != tokenizedLine.end(); token++) {
//                cout<<*token<<endl;
                double value = std::stod(*token);
                bool binaryValue = (value > MainTable::fThreshold.first) && (value < MainTable::fThreshold.second);
                fData[idata++] = binaryValue;

                if (saveAbundancesOccurrences) {
                    A[idata / fNRealizations] += (binaryValue ? value : 0.); //underthreshold shoud not be simulated
                    O[idata / fNRealizations] += (binaryValue ? 1. : 0.) / fNRealizations;
                }
                VocabularySize[idata % fNRealizations] += value;
            }
            actualWords++;
        }

        fNComponents = actualWords;
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
                if (saveAbundancesOccurrences) {
                    SaveTotalArray("vocabulary_size.dat", VocabularySize, fNRealizations);
                    SaveHeapsData(VocabularySize);
                }
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
        printf("\r%llu/%llu", realization+1, fNRealizations);

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

#pragma omp parallel for shared(file)
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

#pragma omp critical
            {
                h = GetEntropy(4, P);
                hx[0] = P[0] + P[1]; //prob of Zeros in first
                hx[1] = 1. - hx[0]; //prob of Ones in first

                hy[0] = P[0] + P[2]; //prob of Zeros in second
                hy[1] = 1. - hy[0]; //prob of Ones in second


                H = GetEntropy(2, hx, firstComponent) + GetEntropy(2, hy, secondComponent) - h;
//            if(abs(H)>1) cerr<<endl<<"{"<<hx[0]<<","<<hx[1]<<"}"<<" {"<<hy[0]<<","<<hy[1]<<"}"<<"\t["<<P[0]<<"\t"<<P[1]<<"\t"<<P[2]<<"\t"<<P[3]<<"]\t"<<GetEntropy(2, hx, firstComponent)<<"\t"<<GetEntropy(2, hy, secondComponent)<<"\t"<<h<<"\t"<<H<<endl;
//components commented to avoid unuseful disk space
                file /*<< firstComponent << "\t" << secondComponent << "\t" */<< H << endl;
            }
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


void MainTable::SaveMeansVariances(const char *filename) {
    cout << "Reading " << filename << endl;
    fstream file(filename, std::ios::in);
    fstream meanVariances("meanVariances.csv", std::ios::out);
    meanVariances<<",mean,variance,type_of_gene"<<endl;
    if (!file.is_open()) {
        cerr << "file not found" << endl;
    } else {
        string line;

        //read header line and generate
        std::string header;
        getline(file, header).good();

        uint64_t ngenes = 0;
        while (getline(file, line).good()) {
            auto tokenizedLine = tokenize(line);

            printf("\rngenes: %llu", ++ngenes);

            // first column is ENSG-id
            auto gene = (*(tokenizedLine.begin())).substr(0,15);
            long double sum = 0.;
            long double sumsquare = 0.;
            uint64_t n = 0;

            for (auto token = tokenizedLine.begin() + 1; token != tokenizedLine.end(); token++) {
//                cout<<*token<<endl;
                double value = std::stod(*token);
                if(value > fThreshold.first && value < fThreshold.first) {
                    //fast algo https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Naive_algorithm
                    sum += value;
                    sumsquare += value * value;
                    n += 1;
                }
            }

            if(n>2){
                long double average = sum / n;
                long double variance = static_cast<long double>(sumsquare - (sum*sum)/n)/(n);// /n is on finite sample /(n-1) if infinite sample
                meanVariances<<gene<<","<<average<<","<<variance<<","<<" "<<endl;
            }
        }
    }
    cout<<endl;
}


void MainTable::MakeGraph() {
    cout<<"Making graph.xml"<<endl;
    fstream file("graph.xml", std::ios::out);

    ptree xmlstructure;
    ptree graphml;
    graphml.put("<xmlattr>.xmlns", "http://graphml.graphdrawing.org/xmlns");
    graphml.put("<xmlattr>.xmlns:xsi","http://www.w3.org/2001/XMLSchema-instance");
    graphml.put("<xmlattr>.xsi:schemaLocation","http://graphml.graphdrawing.org/xmlns http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd");

    ptree key0, key1, key2;
    key0.put("<xmlattr>.id", "key0");
    key0.put("<xmlattr>.for", "node");
    key0.put(boost::property_tree::ptree::path_type("<xmlattr>|attr.name",'|'), "count");
    key0.put(boost::property_tree::ptree::path_type("<xmlattr>|attr.type",'|'), "int");
    graphml.add_child("key", key0);

    key1.put("<xmlattr>.id", "key1");
    key1.put("<xmlattr>.for", "node");
    key1.put(boost::property_tree::ptree::path_type("<xmlattr>|attr.name",'|'), "kind");
    key1.put(boost::property_tree::ptree::path_type("<xmlattr>|attr.type",'|'), "int");
    graphml.add_child("key", key1);

    key2.put("<xmlattr>.id", "key2");
    key2.put("<xmlattr>.for", "node");
    key2.put(boost::property_tree::ptree::path_type("<xmlattr>|attr.name",'|'), "name");
    key2.put(boost::property_tree::ptree::path_type("<xmlattr>|attr.type",'|'), "string");
    graphml.add_child("key", key2);

    ptree graph;
    graph.put("<xmlattr>.id", "G");
    graph.put("<xmlattr>.edgedefault", "undirected");
    graph.put(boost::property_tree::ptree::path_type("<xmlattr>|parse.nodeids",'|'),"canonical");
    graph.put(boost::property_tree::ptree::path_type("<xmlattr>|parse.edgeids",'|'),"canonical");
    graph.put(boost::property_tree::ptree::path_type("<xmlattr>|parse.order",'|'),"nodesfirst");


    ptree worddata;
    worddata.put("<xmlattr>.key", "key1");
    worddata.put("", "0");

    uint64_t nodescount = 0;

    for(int w = 0; w < 10; w++){
        ptree node;
        std::ostringstream nodeid, nodename;
        nodeid << "n" << nodescount++;
        nodename<<"pluto";

        node.put("<xmlattr>.id", nodeid.str());
        node.add_child("data", worddata);

        ptree nodedata;
        nodedata.put("", nodename.str());
        nodedata.put("<xmlattr>.key", "key2");
        node.add_child("data", nodedata);

        graph.add_child("node", node);

    }

    ptree docdata;
    docdata.put("<xmlattr>.key", "key1");
    docdata.put("", "1");

    for(int d = 0; d < 10; d++) {
        std::ostringstream nodeid, nodename;
        ptree node;
        nodeid << "n" << nodescount++;
        nodename<<"pluto_doc";

        node.put("<xmlattr>.id", nodeid.str());
        node.add_child("data", docdata);

        ptree nodedata;
        nodedata.put("", nodename.str());
        nodedata.put("<xmlattr>.key", "key2");
        node.add_child("data", nodedata);

        graph.add_child("node", node);
    }

    ptree edgedata;
    edgedata.put("<xmlattr>.key", "key0");
    edgedata.put("", "1");
    for(int e = 0; e < 10; e++){
        std::ostringstream edgeid, startid, targetid;
        startid<<"n"<<4;
        targetid<<"n"<<5;

        edgeid<<"e"<<e;
        ptree edge;
        edge.put("<xmlattr>.id", edgeid.str());
        edge.put("<xmlattr>.source", startid.str());
        edge.put("<xmlattr>.target", targetid.str());
        edge.add_child("data", edgedata);
        graph.add_child("edge", edge);
    }

    graphml.add_child("graph", graph);
    xmlstructure.add_child("graphml", graphml);

    write_xml(file, xmlstructure, boost::property_tree::xml_writer_make_settings<std::string>(' ', 4));

}