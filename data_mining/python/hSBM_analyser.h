//
// Created by Filippo Valle on 2019-05-28.
//

#include <boost/python.hpp>
#include <fstream>

#include "LabelsReshuffler.h"

namespace hSBM {

    class hsbm{
    public:
        inline void analyse_hsbm_results() {
            using namespace boost::python;


        }


        inline static void sample_kullback_liebler(){
            fstream topicdistr("topic-dist.csv",ios::in);
            fstream klout("topic-kl.csv", ios::out);
            string line, nextline;

            if(!topicdistr.is_open()){
                cerr<<"Error reading topic-dist.csv"<<endl;
                return;
            }

            string header;
            getline(topicdistr, header);

            vector<vector<string>> table;

            while(getline(topicdistr, line).good()){
                table.emplace_back(LabelsReshuffler::tokenize(line));
            }

            topicdistr.close();


            for(uint16_t first = 0; first<table.size(); first++) {
                auto firstline = table[first];

                for (uint16_t second = first + 1; second < table.size(); second++) {
                    printf("\r%d/%d", first, second);
                    double kl = 0.;

                    auto secondline = table[second];

                    for (uint16_t x=2; x<firstline.size();x++){
                        //cout<<x<<firstline[x]<<"\t\t"<<secondline[x]<<endl;
                        auto p = stod(firstline[x]);
                        auto q = stod(secondline[x]);
                        if (p!=0 && q!=0){
                            kl += p*log2(p/q);
                            kl += q*log2(q/p);
                        }

                    }
                    klout<<firstline[1]<<","<<secondline[1]<<","<<kl<<endl;
                }

            }
            cout<<endl;
            klout.close();
        }
    };

}
