#include "MainTable.h"
#include "TelegramWatch.h"
#include "NullModel.h"
#include <omp.h>
#include <iostream>


int main(int argc, const char** argv) // minimal test driver
{
    omp_set_nested(1);
    printf("Running TCGA\n");
#pragma omp parallel
    if(omp_get_thread_num()==0) printf("threads: %d\n", omp_get_num_threads());

    TelegramWatch watch("thesis");

    MainTable* TCGA;

    if(argc < 2){
        cerr<<"Please write some options"<<endl;
        cout<<"0 ---> read and extimate correlation mainTable.csv"<<endl;
        cout<<"1 ---> read mainTable.csv"<<endl;
        cout<<"2 ---> GenerateNullData"<<endl;
        cout<<"3 ---> read nullTable.csv"<<endl;
        cout<<"4 ---> read and extimate correlation nullTable.csv"<<endl;
        cout<<"5 ---> read and extimate means and variances"<<endl;
        cout<<"6 ---> read and makeCorpus"<<endl;
    }else {
        switch (std::atoi(argv[1])) {
            case 0:
                TCGA = new MainTable();
                TCGA->read("mainTable.csv", true);
                TCGA->ExtimateCorrelations();
                TCGA->~MainTable();
                break;
            case 1:
                TCGA = new MainTable();
                TCGA->read("mainTable.csv", true);
                TCGA->~MainTable();
                break;
            case 2:
                TCGA = new NullModel();
                ((NullModel *) (TCGA))->GenerateNullData();
                TCGA->~MainTable();
                break;
            case 3:
                TCGA = new MainTable();
                TCGA->read("nullTable.csv", true);
                TCGA->~MainTable();
                break;
            case 4:
                TCGA = new MainTable();
                TCGA->read("nullTable.csv", true);
                TCGA->ExtimateCorrelations("correlations_null.dat");
                TCGA->~MainTable();
            case 5:
                TCGA = new MainTable();
                TCGA->SaveMeansVariances("mainTable.csv");
                TCGA->~MainTable();
                break;
            case 6:
                TCGA = new MainTable();
                TCGA->read("mainTable.csv", false, true);
                TCGA->MakeCorpus();
                TCGA->~MainTable();
                break;
            default:
                std::cerr << "missing arguments" << std::endl;
                break;
        }

    }

    return 0;
}
