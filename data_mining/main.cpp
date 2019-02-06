


#include "MainTable.h"
#include "TelegramWatch.h"
#include "NullModel.h"
#include <omp.h>
#include <iostream>


int main(int argc, const char** argv) // minimal test driver
{
    omp_set_nested(1);

    TelegramWatch watch("thesis");

    MainTable* TCGA;

    switch(std::atoi(argv[1])){
        case 0:
            TCGA = new MainTable();
            TCGA->read("mainTable.csv", true);
//    TCGA->SaveBinary("binaryTable.csv");
//    TCGA->readBinary();
            TCGA->ExtimateCorrelations();
            TCGA->~MainTable();
            break;
        case 1:
            TCGA = new NullModel();
//    ((NullModel*)(TCGA))->GenerateNullBinary();
            ((NullModel*)(TCGA))->GenerateNullData();
            //     ((NullModel*)(TCGA))->SaveBinary("binaryNull.csv");
//     ((NullModel*)(TCGA))->read("binaryNull.csv", true);
//    ((NullModel*)(TCGA))->ExtimateCorrelations("correlations_null.dat");
            TCGA->~MainTable();
            break;

        case 2:
            TCGA = new MainTable();
            TCGA->read("nullTable.csv", true);
            TCGA->ExtimateCorrelations("correlations_null.dat");
            TCGA->~MainTable();
            break;
        default:
            std::cerr<<"missing arguments"<<std::endl;
            break;
    }



    return 0;
}
