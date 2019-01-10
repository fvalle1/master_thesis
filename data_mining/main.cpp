


#include "MainTable.h"
#include "TelegramWatch.h"


int main() // minimal test driver
{
    TelegramWatch watch("thesis");

    auto TCGAData = MainTable();
//    TCGAData.read();
//    TCGAData.SaveBinary();

    TCGAData.readBinary();
    TCGAData.ExtimateCorrelations();

    return 0;
}