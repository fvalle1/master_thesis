


#include "MainTable.h"
#include "TelegramWatch.h"


int main() // minimal test driver
{
    TelegramWatch watch("thesis");

    auto TCGAData = MainTable();
    TCGAData.read();
    TCGAData.SaveBinary();

    return 0;
}