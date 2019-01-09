//
// Created by Filippo Valle on 04/11/2018.
//

#ifndef NEURALANALYZER_TELEGRAMWATCH_H
#define NEURALANALYZER_TELEGRAMWATCH_H

#include <string>
#include <algorithm>
#include <iostream>
#include <boost/timer/timer.hpp>
#include <boost/chrono/include.hpp>
#include <curl/curl.h>                        // strumento per interagire con il web

using std::string;

class TelegramWatch : public boost::timer::auto_cpu_timer {      // classe che eredita da StopWatch
public:
    TelegramWatch();
    TelegramWatch(const char* name);
    ~TelegramWatch();

private:
    void callIFTTT() const;                   // IFTTT = If This Then That. Servizio per mettere in comunicazione + dispositivi. Es imposto di chiamare il mio cellulare.
    static const char* fURL;               
    const char* fName;
};


#endif //NEURALANALYZER_TELEGRAMWATCH_H
