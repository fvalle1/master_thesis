//
// Created by Filippo Valle on 04/11/2018.
//

#include "TelegramWatch.h"

const char* TelegramWatch::fURL = "https://maker.ifttt.com/trigger/cpp/with/key/95ZnJdPvDZmBx8yfEy6WL";   //puntatore ad una stringa di caratteri costante

TelegramWatch::TelegramWatch(): boost::timer::auto_cpu_timer(),
                                fName("")
{
  
}

TelegramWatch::TelegramWatch(const char *name) : boost::timer::auto_cpu_timer(),
                                                 fName(name) {
  
						 }

TelegramWatch::~TelegramWatch(){             // quando distruggo TelegramWatch faccio la chiamata
  callIFTTT();
  std::cout<<std::endl;
}

void TelegramWatch::callIFTTT() const {      
  CURL *curl;
  CURLcode res;
  
  curl_global_init(CURL_GLOBAL_DEFAULT);   // devo inizializzare la libreria curl prima di poter usare qualsiasi sua funzione
  
  curl = curl_easy_init();                 // inizio una sessione di libcurl
  if(curl) {                               // se curl (che è un ptr) punta a qualcosa , allora:
    char encWall[35];                      // encoded Wall: conterrà il valore tempo REALE, però scritto in una stringa
    char encSystem[35];                       // encoded CPU: tempo della CPU scritto in caratteri
    char url[225];
    
    char* encName = curl_easy_escape(curl, fName, 7);            // This function converts the given input string to a URL encoded string and returns that as a new allocated string.
    sprintf(encWall, "%6.2fsecs[Wall_time]", boost::chrono::duration_cast<boost::chrono::milliseconds>(boost::chrono::nanoseconds(elapsed().wall)).count()/1000.);     // Scrivo in caratteri il tempo reale (troncato alla seconda cifra decimale)
    sprintf(encSystem, "%6.2fsecs[Sys_time]", boost::chrono::duration_cast<boost::chrono::milliseconds>(boost::chrono::nanoseconds(elapsed().user)).count()/1000.);
    sprintf(url, "%s?value1=%s&value2=%s&value3=%s", fURL, encName, curl_easy_escape(curl, encWall, 19), curl_easy_escape(curl, encSystem, 18));
                                                                  // scrivo: URL + nome + tempo REALE (in formato url) + tempo CPU (in formato URL) 
    // printf("%s\n", url);
    curl_easy_setopt(curl, CURLOPT_URL, url);                     // Pass in a pointer to the URL to work with. The parameter should be a char * to a zero terminated string which must be URL-encoded
    
#ifdef SKIP_PEER_VERIFICATION
    /*
     * If you want to connect to a site who isn't using a certificate that is
     * signed by one of the certs in the CA bundle you have, you can skip the
     * verification of the server's certificate. This makes the connection
     * A LOT LESS SECURE.
     *
     * If you have a CA cert for the server stored someplace else than in the
     * default bundle, then the CURLOPT_CAPATH option might come handy for
     * you.
     */
    curl_easy_setopt(curl, CURLOPT_SSL_VERIFYPEER, 0L);
#endif
    
#ifdef SKIP_HOSTNAME_VERIFICATION
    /*
     * If the site you're connecting to uses a different host name that what
     * they have mentioned in their server certificate's commonName (or
     * subjectAltName) fields, libcurl will refuse to connect. You can skip
     * this check, but this will make the connection less secure.
     */
    curl_easy_setopt(curl, CURLOPT_SSL_VERIFYHOST, 0L);            // Invoke this function after curl_easy_init and all the curl_easy_setopt calls are made, and will perform the transfer as described in                                                                   // the options
#endif
    
    /* Perform the request, res will get the return code */
    res = curl_easy_perform(curl);
    /* Check for errors */
    if(res != CURLE_OK)
      printf("curl_easy_perform() failed: %s\n", curl_easy_strerror(res));
    
    /* always cleanup */
    curl_easy_cleanup(curl);
  }
  
  curl_global_cleanup();
}

