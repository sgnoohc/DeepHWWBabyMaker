#ifndef hadoopmap_h
#define hadoopmap_h

#include <map>
#include <iostream>
#include <fstream>
#include <string>
#include "TString.h"

using namespace std;

class HadoopPathMap
{
public:
    std::map<TString, TString> _map;
    HadoopPathMap();
    ~HadoopPathMap();
    TString getHadoopPath( TString );
};

#endif
