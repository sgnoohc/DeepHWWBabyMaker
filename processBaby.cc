#include "TChain.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"

#include "ScanChain.h"

#include "hadoopmap.h"

#include <iostream>

int main(int argc, char **argv)
{
    if (argc < 3)
    {
        std::cout << "USAGE: processBaby <tag> <filename> [<max_num_events>] [index] [type anything to make it verbose]" << std::endl;
        return 1;
    }
    TString infile(argv[1]);
    int max_events = -1;
    int index = 1;
    bool verbose = false;
    if (argc >= 3) { max_events = atoi(argv[2]); }
    if (argc >= 4) { index = atoi(argv[3]); }
    if (argc >= 5) { verbose = true; }
    std::cout << "set max number of events to: " << max_events << std::endl;
    std::cout << "running on files: " << infile.Data() << std::endl;

    TChain *chain = new TChain("Events");

    TString input_path = infile;
    TObjArray* files = input_path.Tokenize(",");
    for (int ifile = 0; ifile < files->GetEntries(); ++ifile)
    {
        TString filepath = ((TObjString*) files->At(ifile))->GetString();

        // Intervention needed for CMS3 data sitting on hadoop but mapped to nfs
        HadoopPathMap hmap;
        if (filepath.Contains("nfs"))
        {
            std::cout << "Adding to TChain: file = " << hmap.getHadoopPath(filepath) << std::endl;
            chain->Add(hmap.getHadoopPath(filepath));
        }
        else
        {
            std::cout << "Adding to TChain: file = " << filepath << std::endl;
            chain->Add(filepath.Data());
        }

    }

    if (chain->GetEntries() == 0)
    {
        std::cout << "ERROR: no entries in chain. filename was: " << infile << std::endl;
        return 2;
    }
    //--------------------------------
    // run
    //--------------------------------
    babyMaker *looper = new babyMaker();
    looper->ScanChain(chain, "dummy", max_events, index, verbose);
    return 0;
}
