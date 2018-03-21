#ifndef ScanChain_h
#define ScanChain_h

// C++
#include <iostream>
#include <vector>
#include <set>
#include <string>
#include <vector>
#include <algorithm>

// ROOT
#include "TChain.h"
#include "TDirectory.h"
#include "TTreeCache.h"
#include "Math/VectorUtil.h"
#include "TVector2.h"
#include "TBenchmark.h"
#include "TLorentzVector.h"
#include "TH2.h"
#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH2.h"
#include "TString.h"
#include "TMVA/Reader.h"
#include "Math/LorentzVector.h"
#include "Math/GenVector/LorentzVector.h"

// CORE
#include "CMS3.h"
#include "Base.h"
#include "OSSelections.h"
#include "SSSelections.h"
#include "VVVSelections.h"
#include "ElectronSelections.h"
#include "IsolationTools.h"
#include "JetSelections.h"
#include "MuonSelections.h"
#include "IsoTrackVeto.h"
#include "PhotonSelections.h"
#include "TriggerSelections.h"
#include "VertexSelections.h"
#include "MCSelections.h"
#include "MetSelections.h"
#include "SimPa.h"
#include "Tools/MT2/MT2.h"
#include "Tools/hemJet.h"
#include "Tools/utils.h"
#include "Tools/goodrun.h"
#include "Tools/btagsf/BTagCalibrationStandalone.h"
#include "Tools/btagsf/BTagCalibrationStandalone.h"
#include "Tools/datasetinfo/getDatasetInfo.h"

// RooUtil
#include "rooutil/looper.h"
#include "rooutil/ttreex.h"

// CoreUtil
#include "coreutil/jec.h"
#include "coreutil/btag.h"
#include "coreutil/puwgt.h"
#include "coreutil/genpart.h"
#include "coreutil/trigger.h"
#include "coreutil/electron.h"
#include "coreutil/muon.h"
#include "coreutil/grl.h"
#include "coreutil/datasetinfo.h"
#include "coreutil/jet.h"
#include "coreutil/met.h"
#include "coreutil/track.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

// `````````````````````````````````````````````````````````````````````````````````````````````````````````````
//
//
//
// Baby Maker Class
//
//
//
// ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

class babyMaker
{

private:

    SimPa simpa;
    CoreUtil::jec coreJec;
    CoreUtil::grl coreGRL;
    CoreUtil::btag coreBtagSF;
    CoreUtil::puwgt corePUWgt;
    CoreUtil::trigger coreTrigger;
    CoreUtil::genpart coreGenPart;
    CoreUtil::electron coreElectron;
    CoreUtil::muon coreMuon;
    CoreUtil::datasetinfo coreDatasetInfo;
    CoreUtil::jet coreJet;
    CoreUtil::met coreMET;
    CoreUtil::track coreTrack;

    TFile* ofile;
    TTree* t;
    RooUtil::TTreeX* tx;

    bool isData;

public:

    babyMaker() {}
    ~babyMaker() {}
    void ScanChain(TChain*, std::string = "testSample", int max_events = -1, int index = 1, bool verbose = false);

    void CreateOutput(int index=1);

    void ConfigureGoodRunsList();

    void ProcessTriggers();
    void ProcessGenParticles();
    void ProcessElectrons();
    void ProcessMuons();
    void ProcessJets();
    void ProcessMET();
    void ProcessTracks();

    bool PassPresel();

    void FillOutput();

    void FillEventInfo();
    void FillElectrons();
    void FillMuons();
    void FillJets();
    void FillMET();
    void FillTracks();
    void SortJetBranches();
    void FillVertexInfo();
    void FillMETFilter();
    void FillTTree();

    bool isLeptonOverlappingWithJet(int ijet);
    bool isLeptonOverlappingWithTrack(int ijet);
    static bool isLooseMuon(int);
    static bool isLooseElectron(int);
    static bool isPreselMuon(int);
    static bool isPreselElectron(int);
    static bool checkMuonTag(int, int);
    static bool checkElectronTag(int, int);

    // Calculator
    static int passCount(const vector<int>& vec);

};


#endif
