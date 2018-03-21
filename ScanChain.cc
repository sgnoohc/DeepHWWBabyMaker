#include "ScanChain.h"

using namespace std;

//##############################################################################################################
void babyMaker::ScanChain(TChain* chain, std::string baby_name, int max_events, int index, bool verbose)
{

    // Looper
    RooUtil::Looper<CMS3> looper(chain, &cms3, max_events);

    // Output root file
    CreateOutput(index);

    int nskipped_batch = 0;
    int nskipped = 0;
    while (looper.nextEvent())
    {
        try
        {
            if (verbose)
                cout << "[verbose] Processed " << looper.getNEventsProcessed() << " out of " << chain->GetEntries() << endl;
            
            coreJec.setJECFor(looper.getCurrentFileName());

            // Loop over gen level particles
            ProcessGenParticles();

            // Fill baby ntuple
            FillOutput();
        }
        catch (const std::ios_base::failure& e)
        {
            cout << endl;
            cout << "[CheckCorrupt] Caught an I/O failure in the ROOT file." << endl;
            cout << "[CheckCorrupt] Possibly corrupted hadoop file." << endl;
            cout << "[CheckCorrupt] event index = " << looper.getCurrentEventIndex() << " out of " << chain->GetEntries() << endl;
            cout << e.what() << endl;
            cout << endl;
            tx->clear(); // clear the TTree of any residual stuff

            nskipped_batch++;

            // If the nskipped is quite large than skip the entire file
            if (nskipped_batch > 500)
            {
                nskipped += nskipped_batch;
                nskipped_batch = 0;
                for (int i = 0; i < 10000; ++i)
                {
                    if (!looper.nextEvent())
                        break;
                    nskipped++;
                }
            }
        }
    }

    nskipped += nskipped_batch;

    looper.getTree()->PrintCacheStats();

    if (nskipped)
    {
        cout << "[CheckCorrupt] Skipped " << nskipped << " events out of " << chain->GetEntries() << " [" << float(nskipped) / float(chain->GetEntries()) << "%% loss]" << endl;
    }

    ofile->cd();
    t->Write();
}

//##############################################################################################################
void babyMaker::CreateOutput(int index)
{
    ofile = new TFile(Form("output_%d.root", index), "recreate");
    t = new TTree("t", "t");
    tx = new RooUtil::TTreeX(t);

    tx->createBranch<vector<LV>>("genPart_p4");
    tx->createBranch<vector<int>>("genPart_pdgId");
    tx->createBranch<vector<int>>("genPart_status");
    tx->createBranch<vector<int>>("genPart_motherId");
    tx->createBranch<vector<int>>("genPart_grandmaId");

    tx->clear();
}

//##############################################################################################################
void babyMaker::ProcessTriggers() { coreTrigger.process(); }

//##############################################################################################################
void babyMaker::ProcessGenParticles() { coreGenPart.process(isPreselGenPart); }

//##############################################################################################################
void babyMaker::ProcessElectrons() { coreElectron.process(isPreselElectron); }

//##############################################################################################################
void babyMaker::ProcessMuons() { coreMuon.process(isPreselMuon); }

//##############################################################################################################
void babyMaker::ProcessJets() { coreJet.process(coreJec); }

//##############################################################################################################
void babyMaker::ProcessMET() { coreMET.process(coreJec); }

//##############################################################################################################
void babyMaker::ProcessTracks() { coreTrack.process(); }

//##############################################################################################################
bool babyMaker::PassPresel()
{
    // Why? I don't know
    if (cms3.evt_isRealData() && cms3.evt_run() <= 251562)
        return false;

    if (cms3.evt_isRealData() && !goodrun(cms3.evt_run(), cms3.evt_lumiBlock()))
        return false;

    // Place your preselection
    return true;
}

//##############################################################################################################
void babyMaker::FillOutput()
{
    // Fill the generator information
    FillGenInfo();

    // Actually fill the TTree
    FillTTree();
}

//##############################################################################################################
void babyMaker::FillEventInfo()
{
}

//##############################################################################################################
void babyMaker::FillElectrons()
{
}

//##############################################################################################################
void babyMaker::FillMuons()
{
}

//##############################################################################################################
void babyMaker::FillJets()
{
}

//##############################################################################################################
void babyMaker::FillMET()
{
}

//##############################################################################################################
void babyMaker::SortJetBranches()
{
//    tx->sortVecBranchesByPt("jets_p4", {"jets_csv"}, {}, {});
}

//##############################################################################################################
void babyMaker::FillVertexInfo()
{
}

//##############################################################################################################
void babyMaker::FillGenInfo()
{
    for (unsigned int ip = 0; ip < coreGenPart.genPart_p4.size(); ++ip)
    {
        tx->pushbackToBranch<LV>("genPart_p4", coreGenPart.genPart_p4[ip]);
        tx->pushbackToBranch<int>("genPart_pdgId", coreGenPart.genPart_pdgId[ip]);
        tx->pushbackToBranch<int>("genPart_status", coreGenPart.genPart_status[ip]);
        tx->pushbackToBranch<int>("genPart_motherId", coreGenPart.genPart_motherId[ip]);
        tx->pushbackToBranch<int>("genPart_grandmaId", coreGenPart.genPart_grandmaId[ip]);
    }
}

//##############################################################################################################
void babyMaker::FillMETFilter()
{
}

//##############################################################################################################
void babyMaker::FillTTree()
{
    tx->fill();
    tx->clear();
}

//##############################################################################################################
bool babyMaker::isLeptonOverlappingWithJet(int ijet)
{
    bool is_overlapping = false;

    int idx = coreJet.index[ijet];

    for (auto& imu : coreMuon.index)
    {
        if (!(isLooseMuon(imu)))
            continue;

        if (ROOT::Math::VectorUtil::DeltaR(cms3.pfjets_p4()[idx], cms3.mus_p4()[imu]) < 0.4)
        {
            is_overlapping = true;
            break;
        }
    }

    if (is_overlapping)
        return true;

    for (auto& iel : coreElectron.index)
    {
        if (!(isLooseElectron(iel)))
            continue;

        if (ROOT::Math::VectorUtil::DeltaR(cms3.pfjets_p4()[idx], cms3.els_p4()[iel]) < 0.4)
        {
            is_overlapping = true;
            break;
        }
    }

    if (is_overlapping)
        return true;

    return false;
}

//##############################################################################################################
// Used to overlap remova against jets
bool babyMaker::isLooseMuon(int idx)
{
    if (!( cms3.mus_p4()[idx].pt() > 20.                     )) return false;
    if (!( passMuonSelection_VVV(idx, VVV_cutbased_3l_fo_v2) )) return false;
    return true;
}

//##############################################################################################################
// Used to overlap remova against jets
bool babyMaker::isLooseElectron(int idx)
{
    if (!( cms3.els_p4()[idx].pt() > 20.                         )) return false;
    if (!( passElectronSelection_VVV(idx, VVV_cutbased_3l_fo_v2) )) return false;
    return true;
}

//##############################################################################################################
// Used to overlap remova against tracks
bool babyMaker::isPreselMuon(int idx)
{
    if (!( cms3.mus_p4()[idx].pt() >= 10. )) return false;
    return true;
}

//##############################################################################################################
// Used to overlap remova against tracks
bool babyMaker::isPreselElectron(int idx)
{
    if (!( cms3.els_p4()[idx].pt() >= 10.        )) return false;
    return true;
}

//##############################################################################################################
// Used to select OK gen particles
bool babyMaker::isPreselGenPart(int idx)
{
    // Good for powheg pythia8 samples
    int status = cms3.genps_status()[idx];
    int id = cms3.genps_id()[idx];
    int mother = cms3.genps_id_mother()[idx];
    if (!(
        /* is H or W         */     (status == 22 && (id == 25 || abs(id) == 24))
        /* is W decay quarks */  || (status == 23 && abs(mother) == 24 && (abs(id) >= 1 && abs(id) <= 6)) // allow bottom and top just for completion.
        /* is W decay leptons*/  || (status ==  1 && abs(mother) == 24 && (abs(id) >= 11 && abs(id) <= 16))
        /* is W decay tau    */  || (status ==  2 && abs(mother) == 24 && (abs(id) == 15))
        )) return false;
    return true;
}

//eof

