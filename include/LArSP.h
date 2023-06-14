//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue May 16 12:52:46 2023 by ROOT version 6.22/08
// from TTree events/events
// found on file: ../2x2Prod/MiniRun3_1E19_RHC.flow_v3.00000.FLOWunMergedhits_withTruth.root
//////////////////////////////////////////////////////////

#ifndef LArSP_H
#define LArSP_H

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

namespace lar_nd_reco
{

class LArSP {
public :
    TTree          *fChain;   //!pointer to the analyzed TTree or TChain
    Int_t           fCurrent; //!current Tree number in a TChain

    // Declaration of leaf types
    Int_t           event;
    Int_t           subrun;
    Int_t           run;
    Int_t           event_start_t;
    Int_t           event_end_t;
    std::vector<float>   *x;
    std::vector<float>   *y;
    std::vector<float>   *z;
    std::vector<float>   *ts;
    std::vector<float>   *charge;
    std::vector<float>   *E;

    // List of branches
    TBranch        *b_eventID;   //!
    TBranch        *b_subrun;   //!
    TBranch        *b_run;   //!
    TBranch        *b_event_start_t;   //!
    TBranch        *b_event_end_t;   //!
    TBranch        *b_x;   //!
    TBranch        *b_y;   //!
    TBranch        *b_z;   //!
    TBranch        *b_ts;   //!
    TBranch        *b_charge;   //!
    TBranch        *b_E;   //!

    LArSP(TTree *tree);
    virtual ~LArSP();
    virtual Int_t    GetEntry(Long64_t entry);
    virtual Long64_t LoadTree(Long64_t entry);
    virtual void     Init(TTree *tree);
    virtual void     InitMC(TTree *tree);
    virtual Bool_t   Notify();
    virtual void     Show(Long64_t entry = -1);
};

LArSP::LArSP(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
    if (tree == 0) {
        TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../2x2Prod/MiniRun3_1E19_RHC.flow_v3.00000.FLOWunMergedhits_withTruth.root");
        if (!f || !f->IsOpen()) {
            f = new TFile("../2x2Prod/MiniRun3_1E19_RHC.flow_v3.00000.FLOWunMergedhits_withTruth.root");
        }
        f->GetObject("events",tree);
    }
   Init(tree);
}

LArSP::~LArSP()
{
    if (!fChain)
        return;
    delete fChain->GetCurrentFile();
}

Int_t LArSP::GetEntry(Long64_t entry)
{
    // Read contents of entry.
    if (!fChain)
        return 0;
    return fChain->GetEntry(entry);
}

Long64_t LArSP::LoadTree(Long64_t entry)
{
    // Set the environment to read one entry
    if (!fChain)
        return -5;
    Long64_t centry = fChain->LoadTree(entry);
    if (centry < 0)
        return centry;
    if (fChain->GetTreeNumber() != fCurrent)
    {
        fCurrent = fChain->GetTreeNumber();
        Notify();
    }
    return centry;
}

void LArSP::Init(TTree *tree)
{
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the branch addresses and branch
    // pointers of the tree will be set.
    // It is normally not necessary to make changes to the generated
    // code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running on PROOF
    // (once per file to be processed).
    
    // Set object pointer
    x = nullptr;
    y = nullptr;
    z = nullptr;
    ts = nullptr;
    charge = nullptr;
    E = nullptr;
    // Set branch addresses and branch pointers
    if (!tree)
        return;
    fChain = tree;
    fCurrent = -1;
    fChain->SetMakeClass(1);
    
    fChain->SetBranchAddress("event", &event, &b_eventID);
    fChain->SetBranchAddress("subrun", &subrun, &b_subrun);
    fChain->SetBranchAddress("run", &run, &b_run);
    fChain->SetBranchAddress("event_start_t", &event_start_t, &b_event_start_t);
    fChain->SetBranchAddress("event_end_t", &event_end_t, &b_event_end_t);
    fChain->SetBranchAddress("x", &x, &b_x);
    fChain->SetBranchAddress("y", &y, &b_y);
    fChain->SetBranchAddress("z", &z, &b_z);
    fChain->SetBranchAddress("ts", &ts, &b_ts);
    fChain->SetBranchAddress("charge", &charge, &b_charge);
    fChain->SetBranchAddress("E", &E, &b_E);
    Notify();
}

void LArSP::InitMC(TTree*)
{
}

Bool_t LArSP::Notify()
{
    // The Notify() function is called when a new file is opened. This
    // can be either for a new TTree in a TChain or when when a new TTree
    // is started when using PROOF. It is normally not necessary to make changes
    // to the generated code, but the routine can be extended by the
    // user if needed. The return value is currently not used.
    
    return kTRUE;
}

void LArSP::Show(Long64_t entry)
{
    // Print contents of entry.
    // If entry is not specified, print current entry
    if (!fChain) return;
    fChain->Show(entry);
}

} // end namespace lar_nd_reco

#endif // #ifdef LArSP_H
