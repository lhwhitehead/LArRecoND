//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue May 16 12:52:46 2023 by ROOT version 6.22/08
// from TTree events/events
// found on file: ../2x2Prod/MiniRun3_1E19_RHC.flow_v3.00000.FLOWunMergedhits_withTruth.root
//////////////////////////////////////////////////////////

#ifndef LArSPMC_H
#define LArSPMC_H

#include "LArSP.h"

#include <TChain.h>
#include <TFile.h>
#include <TROOT.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

namespace lar_nd_reco
{

class LArSPMC : public LArSP
{
public:
    // Hit level truth information
    std::vector<std::vector<int>> *hit_segmentID;
    std::vector<std::vector<int>> *hit_segmentIndex;
    std::vector<std::vector<int>> *hit_particleID;
    std::vector<std::vector<int>> *hit_particleIndex;
    std::vector<std::vector<int>> *hit_interactionIndex;
    std::vector<std::vector<int>> *hit_pdg;
    std::vector<std::vector<float>> *hit_packetFrac;
    TBranch *b_hit_segmentID;
    TBranch *b_hit_segmentIndex;
    TBranch *b_hit_particleID;
    TBranch *b_hit_particleIndex;
    TBranch *b_hit_interactionIndex;
    TBranch *b_hit_pdg;
    TBranch *b_hit_packetFrac;

    // MC Particle information
    std::vector<float> *mcp_energy;
    std::vector<int> *mcp_pdg;
    std::vector<int> *mcp_nuid;
    std::vector<int> *mcp_id;
    std::vector<int> *mcp_mother;
    std::vector<float> *mcp_px;
    std::vector<float> *mcp_py;
    std::vector<float> *mcp_pz;
    std::vector<float> *mcp_startx;
    std::vector<float> *mcp_starty;
    std::vector<float> *mcp_startz;
    std::vector<float> *mcp_endx;
    std::vector<float> *mcp_endy;
    std::vector<float> *mcp_endz;
    TBranch *b_mcp_energy;
    TBranch *b_mcp_pdg;
    TBranch *b_mcp_nuid;
    TBranch *b_mcp_id;
    TBranch *b_mcp_mother;
    TBranch *b_mcp_px;
    TBranch *b_mcp_py;
    TBranch *b_mcp_pz;
    TBranch *b_mcp_startx;
    TBranch *b_mcp_starty;
    TBranch *b_mcp_startz;
    TBranch *b_mcp_endx;
    TBranch *b_mcp_endy;
    TBranch *b_mcp_endz;

    // Neutrino information
    std::vector<int> *nuID;
    std::vector<float> *nue;
    std::vector<int> *nuPDG;
    std::vector<float> *nupx;
    std::vector<float> *nupy;
    std::vector<float> *nupz;
    std::vector<float> *nuvtxx;
    std::vector<float> *nuvtxy;
    std::vector<float> *nuvtxz;
    TBranch *b_nuID;
    TBranch *b_nue;
    TBranch *b_nuPDG;
    TBranch *b_nupx;
    TBranch *b_nupy;
    TBranch *b_nupz;
    TBranch *b_nuvtxx;
    TBranch *b_nuvtxy;
    TBranch *b_nuvtxz;

    LArSPMC(TTree *tree = 0);
    virtual ~LArSPMC();
    virtual void InitMC(TTree *tree);
};

LArSPMC::LArSPMC(TTree *tree) : LArSP(tree)
{
    if (tree == 0)
    {
        std::cout << "Warning: null tree passed to LArSPMC" << std::endl;
    }
    InitMC(tree);
}

LArSPMC::~LArSPMC()
{
}

void LArSPMC::InitMC(TTree *tree)
{
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the branch addresses and branch
    // pointers of the tree will be set.
    // It is normally not necessary to make changes to the generated
    // code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running on PROOF
    // (once per file to be processed).

    // Set object pointer
    hit_segmentID = nullptr;
    hit_segmentIndex = nullptr;
    hit_particleID = nullptr;
    hit_particleIndex = nullptr;
    hit_interactionIndex = nullptr;
    hit_pdg = nullptr;
    hit_packetFrac = nullptr;

    mcp_energy = nullptr;
    mcp_pdg = nullptr;
    mcp_nuid = nullptr;
    mcp_id = nullptr;
    mcp_mother = nullptr;
    mcp_px = nullptr;
    mcp_py = nullptr;
    mcp_pz = nullptr;
    mcp_startx = nullptr;
    mcp_starty = nullptr;
    mcp_startz = nullptr;
    mcp_endx = nullptr;
    mcp_endy = nullptr;
    mcp_endz = nullptr;

    nuID = nullptr;
    nue = nullptr;
    nuPDG = nullptr;
    nupx = nullptr;
    nupy = nullptr;
    nupz = nullptr;
    nuvtxx = nullptr;
    nuvtxy = nullptr;
    nuvtxz = nullptr;

    // Set branch addresses and branch pointers
    if (!tree)
        return;
    fChain = tree;
    fCurrent = -1;
    fChain->SetMakeClass(1);

    fChain->SetBranchAddress("hit_segmentID", &hit_segmentID, &b_hit_segmentID);
    fChain->SetBranchAddress("hit_segmentIndex", &hit_segmentIndex, &b_hit_segmentIndex);
    fChain->SetBranchAddress("hit_particleID", &hit_particleID, &b_hit_particleID);
    fChain->SetBranchAddress("hit_particleIndex", &hit_particleIndex, &b_hit_particleIndex);
    fChain->SetBranchAddress("hit_interactionIndex", &hit_interactionIndex, &b_hit_interactionIndex);
    fChain->SetBranchAddress("hit_pdg", &hit_pdg, &b_hit_pdg);
    fChain->SetBranchAddress("hit_packetFrac", &hit_packetFrac, &b_hit_packetFrac);
    fChain->SetBranchAddress("mcp_energy", &mcp_energy, &b_mcp_energy);
    fChain->SetBranchAddress("mcp_pdg", &mcp_pdg, &b_mcp_pdg);
    fChain->SetBranchAddress("mcp_nuid", &mcp_nuid, &b_mcp_nuid);
    fChain->SetBranchAddress("mcp_id", &mcp_id, &b_mcp_id);
    fChain->SetBranchAddress("mcp_mother", &mcp_mother, &b_mcp_mother);
    fChain->SetBranchAddress("mcp_px", &mcp_px, &b_mcp_px);
    fChain->SetBranchAddress("mcp_py", &mcp_py, &b_mcp_py);
    fChain->SetBranchAddress("mcp_pz", &mcp_pz, &b_mcp_pz);
    fChain->SetBranchAddress("mcp_startx", &mcp_startx, &b_mcp_startx);
    fChain->SetBranchAddress("mcp_starty", &mcp_starty, &b_mcp_starty);
    fChain->SetBranchAddress("mcp_startz", &mcp_startz, &b_mcp_startz);
    fChain->SetBranchAddress("mcp_endx", &mcp_endx, &b_mcp_endx);
    fChain->SetBranchAddress("mcp_endy", &mcp_endy, &b_mcp_endy);
    fChain->SetBranchAddress("mcp_endz", &mcp_endz, &b_mcp_endz);
    fChain->SetBranchAddress("nuID", &nuID, &b_nuID);
    fChain->SetBranchAddress("nue", &nue, &b_nue);
    fChain->SetBranchAddress("nuPDG", &nuPDG, &b_nuPDG);
    fChain->SetBranchAddress("nupx", &nupx, &b_nupx);
    fChain->SetBranchAddress("nupy", &nupy, &b_nupy);
    fChain->SetBranchAddress("nupz", &nupz, &b_nupz);
    fChain->SetBranchAddress("nuvtxx", &nuvtxx, &b_nuvtxx);
    fChain->SetBranchAddress("nuvtxy", &nuvtxy, &b_nuvtxy);
    fChain->SetBranchAddress("nuvtxz", &nuvtxz, &b_nuvtxz);
}

} // namespace lar_nd_reco

#endif // #ifdef LArSPMC_H
