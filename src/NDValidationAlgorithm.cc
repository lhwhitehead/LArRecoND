/**
 *  @file   src/NDValidationAlgorithm.cc
 *
 *  @brief  Implementation of the particle visualisation algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "NDValidationAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArInteractionTypeHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

using namespace pandora;

namespace lar_content
{

NDValidationAlgorithm::MatchInfo::MatchInfo(const int pdg, const int nHits, const int nSharedHits, const float completeness, const float purity):
    m_pdg{pdg},
    m_nHits{nHits},
    m_nSharedHits{nSharedHits},
    m_completeness{completeness},
    m_purity{purity}
{}

//------------------------------------------------------------------------------------------------------------------------------------------

void NDValidationAlgorithm::MatchInfo::Print() const
{
    std::cout << "    - Pfo: pdg = " << m_pdg << ", hits = " << m_nHits << ", shared hits = " << m_nSharedHits << ", completeness = " << m_completeness << ", purity = " << m_purity << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

NDValidationAlgorithm::NDValidationAlgorithm() :
    m_event{-1},
    m_writeTree{false},
    m_foldToPrimaries{false},
    m_foldDynamic{false},
    m_foldToLeadingShowers{false},
    m_validateEvent{false},
    m_validateMC{false},
    m_printToScreen{false}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

NDValidationAlgorithm::~NDValidationAlgorithm()
{
    if (m_writeTree)
    {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treename.c_str(), m_filename.c_str(), "UPDATE"));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode NDValidationAlgorithm::Run()
{
    ++m_event;
    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));
    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));
    const PfoList *pPfoList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_pfoListName, pPfoList));

    m_matchMap.clear();
    m_mcHitsMap.clear();
    std::map<const MCParticle *, MCParticleList> trueNeutrinoMap;

    for (const MCParticle *pMCParticle : *pMCParticleList)
    {
        if (LArMCParticleHelper::IsNeutrino(pMCParticle) && !LArMCParticleHelper::GetHierarchyTier(pMCParticle))
        {
            MCParticleList descendentsList;
            LArMCParticleHelper::GetAllDescendentMCParticles(pMCParticle, descendentsList);
            trueNeutrinoMap[pMCParticle] = descendentsList;

            if(m_printToScreen)
                std::cout << "Found true neutrino " << pMCParticle << " with " << trueNeutrinoMap[pMCParticle].size() << " primaries" << std::endl;
        }
    }

    std::map<const ParticleFlowObject *, PfoList> recoNeutrinoMap;
    for (const ParticleFlowObject *pPfo : *pPfoList)
    {
        if (LArPfoHelper::IsNeutrino(pPfo))
        {
            PfoList descendentsList;
            LArPfoHelper::GetAllDownstreamPfos(pPfo, descendentsList);
            recoNeutrinoMap[pPfo] = descendentsList;
        }
    }

    LArHierarchyHelper::FoldingParameters foldParameters;
    if (m_foldToPrimaries)
        foldParameters.m_foldToTier = true;
    else if (m_foldDynamic)
        foldParameters.m_foldDynamic = true;
    else if (m_foldToLeadingShowers)
        foldParameters.m_foldToLeadingShowers = true;

    for (auto const &rPair : recoNeutrinoMap)
    {
        LArHierarchyHelper::RecoHierarchy recoHierarchy;
        LArHierarchyHelper::FillRecoHierarchy(rPair.second, foldParameters, recoHierarchy);
        for (auto const &pair : trueNeutrinoMap)
        {
            LArHierarchyHelper::MCHierarchy mcHierarchy;
            LArHierarchyHelper::FillMCHierarchy(pair.second, *pCaloHitList, foldParameters, mcHierarchy);
            LArHierarchyHelper::MatchInfo matchInfo;
            LArHierarchyHelper::MatchHierarchies(mcHierarchy, recoHierarchy, matchInfo);
            //matchInfo.Print(mcHierarchy);

            if (m_validateEvent)
                this->EventValidation(matchInfo);
            else if (m_validateMC)
                this->MCValidation(matchInfo);
        }
    }


    if (m_validateMC && m_printToScreen)
    {
        for (auto const &mcNeutrinoList : trueNeutrinoMap)
        {
            std::cout << std::endl;
            std::cout << "===== Matching Information for Neutrino " << static_cast<int>(reinterpret_cast<std::intptr_t>(mcNeutrinoList.first->GetUid())) << " =====" << std::endl;
            const MCParticleList &nuDescendents = mcNeutrinoList.second;
            unsigned int nReconstructableChildren{0};
            for (auto const &mcMatchPair : m_matchMap)
            { 
                const MCParticle *pMC{mcMatchPair.first};
                if (std::find(nuDescendents.begin(),nuDescendents.end(),pMC) == nuDescendents.end())
                    continue;

                ++nReconstructableChildren;
                const int mcId{static_cast<int>(reinterpret_cast<std::intptr_t>(pMC->GetUid()))};
                const int pdg{pMC->GetParticleId()};
                const int mcHits{m_mcHitsMap.at(pMC)};
                std::cout << "  - MC Particle: id = " << mcId << ", pdg = " << pdg << ", nhits = " << mcHits << ", matches " << mcMatchPair.second.size() << ":" << std::endl;
                for (const MatchInfo &match : mcMatchPair.second)
                    match.Print();
            }
            if (!nReconstructableChildren)
                std::cout << "- No true particles to reconstruct in the detector volume" << std::endl;
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NDValidationAlgorithm::EventValidation(const LArHierarchyHelper::MatchInfo &matchInfo) const
{
    if (m_writeTree)
    {
        const LArHierarchyHelper::MCMatchesVector &matches{matchInfo.GetMatches()};
        MCParticleSet primaryMCSet;
        std::set<const LArHierarchyHelper::MCHierarchy::Node *> trackNodeSet, showerNodeSet;
        int nGoodTrackMatches{0}, nGoodShowerMatches{0};
        int nGoodMatches{0}, nPoorMatches{0}, nUnmatched{0};
        int nGoodTier1Matches{0}, nTier1Nodes{0};
        int nGoodTier1TrackMatches{0}, nTier1TrackNodes{0};
        int nGoodTier1ShowerMatches{0}, nTier1ShowerNodes{0};
        int hasLeadingMuon{0}, hasLeadingElectron{0}, isLeadingLeptonCorrect{0};
        // ATTN: Probably want quality cuts here for "good match" definition
        for (const LArHierarchyHelper::MCMatches &mcMatch : matches)
        {
            const LArHierarchyHelper::MCHierarchy::Node *pNode{mcMatch.GetMC()};
            const MCParticle *const pMC{pNode->GetLeadingMCParticle()};
            primaryMCSet.insert(LArMCParticleHelper::GetPrimaryMCParticle(pMC));
            const int nReco{static_cast<int>(mcMatch.GetRecoMatches().size())};
            const bool isQuality{mcMatch.IsQuality(matchInfo.GetQualityCuts())};
            if (nReco == 1 && isQuality)
                ++nGoodMatches;
            else if (nReco == 1)
                ++nPoorMatches;
            else if (nReco > 1)
                ++nPoorMatches;
            else
                ++nUnmatched;
            if (pNode->GetHierarchyTier() == 1)
            {
                ++nTier1Nodes;
                if (nReco == 1 && isQuality)
                    ++nGoodTier1Matches;
            }

            const int pdg{std::abs(pNode->GetParticleId())};
            if (pNode->IsLeadingLepton())
            {
                if (pdg == MU_MINUS)
                    hasLeadingMuon = 1;
                else if (pdg == E_MINUS)
                    hasLeadingElectron = 1;
                isLeadingLeptonCorrect = nReco == 1 ? 1 : 0;
            }

            if (pdg == PHOTON || pdg == E_MINUS)
            {
                showerNodeSet.insert(pNode);
                if (nReco == 1 && isQuality)
                {
                    ++nGoodShowerMatches;
                    if (pNode->GetHierarchyTier() == 1)
                        ++nGoodTier1ShowerMatches;
                }
                if (pNode->GetHierarchyTier() == 1)
                    ++nTier1ShowerNodes;
            }
            else
            {
                trackNodeSet.insert(pNode);
                if (nReco == 1 && isQuality)
                {
                    ++nGoodTrackMatches;
                    if (pNode->GetHierarchyTier() == 1)
                        ++nGoodTier1TrackMatches;
                }
                if (pNode->GetHierarchyTier() == 1)
                    ++nTier1TrackNodes;
            }
        }

        MCParticleList primaryMCList;
        for (const MCParticle *const pMC : primaryMCSet)
            primaryMCList.emplace_back(pMC);
        const int interactionType{static_cast<int>(LArInteractionTypeHelper::GetInteractionType(primaryMCList))};
        const int nNodes{static_cast<int>(matchInfo.GetNMCNodes())};

        const int nTrackNodes{static_cast<int>(trackNodeSet.size())}, nShowerNodes{static_cast<int>(showerNodeSet.size())};
        const CartesianVector &trueVertex{matchInfo.GetMCNeutrino()->GetVertex()};
        const CartesianVector &recoVertex{LArPfoHelper::GetVertex(matchInfo.GetRecoNeutrino())->GetPosition()};
        const float vtxDx{recoVertex.GetX() - trueVertex.GetX()};
        const float vtxDy{recoVertex.GetY() - trueVertex.GetY()};
        const float vtxDz{recoVertex.GetZ() - trueVertex.GetZ()};
        const float vtxDr{std::sqrt(vtxDx * vtxDx + vtxDy * vtxDy + vtxDz * vtxDz)};

        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "event", m_event));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "interactionType", interactionType));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nGoodMatches", nGoodMatches));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nPoorMatches", nPoorMatches));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nUnmatched", nUnmatched));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nNodes", nNodes));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nGoodTier1Matches", nGoodTier1Matches));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nTier1Nodes", nTier1Nodes));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nGoodTrackMatches", nGoodTrackMatches));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nGoodShowerMatches", nGoodShowerMatches));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nTrackNodes", nTrackNodes));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nShowerNodes", nShowerNodes));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nGoodTier1TrackMatches", nGoodTier1TrackMatches));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nTier1TrackNodes", nTier1TrackNodes));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nGoodTier1ShowerMatches", nGoodTier1ShowerMatches));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nTier1ShowerNodes", nTier1ShowerNodes));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "hasLeadingMuon", hasLeadingMuon));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "hasLeadingElectron", hasLeadingElectron));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isLeadingLeptonCorrect", isLeadingLeptonCorrect));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "vtxDx", vtxDx));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "vtxDy", vtxDy));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "vtxDz", vtxDz));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "vtxDr", vtxDr));
        PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treename.c_str()));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NDValidationAlgorithm::MCValidation(const LArHierarchyHelper::MatchInfo &matchInfo)
{
    if (m_writeTree)
    {
        for (const LArHierarchyHelper::MCMatches &match : matchInfo.GetMatches())
            this->Fill(match, matchInfo);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NDValidationAlgorithm::Fill(const LArHierarchyHelper::MCMatches &matches, const LArHierarchyHelper::MatchInfo &matchInfo)
{
    const LArHierarchyHelper::MCHierarchy::Node *pMCNode{matches.GetMC()};
    const MCParticle *pMCParticle{pMCNode->GetLeadingMCParticle()};
    const int isTestBeam{pMCNode->IsTestBeamParticle() ? 1 : 0};
    const int isCosmicRay{!isTestBeam && pMCNode->IsCosmicRay() ? 1 : 0};
    const int isNeutrinoInt{!(isTestBeam || isCosmicRay) ? 1 : 0};
    const int mcId{static_cast<int>(reinterpret_cast<std::intptr_t>(pMCParticle->GetUid()))};
    const int pdg{pMCNode->GetParticleId()};
    const int tier{pMCNode->GetHierarchyTier()};
    const int mcHits{static_cast<int>(pMCNode->GetCaloHits().size())};
    const int isLeadingLepton{pMCNode->IsLeadingLepton() ? 1 : 0};

    const MCParticleList &parentList{pMCParticle->GetParentList()};
    const int isElectron{std::abs(pMCParticle->GetParticleId()) == E_MINUS ? 1 : 0};
    const int hasMuonParent{parentList.size() == 1 && std::abs(parentList.front()->GetParticleId()) == MU_MINUS ? 1 : 0};
    const int isMichel{isElectron && hasMuonParent && LArMCParticleHelper::IsDecay(pMCParticle) ? 1 : 0};

    const LArHierarchyHelper::RecoHierarchy::NodeVector &nodeVector{matches.GetRecoMatches()};
    const int nMatches{static_cast<int>(nodeVector.size())};
    IntVector recoIdVector, nRecoHitsVector, nSharedHitsVector;
    FloatVector purityVector, completenessVector;
    FloatVector purityAdcVector, completenessAdcVector;
    FloatVector purityVectorU, purityVectorV, purityVectorW, completenessVectorU, completenessVectorV, completenessVectorW;
    FloatVector purityAdcVectorU, purityAdcVectorV, purityAdcVectorW, completenessAdcVectorU, completenessAdcVectorV, completenessAdcVectorW;
    const CartesianVector &trueVertex{pMCNode->GetLeadingMCParticle()->GetVertex()};
    float vtxDx{0.f}, vtxDy{0.f}, vtxDz{0.f}, vtxDr{0.f};

    m_mcHitsMap[pMCParticle] = mcHits;

    for (const LArHierarchyHelper::RecoHierarchy::Node *pRecoNode : nodeVector)
    {
        const int recoId{pRecoNode->GetParticleId()};
        const int nRecoHits{static_cast<int>(pRecoNode->GetCaloHits().size())};
        const int nSharedHits{static_cast<int>(matches.GetSharedHits(pRecoNode))};
        const float completeness{matches.GetCompleteness(pRecoNode)};
        const float purity{matches.GetPurity(pRecoNode)};
        m_matchMap[pMCParticle].emplace_back(MatchInfo(recoId, nRecoHits, nSharedHits, completeness, purity));

        recoIdVector.emplace_back(recoId);
        nRecoHitsVector.emplace_back(nRecoHits);
        nSharedHitsVector.emplace_back(nSharedHits);
        purityVector.emplace_back(completeness);
        completenessVector.emplace_back(purity);
        purityAdcVector.emplace_back(matches.GetPurity(pRecoNode, true));
        completenessAdcVector.emplace_back(matches.GetCompleteness(pRecoNode, true));
        purityVectorU.emplace_back(matches.GetPurity(pRecoNode, TPC_VIEW_U));
        purityVectorV.emplace_back(matches.GetPurity(pRecoNode, TPC_VIEW_V));
        purityVectorW.emplace_back(matches.GetPurity(pRecoNode, TPC_VIEW_W));
        completenessVectorU.emplace_back(matches.GetCompleteness(pRecoNode, TPC_VIEW_U));
        completenessVectorV.emplace_back(matches.GetCompleteness(pRecoNode, TPC_VIEW_V));
        completenessVectorW.emplace_back(matches.GetCompleteness(pRecoNode, TPC_VIEW_W));
        purityAdcVectorU.emplace_back(matches.GetPurity(pRecoNode, TPC_VIEW_U, true));
        purityAdcVectorV.emplace_back(matches.GetPurity(pRecoNode, TPC_VIEW_V, true));
        purityAdcVectorW.emplace_back(matches.GetPurity(pRecoNode, TPC_VIEW_W, true));
        completenessAdcVectorU.emplace_back(matches.GetCompleteness(pRecoNode, TPC_VIEW_U, true));
        completenessAdcVectorV.emplace_back(matches.GetCompleteness(pRecoNode, TPC_VIEW_V, true));
        completenessAdcVectorW.emplace_back(matches.GetCompleteness(pRecoNode, TPC_VIEW_W, true));
        if (nMatches == 1)
        {
            // Only makes sense to calculate vertex delta if we have a one-to-one match
            const CartesianVector &recoVertex{LArPfoHelper::GetVertex(matchInfo.GetRecoNeutrino())->GetPosition()};
            vtxDx = recoVertex.GetX() - trueVertex.GetX();
            vtxDy = recoVertex.GetY() - trueVertex.GetY();
            vtxDz = recoVertex.GetZ() - trueVertex.GetZ();
            vtxDr = std::sqrt(vtxDx * vtxDx + vtxDy * vtxDy + vtxDz * vtxDz);
        }

    }

    // Would like to add information on hierarchy matching. Needs some thought, it's extremely complicated

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "event", m_event));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "mcId", mcId));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "mcPDG", pdg));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "mcTier", tier));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "mcNHits", mcHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isNuInteraction", isNeutrinoInt));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isCosmicRay", isCosmicRay));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isTestBeam", isTestBeam));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isLeadingLepton", isLeadingLepton));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isMichel", isMichel));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nMatches", nMatches));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "recoIdVector", &recoIdVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nRecoHitsVector", &nRecoHitsVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nSharedHitsVector", &nSharedHitsVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "purityVector", &purityVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "completenessVector", &completenessVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "purityAdcVector", &purityAdcVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "completenessAdcVector", &completenessAdcVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "purityVectorU", &purityVectorU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "purityVectorV", &purityVectorV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "purityVectorW", &purityVectorW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "completenessVectorU", &completenessVectorU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "completenessVectorV", &completenessVectorV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "completenessVectorW", &completenessVectorW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "purityAdcVectorU", &purityAdcVectorU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "purityAdcVectorV", &purityAdcVectorV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "purityAdcVectorW", &purityAdcVectorW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "completenessAdcVectorU", &completenessAdcVectorU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "completenessAdcVectorV", &completenessAdcVectorV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "completenessAdcVectorW", &completenessAdcVectorW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "vtxDx", vtxDx));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "vtxDy", vtxDy));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "vtxDz", vtxDz));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "vtxDr", vtxDr));
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treename.c_str()));

}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode NDValidationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    if (m_caloHitListName.empty())
        m_caloHitListName = "CaloHitList2D";
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));
    if (m_pfoListName.empty())
        m_pfoListName = "RecreatedPfos";

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ValidateEvent", m_validateEvent));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ValidateMC", m_validateMC));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "WriteTree", m_writeTree));
    if (m_writeTree)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "FileName", m_filename));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TreeName", m_treename));
        if (!(m_validateEvent || m_validateMC))
        {
            std::cout << "Error: WriteTree requested but no tree names found" << std::endl;
            return STATUS_CODE_NOT_FOUND;
        }
        else if (m_validateEvent && m_validateMC)
        {
            std::cout << "Error: Both event-level and MC-level validation requested simulataneously" << std::endl;
            return STATUS_CODE_INVALID_PARAMETER;
        }
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FoldToPrimaries", m_foldToPrimaries));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FoldDynamic", m_foldDynamic));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FoldToLeadingShowers", m_foldToLeadingShowers));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "PrintToScreen", m_printToScreen));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
