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

NDValidationAlgorithm::RecoMatchInfo::RecoMatchInfo(const ParticleFlowObject *const pPfo, const int pdg, const int nHits, const int nSharedHits, const float completeness, const float purity, const float completenessADC, const float purityADC):
    m_pPfo{pPfo},
    m_pdg{pdg},
    m_nHits{nHits},
    m_nSharedHits{nSharedHits},
    m_completeness{completeness}, m_purity{purity}, m_completenessADC{completenessADC}, m_purityADC{purityADC},
    m_completenessU{-1.f}, m_purityU{-1.f}, m_completenessADCU{-1.f}, m_purityADCU{-1.f},
    m_completenessV{-1.f}, m_purityV{-1.f}, m_completenessADCV{-1.f}, m_purityADCV{-1.f},
    m_completenessW{-1.f}, m_purityW{-1.f}, m_completenessADCW{-1.f}, m_purityADCW{-1.f}
{}

//------------------------------------------------------------------------------------------------------------------------------------------

void NDValidationAlgorithm::RecoMatchInfo::Print() const
{
    std::cout << "    - Pfo: pdg = " << m_pdg << ", hits = " << m_nHits << ", shared hits = " << m_nSharedHits << ", completeness = " << m_completeness << ", purity = " << m_purity << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NDValidationAlgorithm::RecoMatchInfo::SetTwoDValues(const float compU, const float purityU, const float compV, const float purityV, const float compW, const float purityW, const float compADCU, const float purityADCU, const float compADCV, const float purityADCV, const float compADCW, const float purityADCW)
{
    m_completenessU = compU; m_purityU = purityU; m_completenessADCU = compADCU; m_purityADCU = purityADCU;
    m_completenessV = compV; m_purityV = purityV; m_completenessADCV = compADCV; m_purityADCV = purityADCV;
    m_completenessW = compW; m_purityW = purityW; m_completenessADCW = compADCW; m_purityADCW = purityADCW;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

NDValidationAlgorithm::NDValidationAlgorithm() :
    m_event{-1},
    m_writeTree{false},
    m_foldToPrimaries{false},
    m_foldDynamic{false},
    m_foldToLeadingShowers{false},
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
        }
    }
    if(m_printToScreen)
        std::cout << "Found " << trueNeutrinoMap.size() << " true neutrinos" << std::endl;

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
    if(m_printToScreen)
        std::cout << "Found " << recoNeutrinoMap.size() << " reconstructed neutrinos" << std::endl;

    LArHierarchyHelper::FoldingParameters foldParameters;
    if (m_foldToPrimaries)
        foldParameters.m_foldToTier = true;
    else if (m_foldDynamic)
        foldParameters.m_foldDynamic = true;
    else if (m_foldToLeadingShowers)
        foldParameters.m_foldToLeadingShowers = true;

    for (auto const &truePair : trueNeutrinoMap)
    {
        LArHierarchyHelper::MCHierarchy mcHierarchy;
        LArHierarchyHelper::FillMCHierarchy(truePair.second, *pCaloHitList, foldParameters, mcHierarchy);
        for (auto const &recoPair : recoNeutrinoMap)
        {
            LArHierarchyHelper::RecoHierarchy recoHierarchy;
            LArHierarchyHelper::FillRecoHierarchy(recoPair.second, foldParameters, recoHierarchy);
            LArHierarchyHelper::MatchInfo matchInfo;
            LArHierarchyHelper::MatchHierarchies(mcHierarchy, recoHierarchy, matchInfo);
            this->ExtractRecoMatches(matchInfo);
        }
    }

    if (m_printToScreen)
    {
        for (auto const &mcNeutrinoList : trueNeutrinoMap)
        {
            const MCParticle *const pNu = mcNeutrinoList.first;
            const int nuId = static_cast<int>(reinterpret_cast<std::intptr_t>(pNu->GetUid()));
            const bool isCC = this->IsCC(pNu);
            const int nuanceCode{(dynamic_cast<const LArMCParticle *>(pNu))->GetNuanceCode()};
            std::cout << std::endl;
            std::cout << "===== Matching Information for Neutrino: id = " << nuId << ", pdg = " << pNu->GetParticleId() << ", interaction = " << (isCC ? "CC " : "NC ") << this->ConvertNuanceCodeToString(nuanceCode) << " =====" << std::endl;
            const MCParticleList &nuDescendents = mcNeutrinoList.second;
            unsigned int nReconstructableChildren{0};
            for (auto const &mcMatchPair : m_matchMap)
            { 
                const MCParticle *const pMC{mcMatchPair.first};
                if (std::find(nuDescendents.begin(),nuDescendents.end(),pMC) == nuDescendents.end())
                    continue;

                ++nReconstructableChildren;
                const int mcId{static_cast<int>(reinterpret_cast<std::intptr_t>(pMC->GetUid()))};
                const int pdg{pMC->GetParticleId()};
                const int mcHits{m_mcHitsMap.at(pMC)};
                std::cout << "  - MC Particle: id = " << mcId << ", pdg = " << pdg << ", nhits = " << mcHits << ", matches " << mcMatchPair.second.size() << ":" << std::endl;
                for (const RecoMatchInfo &match : mcMatchPair.second)
                    match.Print();
            }
            if (!nReconstructableChildren)
                std::cout << "- No true particles to reconstruct in the detector volume" << std::endl;
        }
    }

    if (m_writeTree)
    {
        for (auto const &match : m_matchMap)
            this->Fill(match.first, match.second);
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NDValidationAlgorithm::ExtractRecoMatches(const LArHierarchyHelper::MatchInfo &match)
{
    MCParticleToMCMatchesMap mcParticleMatches;

    // Matches for different primary particles
    for (const LArHierarchyHelper::MCMatches &mcMatches : match.GetMatches())
    {
        const LArHierarchyHelper::MCHierarchy::Node *pMCNode{mcMatches.GetMC()};
        const MCParticle *pMCParticle{pMCNode->GetLeadingMCParticle()};

        if (!mcParticleMatches.count(pMCParticle))
            mcParticleMatches[pMCParticle] = MCMatchesVector();
        mcParticleMatches[pMCParticle].emplace_back(mcMatches);

        // Set the number of hits if we haven't before
        if (!m_mcHitsMap.count(pMCParticle))
            m_mcHitsMap[pMCParticle] = static_cast<int>(pMCNode->GetCaloHits().size());
    }

    for (auto const &pair : mcParticleMatches)
    {
        RecoMatchInfoVector newMatches = this->CreateMatchInfoObjects(pair.second);
        m_matchMap[pair.first].insert(m_matchMap[pair.first].end(), newMatches.begin(), newMatches.end());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

NDValidationAlgorithm::RecoMatchInfoVector NDValidationAlgorithm::CreateMatchInfoObjects(const MCMatchesVector &mcMatches)
{
    RecoMatchInfoVector outputMatchVector;

    for (const LArHierarchyHelper::MCMatches &match : mcMatches)
    {
        for (const LArHierarchyHelper::RecoHierarchy::Node *pRecoNode : match.GetRecoMatches())
        {
            const int recoId{pRecoNode->GetParticleId()};
            const int nRecoHits{static_cast<int>(pRecoNode->GetCaloHits().size())};
            const int nSharedHits{static_cast<int>(match.GetSharedHits(pRecoNode))};
            const float completeness{match.GetCompleteness(pRecoNode)};
            const float purity{match.GetPurity(pRecoNode)};
            const float completenessADC{match.GetCompleteness(pRecoNode, true)};
            const float purityADC{match.GetPurity(pRecoNode, true)};

            const ParticleFlowObject *const pPfo{pRecoNode->GetRecoParticles().front()};
            RecoMatchInfo newRecoMatch(pPfo, recoId, nRecoHits, nSharedHits, completeness, purity, completenessADC, purityADC);
            
            const float compU{match.GetCompleteness(pRecoNode, TPC_VIEW_U)};
            const float purityU{match.GetPurity(pRecoNode, TPC_VIEW_U)};
            const float compADCU{match.GetCompleteness(pRecoNode, TPC_VIEW_U, true)};
            const float purityADCU{match.GetPurity(pRecoNode, TPC_VIEW_U, true)};

            const float compV{match.GetCompleteness(pRecoNode, TPC_VIEW_V)};
            const float purityV{match.GetPurity(pRecoNode, TPC_VIEW_V)};
            const float compADCV{match.GetCompleteness(pRecoNode, TPC_VIEW_V, true)};
            const float purityADCV{match.GetPurity(pRecoNode, TPC_VIEW_V, true)};

            const float compW{match.GetCompleteness(pRecoNode, TPC_VIEW_W)};
            const float purityW{match.GetPurity(pRecoNode, TPC_VIEW_W)};
            const float compADCW{match.GetCompleteness(pRecoNode, TPC_VIEW_W, true)};
            const float purityADCW{match.GetPurity(pRecoNode, TPC_VIEW_W, true)};

            newRecoMatch.SetTwoDValues(compU, purityU, compV, purityV, compW, purityW, compADCU, purityADCU, compADCV, purityADCV, compADCW, purityADCW);
  
            outputMatchVector.emplace_back(newRecoMatch);
        }
    }

    return outputMatchVector;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NDValidationAlgorithm::Fill(const MCParticle *pMCParticle, const RecoMatchInfoVector &matches)
{
    const int mcId{static_cast<int>(reinterpret_cast<std::intptr_t>(pMCParticle->GetUid()))};

    const MCParticle *const pMCParent = LArMCParticleHelper::GetParentMCParticle(pMCParticle);
    if (pMCParent == nullptr)
    {
        std::cout << "Warning: particle has no neutrino parent. Not filling the tree" << std::endl;
        return;
    }
    const int nuId{static_cast<int>(reinterpret_cast<std::intptr_t>(pMCParent->GetUid()))};

    // All MC nodes contain the same true particle information, so we can just use the first one
    const int isTestBeam{LArMCParticleHelper::IsBeamParticle(pMCParticle) ? 1 : 0};
    const int isCosmicRay{!isTestBeam && LArMCParticleHelper::IsCosmicRay(pMCParticle) ? 1 : 0};
    const int isNeutrinoInt{!(isTestBeam || isCosmicRay) ? 1 : 0};
    const int pdg{pMCParticle->GetParticleId()};
    const int tier{LArMCParticleHelper::GetHierarchyTier(pMCParticle)};
    const int mcHits{m_mcHitsMap.at(pMCParticle)};
    const int isLeadingLepton{LArMCParticleHelper::IsLeading(pMCParticle) ? 1 : 0};

    const MCParticleList &parentList{pMCParticle->GetParentList()};
    const int isElectron{std::abs(pMCParticle->GetParticleId()) == E_MINUS ? 1 : 0};
    const int hasMuonParent{parentList.size() == 1 && std::abs(parentList.front()->GetParticleId()) == MU_MINUS ? 1 : 0};
    const int isMichel{isElectron && hasMuonParent && LArMCParticleHelper::IsDecay(pMCParticle) ? 1 : 0};

    const float mcVtxX{pMCParticle->GetVertex().GetX()};
    const float mcVtxY{pMCParticle->GetVertex().GetY()};
    const float mcVtxZ{pMCParticle->GetVertex().GetZ()};

    IntVector recoIdVector, nRecoHitsVector, nSharedHitsVector;
    FloatVector recoVtxXVector, recoVtxYVector, recoVtxZVector;
    FloatVector purityVector, completenessVector;
    FloatVector purityAdcVector, completenessAdcVector;
    FloatVector purityVectorU, purityVectorV, purityVectorW, completenessVectorU, completenessVectorV, completenessVectorW;
    FloatVector purityAdcVectorU, purityAdcVectorV, purityAdcVectorW, completenessAdcVectorU, completenessAdcVectorV, completenessAdcVectorW;

    for (const RecoMatchInfo &match : matches)
    {
        recoIdVector.emplace_back(match.m_pdg);
        nRecoHitsVector.emplace_back(match.m_nHits);
        nSharedHitsVector.emplace_back(match.m_nSharedHits);
        purityVector.emplace_back(match.m_completeness);
        completenessVector.emplace_back(match.m_purity);
        purityAdcVector.emplace_back(match.m_completenessADC);
        completenessAdcVector.emplace_back(match.m_purityADC);

        // Per view information
        purityVectorU.emplace_back(match.m_purityU);
        purityVectorV.emplace_back(match.m_purityV);
        purityVectorW.emplace_back(match.m_purityW);
        completenessVectorU.emplace_back(match.m_completenessU);
        completenessVectorV.emplace_back(match.m_completenessV);
        completenessVectorW.emplace_back(match.m_completenessW);
        purityAdcVectorU.emplace_back(match.m_purityADCU);
        purityAdcVectorV.emplace_back(match.m_purityADCV);
        purityAdcVectorW.emplace_back(match.m_purityADCW);
        completenessAdcVectorU.emplace_back(match.m_completenessADCU);
        completenessAdcVectorV.emplace_back(match.m_completenessADCV);
        completenessAdcVectorW.emplace_back(match.m_completenessADCW);

        // Vertex information
        VertexList vertexList = match.m_pPfo->GetVertexList();
        if (vertexList.empty())
        {
            recoVtxXVector.emplace_back(-999.f);
            recoVtxYVector.emplace_back(-999.f);
            recoVtxZVector.emplace_back(-999.f);
        }
        else
        {
            vertexList.sort([](const Vertex *a, const Vertex *b){return a->GetPosition().GetZ() < b->GetPosition().GetZ();});
            recoVtxXVector.emplace_back(vertexList.front()->GetPosition().GetX());
            recoVtxYVector.emplace_back(vertexList.front()->GetPosition().GetY());
            recoVtxZVector.emplace_back(vertexList.front()->GetPosition().GetZ());
        }
    }
    const int nMatches{static_cast<int>(recoIdVector.size())};

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "event", m_event));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "mcId", mcId));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nuId", nuId));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "mcPDG", pdg));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "mcTier", tier));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "mcNHits", mcHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "mcVertexX", mcVtxX));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "mcVertexY", mcVtxY));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "mcVertexZ", mcVtxZ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isNuInteraction", isNeutrinoInt));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isCosmicRay", isCosmicRay));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isTestBeam", isTestBeam));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isLeadingLepton", isLeadingLepton));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isMichel", isMichel));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nMatches", nMatches));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "recoIdVector", &recoIdVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "recoVertexXVector", &recoVtxXVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "recoVertexYVector", &recoVtxYVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "recoVertexZVector", &recoVtxZVector));
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
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treename.c_str()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool NDValidationAlgorithm::IsCC(const MCParticle *pNeutrino) const
{
    bool isCC{true};
    for (auto const *child : pNeutrino->GetDaughterList())
    {
        if (child->GetParticleId() == pNeutrino->GetParticleId())
        {
            isCC = false;
            break;
        }
    }
    return isCC;
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::string NDValidationAlgorithm::ConvertNuanceCodeToString(int code) const
{
    if (code == 1001 || code == 1002)
        return "QEL";
    else if (code == 10)
        return "MEC";
    else if (code == 1)
        return "RES";
    else if (code == 1091 || code == 1092)
        return "DIS";
    else if (code == 3 || code == 4)
        return "COH";
    else if (code == 1098)
        return "NEE"; // Neutrino - electron elastic scattering
    else if (code == 1099)
        return "IMD";
    else
        return "OTHER";
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

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "WriteTree", m_writeTree));
    if (m_writeTree)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "FileName", m_filename));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TreeName", m_treename));
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FoldToPrimaries", m_foldToPrimaries));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FoldDynamic", m_foldDynamic));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FoldToLeadingShowers", m_foldToLeadingShowers));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "PrintToScreen", m_printToScreen));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
