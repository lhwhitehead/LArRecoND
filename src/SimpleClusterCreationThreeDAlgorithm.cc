/**
 *  @file   src/SimpleClusterCreationThreeDAlgorithm.cc
 *
 *  @brief  Implementation of the cluster creation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

#include "SimpleClusterCreationThreeDAlgorithm.h"

using namespace pandora;

namespace lar_content
{

SimpleClusterCreationThreeDAlgorithm::SimpleClusterCreationThreeDAlgorithm() : m_clusteringWindowSquared(0.25f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SimpleClusterCreationThreeDAlgorithm::Run()
{
    const CaloHitList *pCaloHitList3D{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputCaloHitListName3D, pCaloHitList3D));

    // Build map of associations between selected calo hits
    HitAssociationMap hitAssociationMap;
    this->BuildAssociationMap(pCaloHitList3D, hitAssociationMap);
    if (hitAssociationMap.empty())
        return STATUS_CODE_SUCCESS;

    // Create new clusters
    this->CreateClusters(pCaloHitList3D, hitAssociationMap);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SimpleClusterCreationThreeDAlgorithm::BuildAssociationMap(const CaloHitList *const pCaloHitList, HitAssociationMap &hitAssociationMap) const
{
    for (const CaloHit *const pCaloHitI : *pCaloHitList)
    {
        for (const CaloHit *const pCaloHitJ : *pCaloHitList)
        {
            if (pCaloHitI == pCaloHitJ)
                continue;

            const float distSquared((pCaloHitI->GetPositionVector() - pCaloHitJ->GetPositionVector()).GetMagnitudeSquared());
            if (distSquared < m_clusteringWindowSquared)
            {
                CaloHitList &caloHitListI(hitAssociationMap[pCaloHitI]);

                if (caloHitListI.end() == std::find(caloHitListI.begin(), caloHitListI.end(), pCaloHitJ))
                    caloHitListI.push_back(pCaloHitJ);

                CaloHitList &caloHitListJ(hitAssociationMap[pCaloHitI]);

                if (caloHitListJ.end() == std::find(caloHitListJ.begin(), caloHitListJ.end(), pCaloHitI))
                    caloHitListJ.push_back(pCaloHitI);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SimpleClusterCreationThreeDAlgorithm::CreateClusters(const CaloHitList *const pCaloHitList, const HitAssociationMap &hitAssociationMap) const
{
    CaloHitSet vetoList;
    CaloHitVector caloHitVector(pCaloHitList->begin(), pCaloHitList->end());
    std::sort(caloHitVector.begin(), caloHitVector.end(), LArClusterHelper::SortHitsByPosition);

    std::vector<PandoraContentApi::Cluster::Parameters> clustersU;
    std::vector<PandoraContentApi::Cluster::Parameters> clustersV;
    std::vector<PandoraContentApi::Cluster::Parameters> clustersW;

    CaloHitList usedHitsU, usedHitsV, usedHitsW;

    std::cout << "Making clusters from " << caloHitVector.size() << " 3D hits" << std::endl;
    for (const CaloHit *const pSeedCaloHit : caloHitVector)
    {
        if (vetoList.count(pSeedCaloHit))
            continue;

        CaloHitList mergeList;
        mergeList.emplace_back(pSeedCaloHit);
        this->CollectAssociatedHits(pSeedCaloHit, pSeedCaloHit, hitAssociationMap, vetoList, mergeList);

        // Now we need to find the matching 2D hits and build the 2D clusters
        CaloHitList associatedHitsU;
        this->GetAssociatedTwoDHits(mergeList, associatedHitsU, usedHitsU, m_inputCaloHitListNames2D.at(0), TPC_VIEW_U);
        if (!associatedHitsU.empty())
        {
            PandoraContentApi::Cluster::Parameters parametersU;
            parametersU.m_caloHitList = associatedHitsU;
            clustersU.emplace_back(parametersU);
        }

        CaloHitList associatedHitsV;
        this->GetAssociatedTwoDHits(mergeList, associatedHitsV, usedHitsV, m_inputCaloHitListNames2D.at(1), TPC_VIEW_V);
        if (!associatedHitsV.empty())
        {
            PandoraContentApi::Cluster::Parameters parametersV;
            parametersV.m_caloHitList = associatedHitsV;
            clustersV.emplace_back(parametersV);
        }

        CaloHitList associatedHitsW;
        this->GetAssociatedTwoDHits(mergeList, associatedHitsW, usedHitsW, m_inputCaloHitListNames2D.at(2), TPC_VIEW_W); 
        if (!associatedHitsW.empty())
        {
            PandoraContentApi::Cluster::Parameters parametersW;
            parametersW.m_caloHitList = associatedHitsW;
            clustersW.emplace_back(parametersW);
        }

//        std::cout << "Building 3D cluster with " << mergeList.size() << std::endl;
//        std::cout << "Building 2D clusters: " << associatedHitsU.size() << ", " << associatedHitsV.size() << ", " << associatedHitsW.size() << std::endl;

        vetoList.insert(pSeedCaloHit);
        for (const CaloHit *const pAssociatedCaloHit : mergeList)
        {
            vetoList.insert(pAssociatedCaloHit);
        }
    }
   
    // U View 
    const ClusterList *pClusterListU{nullptr};
    std::string tempClusterListName;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pClusterListU, tempClusterListName));
    for (auto parameters : clustersU)
    {
        const Cluster *pClusterU{nullptr};
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pClusterU));
    }
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Cluster>(*this, m_outputClusterListNames.at(0)));

    // V View
    const ClusterList *pClusterListV{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pClusterListV, tempClusterListName));
    for (auto parameters : clustersV)
    {
        const Cluster *pClusterV{nullptr};
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pClusterV));
    }
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Cluster>(*this, m_outputClusterListNames.at(1)));

    // W View
    const ClusterList *pClusterListW{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pClusterListW, tempClusterListName));
    for (auto parameters : clustersW)
    {
        const Cluster *pClusterW{nullptr};
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pClusterW));
    }
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Cluster>(*this, m_outputClusterListNames.at(2)));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SimpleClusterCreationThreeDAlgorithm::CollectAssociatedHits(const CaloHit *const pSeedCaloHit, const CaloHit *const pCurrentCaloHit,
    const HitAssociationMap &hitAssociationMap, const CaloHitSet &vetoList, CaloHitList &mergeList) const
{
    if (vetoList.count(pCurrentCaloHit))
        return;

    HitAssociationMap::const_iterator iter1 = hitAssociationMap.find(pCurrentCaloHit);
    if (iter1 == hitAssociationMap.end())
        return;

    CaloHitVector caloHitVector(iter1->second.begin(), iter1->second.end());
    std::sort(caloHitVector.begin(), caloHitVector.end(), LArClusterHelper::SortHitsByPosition);

    for (const CaloHit *const pAssociatedCaloHit : caloHitVector)
    {
        if (pAssociatedCaloHit == pSeedCaloHit)
            continue;

        if (mergeList.end() != std::find(mergeList.begin(), mergeList.end(), pAssociatedCaloHit))
            continue;

        mergeList.push_back(pAssociatedCaloHit);

        this->CollectAssociatedHits(pSeedCaloHit, pAssociatedCaloHit, hitAssociationMap, vetoList, mergeList);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SimpleClusterCreationThreeDAlgorithm::GetAssociatedTwoDHits(const CaloHitList &cluster3DHits, CaloHitList &associatedHits,
    CaloHitList &usedHits2D, const std::string &caloHitListName2D, const HitType &hitType) const
{
    const CaloHitList *pCaloHitList2D{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, caloHitListName2D, pCaloHitList2D));

    for (const CaloHit *const pCaloHit3D : cluster3DHits)
    {
        const CartesianVector posThreeD = pCaloHit3D->GetPositionVector();
        float wirePos{0.f};

        if (hitType == TPC_VIEW_U)
            wirePos = PandoraContentApi::GetPlugins(*this)->GetLArTransformationPlugin()->YZtoU(posThreeD.GetY(), posThreeD.GetZ());
        if (hitType == TPC_VIEW_V)
            wirePos = PandoraContentApi::GetPlugins(*this)->GetLArTransformationPlugin()->YZtoV(posThreeD.GetY(), posThreeD.GetZ());
        if (hitType == TPC_VIEW_W)
            wirePos = PandoraContentApi::GetPlugins(*this)->GetLArTransformationPlugin()->YZtoW(posThreeD.GetY(), posThreeD.GetZ());
        
        for (const CaloHit *const pCaloHit2D : *pCaloHitList2D)
        {
            if (!PandoraContentApi::IsAvailable(*this,pCaloHit2D))
                continue;

            if (std::find(usedHits2D.begin(),usedHits2D.end(),pCaloHit2D) != usedHits2D.end())
                continue;

            const CartesianVector posTwoD = pCaloHit2D->GetPositionVector();

            if (std::fabs(wirePos - posTwoD.GetZ()) < std::numeric_limits<float>::epsilon() &&
                std::fabs(posThreeD.GetX() - posTwoD.GetX()) < std::numeric_limits<float>::epsilon())// &&
//                std::fabs(pCaloHit3D->GetInputEnergy() - pCaloHit2D->GetInputEnergy()) < std::numeric_limits<float>::epsilon())
            {
                associatedHits.emplace_back(pCaloHit2D);
                usedHits2D.emplace_back(pCaloHit2D);
                break;
            }
        }
    } 

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SimpleClusterCreationThreeDAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    float clusteringWindow = std::sqrt(m_clusteringWindowSquared);
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ClusteringWindow", clusteringWindow));
    m_clusteringWindowSquared = clusteringWindow * clusteringWindow;

    std::string tempCaloHitName;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputCaloHitListNameU", tempCaloHitName));
    m_inputCaloHitListNames2D.push_back(tempCaloHitName);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputCaloHitListNameV", tempCaloHitName));
    m_inputCaloHitListNames2D.push_back(tempCaloHitName);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputCaloHitListNameW", tempCaloHitName));
    m_inputCaloHitListNames2D.push_back(tempCaloHitName);

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputCaloHitListName3D", m_inputCaloHitListName3D));

    std::string tempClusterName;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputClusterListNameU", tempClusterName));
    m_outputClusterListNames.push_back(tempClusterName);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputClusterListNameV", tempClusterName));
    m_outputClusterListNames.push_back(tempClusterName);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputClusterListNameW", tempClusterName));
    m_outputClusterListNames.push_back(tempClusterName);
    

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
