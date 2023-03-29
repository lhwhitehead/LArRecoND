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
        const CartesianVector positionI = pCaloHitI->GetPositionVector();
        for (const CaloHit *const pCaloHitJ : *pCaloHitList)
        {
            if (pCaloHitI == pCaloHitJ)
                continue;

            const float distSquared((positionI - pCaloHitJ->GetPositionVector()).GetMagnitudeSquared());
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

    std::vector<PandoraContentApi::Cluster::Parameters> clusters3D;

    std::cout << "Making clusters from " << caloHitVector.size() << " 3D hits" << std::endl;
    for (const CaloHit *const pSeedCaloHit : caloHitVector)
    {
        if (vetoList.count(pSeedCaloHit))
            continue;

        CaloHitList mergeList;
        mergeList.emplace_back(pSeedCaloHit);
        this->CollectAssociatedHits(pSeedCaloHit, pSeedCaloHit, hitAssociationMap, vetoList, mergeList);

        if (!mergeList.empty())
        {
            PandoraContentApi::Cluster::Parameters parameters3D;
            parameters3D.m_caloHitList = mergeList;
            clusters3D.emplace_back(parameters3D);
        }

        vetoList.insert(pSeedCaloHit);
        for (const CaloHit *const pAssociatedCaloHit : mergeList)
        {
            vetoList.insert(pAssociatedCaloHit);
        }
    }
 
    // 3D clusters
    const ClusterList *pClusterList3D{nullptr};
    std::string tempClusterListName;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pClusterList3D, tempClusterListName));
    for (auto parameters : clusters3D)
    {
        const Cluster *pCluster3D{nullptr};
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pCluster3D));
    }
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Cluster>(*this, m_outputClusterListName3D));

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

StatusCode SimpleClusterCreationThreeDAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    float clusteringWindow = std::sqrt(m_clusteringWindowSquared);
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ClusteringWindow", clusteringWindow));
    m_clusteringWindowSquared = clusteringWindow * clusteringWindow;

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputCaloHitListName3D", m_inputCaloHitListName3D));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputClusterListName3D", m_outputClusterListName3D));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
