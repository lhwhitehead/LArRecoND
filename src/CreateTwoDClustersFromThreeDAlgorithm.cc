/**
 *  @file   src/CreateTwoDClustersFromThreeDAlgorithm.cc
 *
 *  @brief  Implementation of the 3D to 2D cluster creation algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

#include "CreateTwoDClustersFromThreeDAlgorithm.h"

using namespace pandora;

namespace lar_content
{

CreateTwoDClustersFromThreeDAlgorithm::CreateTwoDClustersFromThreeDAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CreateTwoDClustersFromThreeDAlgorithm::Run()
{
    const ClusterList *pClusterList3D{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputClusterListName3D, pClusterList3D));

    const CaloHitList *pCaloHitListU{nullptr};
    const CaloHitList *pCaloHitListV{nullptr};
    const CaloHitList *pCaloHitListW{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputCaloHitListNames2D.at(0), pCaloHitListU));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputCaloHitListNames2D.at(1), pCaloHitListV));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputCaloHitListNames2D.at(2), pCaloHitListW));

    HitUsedMap usedHitsU, usedHitsV, usedHitsW;

    for (const CaloHit *const pCaloHit : *pCaloHitListU)
        usedHitsU[pCaloHit] = false;
    for (const CaloHit *const pCaloHit : *pCaloHitListV)
        usedHitsV[pCaloHit] = false;
    for (const CaloHit *const pCaloHit : *pCaloHitListW)
        usedHitsW[pCaloHit] = false;

    std::vector<PandoraContentApi::Cluster::Parameters> clustersU;
    std::vector<PandoraContentApi::Cluster::Parameters> clustersV;
    std::vector<PandoraContentApi::Cluster::Parameters> clustersW;

    for (const Cluster *pCluster : *pClusterList3D)
    {
        const OrderedCaloHitList &clusterHits3DOrdered = pCluster->GetOrderedCaloHitList();
        CaloHitList clusterHits3D;
        clusterHits3DOrdered.FillCaloHitList(clusterHits3D);

        // Now we need to find the matching 2D hits and build the 2D clusters
        CaloHitList associatedHitsU, associatedHitsV, associatedHitsW;

        for (const CaloHit *pCaloHit3D : clusterHits3D)
        {
            this->GetAssociatedTwoDHit(pCaloHit3D, pCaloHitListU, associatedHitsU, usedHitsU, TPC_VIEW_U);
            this->GetAssociatedTwoDHit(pCaloHit3D, pCaloHitListV, associatedHitsV, usedHitsV, TPC_VIEW_V);
            this->GetAssociatedTwoDHit(pCaloHit3D, pCaloHitListW, associatedHitsW, usedHitsW, TPC_VIEW_W);
        }

        if (!associatedHitsU.empty())
        {
            PandoraContentApi::Cluster::Parameters parametersU;
            parametersU.m_caloHitList = associatedHitsU;
            clustersU.emplace_back(parametersU);
        }

        if (!associatedHitsV.empty())
        {
            PandoraContentApi::Cluster::Parameters parametersV;
            parametersV.m_caloHitList = associatedHitsV;
            clustersV.emplace_back(parametersV);
        }

        if (!associatedHitsW.empty())
        {
            PandoraContentApi::Cluster::Parameters parametersW;
            parametersW.m_caloHitList = associatedHitsW;
            clustersW.emplace_back(parametersW);
        }
    }

    // Create new clusters
    std::string tempClusterListName;
  
    // U View 
    const ClusterList *pClusterListU{nullptr};
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


    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CreateTwoDClustersFromThreeDAlgorithm::GetAssociatedTwoDHit(const CaloHit *const pCaloHit3D, const CaloHitList *const pCaloHitList2D, 
    CaloHitList &associatedHits, HitUsedMap &usedHits2D, const HitType &hitType) const
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
        if (usedHits2D.at(pCaloHit2D) == true)
            continue;

        if (!PandoraContentApi::IsAvailable(*this,pCaloHit2D))
            continue;

        const CartesianVector posTwoD = pCaloHit2D->GetPositionVector();

        if (std::fabs(wirePos - posTwoD.GetZ()) < std::numeric_limits<float>::epsilon() &&
            std::fabs(posThreeD.GetX() - posTwoD.GetX()) < std::numeric_limits<float>::epsilon())
        {
            associatedHits.emplace_back(pCaloHit2D);
            usedHits2D.at(pCaloHit2D) = true;
            break;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CreateTwoDClustersFromThreeDAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{

    std::string tempCaloHitName;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputCaloHitListNameU", tempCaloHitName));
    m_inputCaloHitListNames2D.push_back(tempCaloHitName);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputCaloHitListNameV", tempCaloHitName));
    m_inputCaloHitListNames2D.push_back(tempCaloHitName);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputCaloHitListNameW", tempCaloHitName));
    m_inputCaloHitListNames2D.push_back(tempCaloHitName);

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListName3D", m_inputClusterListName3D));

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
