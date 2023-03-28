/**
 *  @file   PfoThreeDHitAssignmentAlgorithm.cc
 *
 *  @brief  
 *
 *  $Log: $
 */

#include "Api/PandoraApi.h"
#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "PfoThreeDHitAssignmentAlgorithm.h"

using namespace pandora;

namespace lar_content
{

PfoThreeDHitAssignmentAlgorithm::PfoThreeDHitAssignmentAlgorithm() : m_inputCaloHitList3DName{""}
{
}

pandora::StatusCode PfoThreeDHitAssignmentAlgorithm::Run()
{
    if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    std::cout << "Running 3D hit assignment" << std::endl;

    const CaloHitList *pCaloHits3D{nullptr};
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputCaloHitList3DName, pCaloHits3D));

    CaloHitList availableHits;
    std::map<const CaloHit*, float> availableHitUPos, availableHitVPos, availableHitWPos;
    for (const CaloHit *pCaloHit : (*pCaloHits3D))
    {
        if (!PandoraContentApi::IsAvailable(*this, pCaloHit))
            continue;

        availableHits.emplace_back(pCaloHit);
        const CartesianVector pos3D = pCaloHit->GetPositionVector();
        availableHitUPos[pCaloHit] = PandoraContentApi::GetPlugins(*this)->GetLArTransformationPlugin()->YZtoU(pos3D.GetY(), pos3D.GetZ());
        availableHitVPos[pCaloHit] = PandoraContentApi::GetPlugins(*this)->GetLArTransformationPlugin()->YZtoV(pos3D.GetY(), pos3D.GetZ());
        availableHitWPos[pCaloHit] = PandoraContentApi::GetPlugins(*this)->GetLArTransformationPlugin()->YZtoW(pos3D.GetY(), pos3D.GetZ());
    }

    // Maps to keep track of which 3D hit matches to a 2D hit in a given pfo
    std::map<const ParticleFlowObject*, std::string> pfoToClusterListName;
    std::map<const CaloHit*, const ParticleFlowObject*> availableHitToPfoU, availableHitToPfoV, availableHitToPfoW;
    for (unsigned int i = 0; i < m_inputPfoListNames.size(); ++i)
    {
        const std::string pfoListName(m_inputPfoListNames.at(i));
        const PfoList *pPfoList{nullptr};
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, pfoListName, pPfoList));

        for (const ParticleFlowObject *pPfo : (*pPfoList))
        {
            pfoToClusterListName[pPfo] = m_outputClusterListNames.at(i);

            CaloHitList theseCaloHitsU, theseCaloHitsV, theseCaloHitsW;
            LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_U, theseCaloHitsU);
            LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_V, theseCaloHitsV);
            LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_W, theseCaloHitsW);

            for (const CaloHit *const pHit3D : availableHits)
            {
                const CartesianVector pos3D = pHit3D->GetPositionVector();
                for (const CaloHit *const pHit2D : theseCaloHitsU)
                {
                    const CartesianVector pos2D = pHit2D->GetPositionVector();
                    if (std::fabs(pos3D.GetX() - pos2D.GetX()) > std::numeric_limits<float>::epsilon())
                        continue;
                    if (std::fabs(availableHitUPos.at(pHit3D) - pos2D.GetZ()) > std::numeric_limits<float>::epsilon())
                        continue;
                   availableHitToPfoU[pHit3D] = pPfo;
                   break;
                }

                for (const CaloHit *const pHit2D : theseCaloHitsV)
                {
                    const CartesianVector pos2D = pHit2D->GetPositionVector();
                    if (std::fabs(pos3D.GetX() - pos2D.GetX()) > std::numeric_limits<float>::epsilon())
                        continue;
                    if (std::fabs(availableHitVPos.at(pHit3D) - pos2D.GetZ()) > std::numeric_limits<float>::epsilon())
                        continue;
                   availableHitToPfoV[pHit3D] = pPfo;
                   break;
                }

                for (const CaloHit *const pHit2D : theseCaloHitsW)
                {
                    const CartesianVector pos2D = pHit2D->GetPositionVector();
                    if (std::fabs(pos3D.GetX() - pos2D.GetX()) > std::numeric_limits<float>::epsilon())
                        continue;
                    if (std::fabs(availableHitWPos.at(pHit3D) - pos2D.GetZ()) > std::numeric_limits<float>::epsilon())
                        continue;
                   availableHitToPfoW[pHit3D] = pPfo;
                   break;
                }
            }
        }
    }

    // Assign the hits to the correct Pfos
    std::map<const ParticleFlowObject*, CaloHitList> pfoToHits;
    for (const CaloHit *const pCaloHit : availableHits)
    {
        std::vector<const ParticleFlowObject*> candidatePfos;
        if (availableHitToPfoU.count(pCaloHit))
            candidatePfos.emplace_back(availableHitToPfoU.at(pCaloHit));
        if (availableHitToPfoV.count(pCaloHit))
            candidatePfos.emplace_back(availableHitToPfoV.at(pCaloHit));
        if (availableHitToPfoW.count(pCaloHit))
            candidatePfos.emplace_back(availableHitToPfoW.at(pCaloHit));

        const size_t nPfos = candidatePfos.size();
        if (0 == nPfos)
            continue;

        unsigned int bestPfoIndex{0};
        if (1 == nPfos)
            bestPfoIndex = 0;
        else if (2 == nPfos)
        {
            if (candidatePfos.at(0) == candidatePfos.at(1))
                bestPfoIndex = 0;
            else
            {
                // If the hits are attached to different pfos then add the 3D hit to the biggest one
                const unsigned int nHitsPfo0(LArPfoHelper::GetNumberOfTwoDHits(candidatePfos.at(0)));
                const unsigned int nHitsPfo1(LArPfoHelper::GetNumberOfTwoDHits(candidatePfos.at(1)));
                bestPfoIndex = (nHitsPfo0 >= nHitsPfo1) ? 0 : 1;
            }
        }
        else if (3 == nPfos)
        {
            if (candidatePfos.at(0) == candidatePfos.at(1) && candidatePfos.at(0) == candidatePfos.at(2))
                bestPfoIndex = 0;
            else if (candidatePfos.at(0) != candidatePfos.at(1) && candidatePfos.at(0) != candidatePfos.at(1))
            {
                const unsigned int nHitsPfo0(LArPfoHelper::GetNumberOfTwoDHits(candidatePfos.at(0)));
                const unsigned int nHitsPfo1(LArPfoHelper::GetNumberOfTwoDHits(candidatePfos.at(1)));
                const unsigned int nHitsPfo2(LArPfoHelper::GetNumberOfTwoDHits(candidatePfos.at(2)));
                if (nHitsPfo0 >= nHitsPfo1 && nHitsPfo0 >= nHitsPfo2)
                    bestPfoIndex = 0;
                else if (nHitsPfo1 >= nHitsPfo0 && nHitsPfo1 >= nHitsPfo2)
                    bestPfoIndex = 1;
                else
                    bestPfoIndex = 2;   
            }
            else
            {
                // We are in the best two of three situation - pick the pfo with two hits
                if (candidatePfos.at(0) == candidatePfos.at(1) || candidatePfos.at(0) == candidatePfos.at(2))
                    bestPfoIndex = 0;
                else
                    bestPfoIndex = 1;
            }
        }
        
        if (!pfoToHits.count(candidatePfos.at(bestPfoIndex)))
        {
            pfoToHits[candidatePfos.at(bestPfoIndex)] = CaloHitList();
        }
        pfoToHits[candidatePfos.at(bestPfoIndex)].emplace_back(pCaloHit);
    }

    // Assign the hits
    for (auto const &pfoPair : pfoToHits)
        AddHitsToPfo(pfoPair.first,pfoPair.second,pfoToClusterListName.at(pfoPair.first));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PfoThreeDHitAssignmentAlgorithm::AddHitsToPfo(const ParticleFlowObject *pPfo, const CaloHitList &hits, const std::string listName) const
{
    ClusterList clusters3D;
    LArPfoHelper::GetThreeDClusterList(pPfo,clusters3D);
    const unsigned int nClusters(clusters3D.size());

    if (0 == nClusters)
    {
        const ClusterList *pClusterList(nullptr);
        std::string clusterListName;
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pClusterList, clusterListName));

        PandoraContentApi::Cluster::Parameters clusterParams;
        clusterParams.m_caloHitList.insert(clusterParams.m_caloHitList.end(), hits.begin(), hits.end());

        const Cluster *pCluster3D(nullptr);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, clusterParams, pCluster3D));
        if (!pCluster3D || !pClusterList || pClusterList->empty())
            throw StatusCodeException(STATUS_CODE_FAILURE);

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Cluster>(*this, listName));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToPfo(*this, pPfo, pCluster3D));
        std::cout << "Pfo " << pPfo << ": created cluster " << pCluster3D << " with " << hits.size() << " hits" << std::endl;
    }
    else if (1 == nClusters)
    {
        const Cluster *const pCluster3D = *(clusters3D.begin());
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, pCluster3D, &hits));
        std::cout << "Pfo " << pPfo << ": added " << hits.size() << " hits to existing cluster " << pCluster3D << std::endl;
    }
    else
        throw StatusCodeException(STATUS_CODE_FAILURE);

}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode PfoThreeDHitAssignmentAlgorithm::ReadSettings(const pandora::TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "InputCaloHitList3DName", m_inputCaloHitList3DName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "InputPfoListNames", m_inputPfoListNames));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "OutputClusterListNames", m_outputClusterListNames));

    if (m_inputPfoListNames.size() != m_outputClusterListNames.size())
    {
        std::cout << "LArPfoThreeDHitAssignment: number of input pfo list names must be the same as the number of output cluster list names" << std::endl;
        return STATUS_CODE_FAILURE;
    }

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
