/**
 *  @file   larpandoracontent/LArThreeDReco/LArHitCreation/CompareToInputThreeDHitsTool.cc
 *
 *  @brief  Implementation of the delta ray shower hits tool.
 *
 *  $Log: $
 */

#include "Api/PandoraApi.h"
#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "CompareToInputThreeDHitsTool.h"
#include "larpandoracontent/LArThreeDReco/LArHitCreation/ThreeDHitCreationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

CompareToInputThreeDHitsTool::CompareToInputThreeDHitsTool() : m_inputCaloHitList3DName{""}
{
}

void CompareToInputThreeDHitsTool::Run(ThreeDHitCreationAlgorithm *const pAlgorithm, const ParticleFlowObject *const pPfo,
    const CaloHitVector &inputTwoDHits, ProtoHitVector &protoHitVector)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    std::cout << m_inputCaloHitList3DName << " :: hit creation for pfo " << pPfo << " has " << inputTwoDHits.size() << " 2d hits to play with" << std::endl;
    try
    {
        const CaloHitList *pCaloHits3D{nullptr};
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*pAlgorithm, m_inputCaloHitList3DName, pCaloHits3D));

        std::map<const CaloHit*, bool> isHitUsed;
        for (const CaloHit *const pCaloHit3D : (*pCaloHits3D))
            isHitUsed[pCaloHit3D] = false;
  
        for (const CaloHit *const pCaloHit2D : inputTwoDHits)
        {
            const float xPos2D(pCaloHit2D->GetPositionVector().GetX());
            const float zPos2D(pCaloHit2D->GetPositionVector().GetZ());
            const HitType hitType(pCaloHit2D->GetHitType());

            for (const CaloHit *const pCaloHit3D : (*pCaloHits3D))
            {
                if (isHitUsed.at(pCaloHit3D))
                    continue;

                if (!PandoraContentApi::IsAvailable(*pAlgorithm,pCaloHit3D))
                    continue;

                const CartesianVector pos3D(pCaloHit3D->GetPositionVector());
                const float uPos(PandoraContentApi::GetPlugins(*pAlgorithm)->GetLArTransformationPlugin()->YZtoU(pos3D.GetY(), pos3D.GetZ()));
                const float vPos(PandoraContentApi::GetPlugins(*pAlgorithm)->GetLArTransformationPlugin()->YZtoV(pos3D.GetY(), pos3D.GetZ()));
                const float wPos(PandoraContentApi::GetPlugins(*pAlgorithm)->GetLArTransformationPlugin()->YZtoW(pos3D.GetY(), pos3D.GetZ()));
                const float wirePos(hitType == TPC_VIEW_U ? uPos : (hitType == TPC_VIEW_V ? vPos : wPos));

                if ((std::fabs(xPos2D - pos3D.GetX()) < std::numeric_limits<float>::epsilon())
                    && (std::fabs(zPos2D - wirePos) < std::numeric_limits<float>::epsilon()))
//                    && (std::fabs(pCaloHit2D->GetInputEnergy() - pCaloHit3D->GetInputEnergy()) < std::numeric_limits<float>::epsilon()))
                {
                    ProtoHit protoHit(pCaloHit2D);
                    protoHit.SetPosition3D(pos3D,0); // Second argument is the chi2, set to zero as this is a guaranteed match
                    if (protoHit.IsPositionSet())
                        protoHitVector.push_back(protoHit);
                    isHitUsed[pCaloHit3D] = true;
                    break;
                }
            }   
        }
        unsigned int nUsed{0};
        for (auto const &pair : isHitUsed)
            if (pair.second) ++nUsed;
        std::cout << "Assigned " << nUsed << " 3D hits to Pfo " << pPfo << std::endl;
    }
    catch (StatusCodeException &)
    {
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode CompareToInputThreeDHitsTool::ReadSettings(const pandora::TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "InputCaloHitList3DName", m_inputCaloHitList3DName));
    return HitCreationBaseTool::ReadSettings(xmlHandle); 
}

} // namespace lar_content
