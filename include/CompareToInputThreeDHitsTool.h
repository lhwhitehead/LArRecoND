/**
 *  @file   larpandoracontent/LArThreeDReco/LArHitCreation/CompareToInputThreeDHitsTool.h
 *
 *  @brief  Header file for the 3d hit creation tool using input 3D hits
 *
 *  $Log: $
 */
#ifndef COMPARE_TO_INPUT_THREE_D_HITS_TOOL_H
#define COMPARE_TO_INPUT_THREE_D_HITS_TOOL_H 1

#include "larpandoracontent/LArThreeDReco/LArHitCreation/HitCreationBaseTool.h"

namespace lar_content
{

/**
 *  @brief  CompareToInputThreeDHitsTool class
 */
class CompareToInputThreeDHitsTool : public HitCreationBaseTool
{
public:
    CompareToInputThreeDHitsTool();

    virtual void Run(ThreeDHitCreationAlgorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pPfo,
        const pandora::CaloHitVector &inputTwoDHits, ProtoHitVector &protoHitVector);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_inputCaloHitList3DName;
};

} // namespace lar_content

#endif // #ifndef COMPARE_TO_INPUT_THREE_D_HITS_TOOL_H
