/**
 *  @file   larpandoracontent/LArNDContent.cc
 *
 *  @brief  Factory implementations for content intended for use with particle flow reconstruction at liquid argon time projection chambers
 *
 *  $Log: $
 */

#include "Api/PandoraApi.h"

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmTool.h"
#include "Pandora/Pandora.h"

#include "MasterThreeDAlgorithm.h"
#include "SlicingThreeDAlgorithm.h"
#include "PrepareClusteringTwoDAlgorithm.h" 
#include "PreProcessingThreeDAlgorithm.h"
#include "SimpleClusterCreationThreeDAlgorithm.h"

#include "EventSlicingThreeDTool.h"
#include "CompareToInputThreeDHitsTool.h"
#include "LArNDContent.h"

// clang-format off
#define LAR_ND_ALGORITHM_LIST(d)                                                                                                   \
    d("LArMasterThreeD",                        MasterThreeDAlgorithm)                                                             \
    d("LArSimpleClusterCreationThreeD",         SimpleClusterCreationThreeDAlgorithm)                                              \
    d("LArSlicingThreeD",                       SlicingThreeDAlgorithm)                                                            \
    d("LArPrepareClusteringTwoD",               PrepareClusteringTwoDAlgorithm)                                                    \
    d("LArPreProcessingThreeD",                 PreProcessingThreeDAlgorithm)                                                      \

#define LAR_ND_ALGORITHM_TOOL_LIST(d)                                                                                              \
    d("LArEventSlicingThreeD",                  EventSlicingThreeDTool)                                                            \
    d("LArCompareToInputThreeDHits",            CompareToInputThreeDHitsTool)                                                      \

#define FACTORY Factory

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_content
{

#define LAR_ND_CONTENT_CREATE_ALGORITHM_FACTORY(a, b)                                                                              \
class b##FACTORY : public pandora::AlgorithmFactory                                                                             \
{                                                                                                                               \
public:                                                                                                                         \
    pandora::Algorithm *CreateAlgorithm() const {return new b;};                                                                \
};

LAR_ND_ALGORITHM_LIST(LAR_ND_CONTENT_CREATE_ALGORITHM_FACTORY)

//------------------------------------------------------------------------------------------------------------------------------------------

#define LAR_ND_CONTENT_CREATE_ALGORITHM_TOOL_FACTORY(a, b)                                                                         \
class b##FACTORY : public pandora::AlgorithmToolFactory                                                                         \
{                                                                                                                               \
public:                                                                                                                         \
    pandora::AlgorithmTool *CreateAlgorithmTool() const {return new b;};                                                        \
};

LAR_ND_ALGORITHM_TOOL_LIST(LAR_ND_CONTENT_CREATE_ALGORITHM_TOOL_FACTORY)

} // namespace lar_content

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

#define LAR_ND_CONTENT_REGISTER_ALGORITHM(a, b)                                                                                    \
{                                                                                                                               \
    const pandora::StatusCode statusCode(PandoraApi::RegisterAlgorithmFactory(pandora, a, new lar_content::b##FACTORY));        \
    if (pandora::STATUS_CODE_SUCCESS != statusCode)                                                                             \
        return statusCode;                                                                                                      \
}

#define LAR_ND_CONTENT_REGISTER_ALGORITHM_TOOL(a, b)                                                                               \
{                                                                                                                               \
    const pandora::StatusCode statusCode(PandoraApi::RegisterAlgorithmToolFactory(pandora, a, new lar_content::b##FACTORY));    \
    if (pandora::STATUS_CODE_SUCCESS != statusCode)                                                                             \
        return statusCode;                                                                                                      \
}

pandora::StatusCode LArNDContent::RegisterAlgorithms(const pandora::Pandora &pandora)
{
    LAR_ND_ALGORITHM_LIST(LAR_ND_CONTENT_REGISTER_ALGORITHM);
    LAR_ND_ALGORITHM_TOOL_LIST(LAR_ND_CONTENT_REGISTER_ALGORITHM_TOOL);
    return pandora::STATUS_CODE_SUCCESS;
}

// clang-format on
