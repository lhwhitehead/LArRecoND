/**
 *  @file   include/PfoThreeDHitAssignmentAlgorithm.h
 *
 *  @brief  Header file for the cluster creation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_PFO_THREE_D_HIT_ASSIGNMENT_ALGORITHM_H
#define LAR_PFO_THREE_D_HIT_ASSIGNMENT_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Objects/ParticleFlowObject.h"
#include "Objects/CaloHit.h"

namespace lar_content
{

/**
 *  @brief  PfoThreeDHitAssignmentAlgorithm class
 */
class PfoThreeDHitAssignmentAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    PfoThreeDHitAssignmentAlgorithm();

private:
    pandora::StatusCode Run();

    void AddHitsToPfo(const pandora::ParticleFlowObject *pPfo, const pandora::CaloHitList &hits, const std::string listName) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_inputCaloHitList3DName;
    std::vector<std::string> m_inputPfoListNames;
    std::vector<std::string> m_outputClusterListNames;
};

} // namespace lar_content

#endif // #ifndef LAR_PFO_THREE_D_HIT_ASSIGNMENT_ALGORITHM_H
