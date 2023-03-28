/**
 *  @file   include/ReplaceHitAndClusterListsAlgorithm.h
 *
 *  @brief  Header file for the cluster creation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_REPLACE_HIT_AND_CLUSTER_LISTS_ALGORITHM_H
#define LAR_REPLACE_HIT_AND_CLUSTER_LISTS_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include <unordered_map>

namespace lar_content
{

/**
 *  @brief  ReplaceHitAndClusterListsAlgorithm class
 */
class ReplaceHitAndClusterListsAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    ReplaceHitAndClusterListsAlgorithm();

private:
    pandora::StatusCode Run();

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_inputClusterListName;
    std::string m_inputCaloHitListName;
};

} // namespace lar_content

#endif // #ifndef LAR_REPLACE_HIT_AND_CLUSTER_LISTS_ALGORITHM_H
