/**
 *  @file   include/PrepareClusteringTwoDAlgorithm.h
 *
 *  @brief  Header file for the cluster creation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_PREPARE_CLUSTERING_TWO_D_ALGORITHM_H
#define LAR_PREPARE_CLUSTERING_TWO_D_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include <unordered_map>

namespace lar_content
{

/**
 *  @brief  PrepareClusteringTwoDAlgorithm class
 */
class PrepareClusteringTwoDAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    PrepareClusteringTwoDAlgorithm();

private:
    pandora::StatusCode Run();

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_inputClusterListName;
    std::string m_inputCaloHitListName;
};

} // namespace lar_content

#endif // #ifndef LAR_PREPARE_CLUSTERING_TWO_D_ALGORITHM_H
