/**
 *  @file   include/SimpleClusterCreationThreeDAlgorithm.h
 *
 *  @brief  Header file for the cluster creation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_SIMPLE_CLUSTER_CREATION_THREE_D_ALGORITHM_H
#define LAR_SIMPLE_CLUSTER_CREATION_THREE_D_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include <unordered_map>

namespace lar_content
{

/**
 *  @brief  SimpleClusterCreationThreeDAlgorithm class
 */
class SimpleClusterCreationThreeDAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    SimpleClusterCreationThreeDAlgorithm();

private:
    pandora::StatusCode Run();

    typedef std::unordered_map<const pandora::CaloHit *, pandora::CaloHitList> HitAssociationMap;

    /**
     *  @brief Create map of associations between calo hits
     *
     *  @param caloHitList The pointer to the input list of calo hits
     *  @param hitAssociationMap The map of associations between calo hits
     */
    void BuildAssociationMap(const pandora::CaloHitList *const pCaloHitList, HitAssociationMap &hitAssociationMap) const;

    /**
     *  @brief Create clusters from selected calo hits and their associations
     *
     *  @param caloHitList The input list of calo hits
     *  @param hitAssociationMap The map of associations between calo hits
     *
     *  @return status code
     */
    pandora::StatusCode CreateClusters(const pandora::CaloHitList *const pCaloHitList, const HitAssociationMap &hitAssociationMap) const;

    /**
     *  @brief For a given seed calo hits, collect up all the associated calo hits
     *
     *  @param pSeedCaloHit the seed calo hits
     *  @param pCurrentCaloHit a possible associated calo hit
     *  @param hitAssociationMap the map of associations between hits
     *  @param vetoList the list of used calo hits
     *  @param mergeList the list of hits associated with the seed hit
     */
    void CollectAssociatedHits(const pandora::CaloHit *const pSeedCaloHit, const pandora::CaloHit *const pCurrentCaloHit,
        const HitAssociationMap &hitAssociationMap, const pandora::CaloHitSet &vetoList, pandora::CaloHitList &mergeList) const;

    /**
     *  @brief Get the hits in each view matching to the clustered 3D hits
     *
     *  @param pCaloHit3D a pointer to the 3D hit
     *  @param pCaloHitList2D a pointer to the 2D hit list
     *  @param associatedHits reference to empty hit list filled with 2D hits
     *  @param usedHits2D list of already used 2D hits
     *  @param hitType the type of hits in the hit list
     */  
    void GetAssociatedTwoDHit(const pandora::CaloHit *const pCaloHit3D, const pandora::CaloHitList *const pCaloHitList2D,
        pandora::CaloHitList &associatedHits, std::map<const pandora::CaloHit*, bool> &usedHits2D, const pandora::HitType &hitType) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float m_clusteringWindowSquared; ///< Maximum distance (squared) for two hits to be joined

    std::vector<std::string> m_inputCaloHitListNames2D;
    std::string              m_inputCaloHitListName3D;
    std::vector<std::string> m_outputClusterListNames;
};

} // namespace lar_content

#endif // #ifndef LAR_SIMPLE_CLUSTER_CREATION_THREE_D_ALGORITHM_H
