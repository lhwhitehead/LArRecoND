/**
 *  @file   include/MergeClearTracksThreeDAlgorithm.h
 *
 *  @brief  Header file for the clear 3D track merging algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_MERGE_CLEAR_TRACKS_THREE_D_ALGORITHM_H
#define LAR_MERGE_CLEAR_TRACKS_THREE_D_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"

namespace lar_content
{

/**
 *  @brief  MergeClearTracksThreeDAlgorithm class
 */
class MergeClearTracksThreeDAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    MergeClearTracksThreeDAlgorithm();

    pandora::StatusCode Run();

private:
    typedef std::map<const pandora::Cluster *const, float> ClusterFloatMap;
    typedef std::map<const pandora::Cluster *const, const pandora::Cluster *> ClusterMergeMap;

    bool FindMerges(const pandora::ClusterList *const pClusterList) const;
    void CanMergeClusters(const pandora::Cluster *const pLargeCluster, const pandora::Cluster *const pSmallCluster,
        ClusterMergeMap &mergeCandidates, ClusterFloatMap &mergeDistances) const;
    void MergeClusters(const pandora::Cluster *const pLargeCluster, const pandora::Cluster *const pSmallCluster) const;
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float m_slidingFitWindow;           ///< Number of layers in the sliding fit window
    float m_maxGapLengthCut;            ///< The maximum distance between to tracks to allow merging
    float m_maxGapTransverseCut;        ///< The maximum transverse distance between tracks to allow merging
    float m_minCosThetaCut;             ///< The minimum angle cosine between track segment directions
    std::string m_inputClusterListName; ///< The name of the input 3D cluster list
};

} // namespace lar_content

#endif // #ifndef LAR_MERGE_CLEAR_TRACKS_THREE_D_ALGORITHM_H
