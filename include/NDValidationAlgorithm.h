/**
 *  @file   larpandoracontent/LArMonitoring/NDValidationAlgorithm.h
 *
 *  @brief  Header file for the hierarchy validation algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_ND_VALIDATION_ALGORITHM_H
#define LAR_ND_VALIDATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArHelpers/LArHierarchyHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

namespace lar_content
{

/**
 *  @brief  NDValidationAlgorithm class
 */
class NDValidationAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    NDValidationAlgorithm();

    virtual ~NDValidationAlgorithm();

private:

    class RecoMatchInfo
    {
    public:
        RecoMatchInfo(const pandora::ParticleFlowObject *const pPfo, const int pdg, const int nHits, const int nSharedHits, const float completeness, const float purity, const float completenessADC, const float purityADC);

        void Print() const;

        void SetTwoDValues(const float compU, const float purityU, const float compV, const float purityV, const float compW, const float purityW, const float compADCU, const float purityADCU, const float compADCV, const float purityADCV, const float compADCW, const float purityADCW);

        const pandora::ParticleFlowObject *m_pPfo;
        int m_pdg;
        int m_nHits;
        int m_nSharedHits;
        float m_completeness;
        float m_purity;
        float m_completenessADC;
        float m_purityADC;

        float m_completenessU;
        float m_purityU;
        float m_completenessADCU;
        float m_purityADCU;
        float m_completenessV;
        float m_purityV;
        float m_completenessADCV;
        float m_purityADCV;
        float m_completenessW;
        float m_purityW;
        float m_completenessADCW;
        float m_purityADCW;
    };
    typedef std::vector<RecoMatchInfo> RecoMatchInfoVector;
    typedef std::map<const pandora::MCParticle *, RecoMatchInfoVector> RecoMatchInfoMap;
    typedef std::map<const pandora::MCParticle *, int> MCHitsMap;

    typedef std::vector<LArHierarchyHelper::MCMatches> MCMatchesVector;
    typedef std::map<const pandora::MCParticle *, MCMatchesVector> MCParticleToMCMatchesMap;

    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Extract the reconstructed matches from the MatchInfo object
     *
     *  @param  matchInfo The match info object that stores the hierarchy matching results
     */
    void ExtractRecoMatches(const LArHierarchyHelper::MatchInfo &matchInfo);

    /**
     *  @brief  Extract the reconstructed matches from the MatchInfo object
     *
     *  @param  matchInfo The match info object that stores the hierarchy matching results
     *
     *  @return vector of reco matches from a vector of LArHierarchyHelper::MCMatches objects
     */
    RecoMatchInfoVector CreateMatchInfoObjects(const MCMatchesVector &mcMatches);

    /**
     *  @brief  Collates variables and fills ROOT tree for MC particles with matches
     *
     *  @param  pMCParticle pointer to the MCParticle with matches
     *  @param  matches vector of summarised matched pfos to the MCParticle
     */
    void Fill(const pandora::MCParticle *pMCParticle, const RecoMatchInfoVector &matches);

    /**
     *  @brief  Check if the interaction is CC or NC
     *
     *  @param  pNeutrino the MC neutrino particle
     *
     *  @return true for cc, false for nc
     */
    bool IsCC(const pandora::MCParticle *pNeutrino) const;

    /**
     *  @brief  Get a simple interaction type string from the nuance code
     *
     *  @param  code the nuance code for the interaction
     *
     *  @return either "QE", "MEC", "RES", "DIS", "COH", "NEE", "IMD" or "OTHER"
     */
    std::string ConvertNuanceCodeToString(int code) const;

    int m_event;                   ///< The current event
    std::string m_caloHitListName; ///< Name of input calo hit list
    std::string m_pfoListName;     ///< Name of input PFO list
    bool m_writeTree;              ///< Whether or not to output validation information to a ROOT file
    std::string m_filename;        ///< The name of the ROOT file to write
    std::string m_treename;        ///< The name of the ROOT tree to write
    bool m_foldToPrimaries;        ///< Whether or not to fold the hierarchy back to primary particles
    bool m_foldDynamic;            ///< Whether or not to fold the hierarchy dynamically
    bool m_foldToLeadingShowers;   ///< Whether or not to fold the hierarchy back to leading shower particles
    bool m_printToScreen;          ///< Whether to print information to the terminal
    RecoMatchInfoMap m_matchMap;   ///< Map to consolidate match information across slices
    MCHitsMap m_mcHitsMap;         ///< Map of the number of hits for a given MCParticle
};

} // namespace lar_content

#endif // LAR_ND_VALIDATION_ALGORITHM_H
