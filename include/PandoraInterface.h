/**
 *  @file   LArReco/include/PandoraInterface.h
 *
 *  @brief  Header file for PandoraInterface.
 *
 *  $Log: $
 */
#ifndef PANDORA_ND_INTERFACE_H
#define PANDORA_ND_INTERFACE_H 1

#include "Pandora/PandoraInputTypes.h"
#include "TG4Event.h"
#include "TGeoManager.h"
#include "TGeoNode.h"

#include "LArGrid.h"
#include "LArHitInfo.h"
#include "LArNDGeomSimple.h"
#include "LArSED.h"
#include "LArVoxel.h"

namespace pandora
{
class Pandora;
}

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_nd_reco
{

typedef std::map<int, float> MCParticleEnergyMap;
typedef std::vector<LArVoxel> LArVoxelList;

/**
 *  @brief  Parameters class
 */
class Parameters
{
public:
    /**
     *  @brief  Default constructor
     */
    Parameters();

    enum LArNDFormat
    {
        EDepSim = 0,
        SED = 1,
        SP = 2
    };

    LArNDFormat m_dataFormat; ///< The expected input data format (EDepSim rooTracker or SED ROOT)

    std::string m_settingsFile;  ///< The path to the pandora settings file
                                 ///< (mandatory parameter)
    std::string m_inputFileName; ///< The path to the input file containing events
                                 ///< and/or geometry information
    std::string m_inputTreeName; ///< The optional name of the event TTree

    std::string m_geomFileName;    ///< The ROOT file name containing the TGeoManager info
    std::string m_geomManagerName; ///< The name of the TGeoManager

    std::string m_geometryVolName;  ///< The name of the Geant4 detector placement volume
    std::string m_sensitiveDetName; ///< The name of the Geant4 sensitive hit detector
    bool m_useModularGeometry;      ///< Include each TPC as a separate volume in the geometry

    int m_nEventsToProcess;          ///< The number of events to process (default all
                                     ///< events in file)
    bool m_shouldDisplayEventNumber; ///< Whether event numbers should be
                                     ///< displayed (default false)

    bool m_shouldRunAllHitsCosmicReco;  ///< Whether to run all hits cosmic-ray reconstruction
    bool m_shouldRunStitching;          ///< Whether to stitch cosmic-ray muons crossing between volumes
    bool m_shouldRunCosmicHitRemoval;   ///< Whether to remove hits from tagged cosmic-rays
    bool m_shouldRunSlicing;            ///< Whether to slice events into separate regions for processing
    bool m_shouldRunNeutrinoRecoOption; ///< Whether to run neutrino reconstruction for each slice
    bool m_shouldRunCosmicRecoOption;   ///< Whether to run cosmic-ray reconstruction for each slice
    bool m_shouldPerformSliceId;        ///< Whether to identify slices and select most appropriate pfos
    bool m_printOverallRecoStatus;      ///< Whether to print current operation status messages

    int m_nEventsToSkip;       ///< The number of events to skip
    int m_maxMergedVoxels;     ///< The max number of merged voxels to process (default all)
    float m_minVoxelMipEquivE; ///< The minimum required voxel equivalent MIP energy (default = 0.3)

    bool m_use3D;     ///< Create 3D LArCaloHits
    bool m_useLArTPC; ///< Create LArTPC LArCaloHits with u,v,w views

    float m_voxelWidth;  ///< Voxel box width (cm)
    float m_lengthScale; ///< The scaling factor to set all lengths to cm
    float m_energyScale; ///< The scaling factor to set all energies to GeV

    const float m_mm2cm{0.1f};          ///< mm to cm conversion
    const float m_MeV2GeV{1e-3};        ///< Geant4 MeV to GeV conversion
    const float m_voxelPathShift{1e-3}; ///< Small path shift to find next voxel
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline Parameters::Parameters() :
    m_dataFormat(Parameters::LArNDFormat::EDepSim),
    m_settingsFile(""),
    m_inputFileName(""),
    m_inputTreeName("EDepSimEvents"),
    m_geomFileName(""),
    m_geomManagerName("EDepSimGeometry"),
    m_geometryVolName("volArgonCubeDetector_PV_0"),
    m_sensitiveDetName("ArgonCube"),
    m_useModularGeometry(false),
    m_nEventsToProcess(-1),
    m_shouldDisplayEventNumber(false),
    m_shouldRunAllHitsCosmicReco(true),
    m_shouldRunStitching(true),
    m_shouldRunCosmicHitRemoval(true),
    m_shouldRunSlicing(true),
    m_shouldRunNeutrinoRecoOption(true),
    m_shouldRunCosmicRecoOption(true),
    m_shouldPerformSliceId(true),
    m_printOverallRecoStatus(false),
    m_nEventsToSkip(0),
    m_maxMergedVoxels(-1),
    m_minVoxelMipEquivE(0.3f),
    m_use3D(true),
    m_useLArTPC(false),
    m_voxelWidth(0.4f),
    m_lengthScale(m_mm2cm),
    m_energyScale(m_MeV2GeV)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Create the detector geometry based on the C++ root file
 *
 *  @param  parameters The application parameters
 *  @param  pPrimaryPandora The address of the primary pandora instance
 *  @param  geom Simple representation of the geometry for assigning TPC numbers
 */
void CreateGeometry(const Parameters &parameters, const pandora::Pandora *const pPrimaryPandora, LArNDGeomSimple &geom);

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Recursively search for volumes with the target name
 *
 *  @param  pSimGeom pointer to the input geometry
 *  @param  targetName the volume name that we want to find
 *  @param  nodeVector daughter indices to recreate the path to the target nodes
 *  @param  currentPath path to the current position in the geometry
 */
void RecursiveGeometrySearch(TGeoManager *pSimGeom, const std::string &targetName, std::vector<std::vector<unsigned int>> &nodePaths,
    std::vector<unsigned int> &currentPath);

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Create and register a tpc in pandora
 *
 *  @param  pPrimaryPandora The address of the primary pandora instance
 *  @param  parameters The application parameters
 *  @param  geom Simple representation of the geometry for assigning TPC numbers
 *  @param  pVolMatrix matrix required to convert TPC coordinates to world
 *  @param  targetNode pointer to the TPC geometry node
 *  @param  tpcNumber the number for the TPC volume
 */
void MakePandoraTPC(const pandora::Pandora *const pPrimaryPandora, const Parameters &parameters, LArNDGeomSimple &geom,
    const std::unique_ptr<TGeoHMatrix> &pVolMatrix, const TGeoNode *targetNode, const unsigned int tpcNumber);

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Process events using the supplied pandora instance
 *
 *  @param  parameters The application parameters
 *  @param  pPrimaryPandora The address of the primary pandora instance
 *  @param  geom Simple representation of the geometry for assigning TPC numbers
 */
void ProcessEvents(const Parameters &parameters, const pandora::Pandora *const pPrimaryPandora, const LArNDGeomSimple &geom);

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Process events using the supplied pandora instance, assuming EDepSim format
 *
 *  @param  parameters The application parameters
 *  @param  pPrimaryPandora The address of the primary pandora instance
 *  @param  geom Simple representation of the geometry for assigning TPC numbers
 */
void ProcessEDepSimEvents(const Parameters &parameters, const pandora::Pandora *const pPrimaryPandora, const LArNDGeomSimple &geom);

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Process events using the supplied pandora instance, assuming SED format
 *
 *  @param  parameters The application parameters
 *  @param  pPrimaryPandora The address of the primary pandora instance
 *  @param  geom Simple representation of the geometry for assigning TPC numbers
 */
void ProcessSEDEvents(const Parameters &parameters, const pandora::Pandora *const pPrimaryPandora, const LArNDGeomSimple &geom);

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Process events using the supplied pandora instance, assuming SP format
 *
 *  @param  parameters The application parameters
 *  @param  pPrimaryPandora The address of the primary pandora instance
 */
void ProcessSPEvents(const Parameters &parameters, const pandora::Pandora *const pPrimaryPandora);

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Create MC particles from the Geant4 trajectories, assuming EDepSim format
 *
 *  @param  event The Geant4 event
 *  @param  pPrimaryPandora The address of the primary pandora instance
 *  @param  parameters The application parameters
 *
 *  @return Map of <trackID, energy> for the MC particles
 */
MCParticleEnergyMap CreateEDepSimMCParticles(const TG4Event &event, const pandora::Pandora *const pPrimaryPandora, const Parameters &parameters);

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Create MC particles from the Geant4 trajectories, assuming SED format
 *
 *  @param  larsed The LArSED data object
 *  @param  pPrimaryPandora The address of the primary pandora instance
 *  @param  parameters The application parameters
 */
void CreateSEDMCParticles(const LArSED &larsed, const pandora::Pandora *const pPrimaryPandora, const Parameters &parameters);

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Convert the GENIE neutrino reaction string to a Nuance-like integer code
 *
 *  @param  reaction The neutrino reaction string
 *
 *  @return The reaction integer code
 */
int GetNuanceCode(const std::string &reaction);

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Get Nuance reaction string for SED format
 *
 *  @param  ccnc Integer specifying CCNC reaction
 *  @param  mode Integer specifying general physics mode
 *
 *  @return Name of the Nuance reaction
 */
std::string GetNuanceReaction(const int ccnc, const int mode);

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Create the LArGrid for making voxels using the geometry
 *
 *  @param  pPrimaryPandora address of the primary pandora instance
 *  @param  parameters the application parameters
 *
 *  @return the LArGrid object
 */
LArGrid MakeVoxelisationGrid(const pandora::Pandora *const pPrimaryPandora, const Parameters &parameters);

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Create the LArGrid for making voxels using the geometry
 *
 *  @param  geom simple representation of the geometry
 *  @param  parameters the application parameters
 *
 *  @return the LArGrid object
 */
LArGrid MakeVoxelisationGrid(const LArNDGeomSimple &geom, const Parameters &parameters);

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Make voxels from a given Geant4 energy deposition step
 *
 *  @param  hitInfo Information about the hit
 *  @param  grid Voxelisation grid
 *  @param  parameters The application parameters
 *  @param  simple geometry information
 *
 *  @return vector of LArVoxels
 */
LArVoxelList MakeVoxels(const LArHitInfo &hitInfo, const LArGrid &grid, const Parameters &parameters, const LArNDGeomSimple &geom);

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Combine energies for voxels with the same ID
 *
 *  @param  voxelList The unmerged list (vector) of voxels
 *
 *  @return vector of merged LArVoxels
 */
LArVoxelList MergeSameVoxels(const LArVoxelList &voxelList);

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Combine energies for voxel projections with the same (wire,drift) position
 *
 *  @param  hits The unmerged list (vector) of voxel projections
 *
 *  @return vector of merged LArVoxelProjections
 */
LArVoxelProjectionList MergeSameProjections(const LArVoxelProjectionList &hits);

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Create the pandora calohits from voxels
 *
 *  @param  voxels the voxels to use to create the hits
 *  @param  mcEnergyMap map of mc particle to its energy
 *  @param  pPrimaryPandora address of the primary pandora instance
 *  @param  parameters the application parameters
 *  @param  hitCounter reference to keep track of the number of hits
 */
void MakeCaloHitsFromVoxels(const LArVoxelList &voxels, const MCParticleEnergyMap &mcEnergyMap,
    const pandora::Pandora *const pPrimaryPandora, const Parameters &parameters, int &hitCounter);

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Return the fraction of the MC particle energy in this voxel
 *
 *  @param  mcEnergyMap map of true particle to the energy
 *  @param  voxelE energy deposited in the voxel
 *  @param  trackID true id of the track that deposited the energy
 *
 *  @return fraction of true particle energy in the voxel as a float
 */
float GetMCEnergyFraction(const MCParticleEnergyMap &mcEnergyMap, const float voxelE, const int trackID);

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  construct a generic set of calohit parameters
 *
 *  @param  voxelWidth size of the voxelisation used
 *
 *  @return LArCaloHitParameters object
 */
lar_content::LArCaloHitParameters MakeDefaultCaloHitParams(float voxelWidth);

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Parse the command line arguments, setting the application parameters
 *
 *  @param  argc argument count
 *  @param  argv argument vector
 *  @param  parameters to receive the application parameters
 *
 *  @return success
 */
bool ParseCommandLine(int argc, char *argv[], Parameters &parameters);

/**
 *  @brief  Print the list of configurable options
 *
 *  @return false, to force abort
 */
bool PrintOptions();

/**
 *  @brief  Process view option so that 3x2D and/or 3D hit lists are created
 *
 *  @param  viewOption the view option string
 *  @param  parameters to receive the application parameters
 *
 */
void ProcessViewOption(const std::string &viewOption, Parameters &parameters);

/**
 *  @brief  Process the provided reco option string to perform high-level steering
 *
 *  @param  recoOption the reco option string
 *  @param  parameters to receive the application parameters
 *
 *  @return success
 */
bool ProcessRecoOption(const std::string &recoOption, Parameters &parameters);

/**
 *  @brief  Process the optional data format and geometry input file
 *
 *  @param  formatOption the data format option string
 *  @param  inputTreeName the name of the input TTree containing the hits
 *  @param  geomFileName the name of the file containing the TGeoManager info
 *  @param  parameters to receive the application parameters
 */
void ProcessFormatOption(const std::string &formatOption, const std::string &inputTreeName, const std::string &geomFileName,
                         const std::string &geomVolName, const std::string &sensDetName, Parameters &parameters);

/**
 *  @brief  Process list of external, commandline parameters to be passed to specific algorithms
 *
 *  @param  parameters the parameters
 *  @param  pPandora the address of the pandora instance
 */
void ProcessExternalParameters(const Parameters &parameters, const pandora::Pandora *const pPandora);

} // namespace lar_nd_reco

#endif // #ifndef PANDORA_ND_INTERFACE_H
