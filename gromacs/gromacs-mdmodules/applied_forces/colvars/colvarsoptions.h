#ifndef GMX_APPLIED_FORCES_COLVARSOPTIONS_H
#define GMX_APPLIED_FORCES_COLVARSOPTIONS_H

#include <string>
#include <vector>
#include <map>

#include "gromacs/mdtypes/imdpoptionprovider.h"
#include "gromacs/utility/real.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/utility/logger.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/utility/mdmodulesnotifiers.h"




namespace gmx
{

class KeyValueTreeObject;
class KeyValueTreeBuilder;
struct CoordinatesAndBoxPreprocessed;
struct MdRunInputFilename;


/*! \internal
 * \brief Input data storage for colvars
 */
class ColvarsOptions final : public IMdpOptionProvider
{
public:
    //! From IMdpOptionProvider
    void initMdpTransform(IKeyValueTreeTransformRules* rules) override;

    /*! \brief
     * Build mdp parameters for colvars to be output after pre-processing.
     * \param[in, out] builder the builder for the mdp options output KV-tree.
     */
    void buildMdpOutput(KeyValueTreeObjectBuilder* builder) const override;

    /*! \brief
     * Connect option name and data.
     */
    void initMdpOptions(IOptionsContainerWithSections* options) override;

    //! Store the paramers that are not mdp options in the tpr file
    void writeInternalParametersToKvt(KeyValueTreeObjectBuilder treeBuilder);

    //! Set the internal parameters that are stored in the tpr file
    void readInternalParametersFromKvt(const KeyValueTreeObject& tree);

    /*! \brief Store the topology of the system.
     * \param[in,out] mtop topology object
     */
    void getTopology(gmx_mtop_t* mtop);

    /*! \brief Process coordinates, PbcType and Box in order to validate the colvars input.
     * \param[in] coord structure with coordinates and box dimensions
     */
    void processCoordinates(const CoordinatesAndBoxPreprocessed& coord);

    //! Set the MDLogger instance
    void setLogger(const MDLogger& logger);

    /*! \brief Process MdRunInputFilename notification during mdrun.
     * In case output_prefix is empty sets it to tpr name
     * \param[in] tprFilename name of the *.tpr file that mdrun simulates
     */
    void processTprFilename(const MdRunInputFilename& tprFilename);

    //! Report if this colvars module is active
    bool isActive() const;

    //! Return the file name of the colvars input
    const std::string& colvarsFileName() const;

    //! Return the file name of the colvars restart input
    const std::string& colvarsRestartFileName() const;

    //! Return the content of the colvars input file
    const std::string& colvarsInputContent() const;

    //! Return the colvars atoms coordinates
    const std::vector<RVec>& colvarsAtomCoords() const;

    //! Return the prefix for output colvars files
    const std::string & colvarsOutputPrefix() const;


    const std::map<std::string, std::string> & colvarsInputFiles() const;


private:

    //! Indicate if colvars module is active
    bool active_ = false;

    /*! \brief Following Tags denotes names of parameters from .mdp file
     * \note Changing this strings will break .tpr backwards compability
     */
    //! \{
    const std::string c_activeTag_ = "active";
    const std::string colvarsFileNameTag_ = "filename";
    const std::string colvarsRestartFileNameTag_ = "filename-restart";
    //! \}

    //! Colvars input filename, default colvars.dat
    std::string colvarsFileName_ = "colvars.dat";
    //! Colvars input restart filename
    std::string colvarsRestartFileName_ = "";



    //! Content of the colvars input file
    std::string colvarsConfigString = "";
    //! Topology of the system
    t_atoms gmx_atoms;
    //! Coordinates
    ArrayRef<const RVec> x;
    //! PBC Type
    PbcType pbc;
    //! Box
    matrix box;
    //! Vector with colvars atoms coordinates
    std::vector<RVec> colvarsAtomCoords_;
    //!Inputs files saved as strings inside KVT
    std::map<std::string, std::string> inputFiles;


    //! Logger instance
    const MDLogger* logger_ = nullptr;

    /*! \brief String containing the prefix for output colvars files
     * default value empty, means will be deduced from *.tpr name during mdrun
     */
    std::string output_prefix_;

};

} // namespace gmx

#endif

