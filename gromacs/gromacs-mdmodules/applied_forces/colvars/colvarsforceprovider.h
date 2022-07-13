#ifndef GMX_APPLIED_FORCES_COLVARSFORCEPROVIDER_H
#define GMX_APPLIED_FORCES_COLVARSFORCEPROVIDER_H

#include <memory>


#include "gromacs/fileio/checkpoint.h"
#include "gromacs/mdspan/extensions.h"
#include "gromacs/mdtypes/iforceprovider.h"
#include "gromacs/utility/classhelpers.h"
#include "gromacs/domdec/localatomset.h"
#include "gromacs/domdec/localatomsetmanager.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/random/tabulatednormaldistribution.h"
#include "gromacs/random/threefry.h"

#include "external/colvars/colvarmodule.h"
#include "external/colvars/colvaratoms.h"

#include "colvarproxygromacs.h"


namespace gmx
{


/*! \internal \brief
 * Implements IForceProvider for colvars.
 * Override the ColvarProxyGromacs generic class for the communication.
 */
class ColvarsForceProvider final :  public ColvarProxyGromacs, public IForceProvider
{

private:

    //! colvars data
    bool total_force_requested;
    double bias_energy;

    bool gmx_bNS; // Is this a neighbor-search step?


    // Node-local bookkepping data
    //! The colvars atom indices
    std::unique_ptr<gmx::LocalAtomSet> colvars_atoms;
    //! Total number of Colvars atoms
    int        n_colvars_atoms = 0;
    //! Unwrapped positions for all Colvars atoms, communicated to all nodes.
    rvec      *x_colvars_unwrapped = nullptr;
    //! Shifts for all Colvars atoms, to make molecule(s) whole.
    ivec      *xa_shifts = nullptr;
    //! Extra shifts since last DD step.
    ivec      *xa_eshifts = nullptr;
    //! Old positions for all Colvars atoms on master.
    rvec      *xa_old_whole = nullptr;
    //! Position of each local atom in the collective array.
    int       *xa_ind = nullptr;
    //! Bias forces on all Colvars atoms
    rvec      *f_colvars = nullptr;


public:
    friend class cvm::atom;
    //! Construct force provider for colvars from its parameters
    ColvarsForceProvider(const std::string&    colvarsConfigString,
                         LocalAtomSetManager*  localAtomSetManager,
                         PbcType               pbcType,
                         double                simulationTimeStep,
                         t_atoms               atoms,
                         const t_commrec*      cr,
                         const MDLogger*       logger,
                         const std::vector<RVec>&     colvarsCoords,
                         const std::string&    outputPrefix,
                         const std::map<std::string, std::string> & KVTInputs);

    /*! \brief Calculate colvars forces
     * \param[in] forceProviderInput input for force provider
     * \param[out] forceProviderOutput output for force provider
     */
    void calculateForces(const ForceProviderInput& forceProviderInput,
                         ForceProviderOutput*      forceProviderOutput) override;

    // Compute virial tensor for position r and force f, and add to matrix vir
    void add_virial_term(matrix vir, rvec const f, gmx::RVec const r);

    /*! \brief Write internal colvars data to checkpoint file.
     * \param[in] checkpointWriting enables writing to the Key-Value-Tree
     *                              that is used for storing the checkpoint
     *                              information
     * \param[in] moduleName names the module that is checkpointing this force-provider
     *
     * \note The provided state to checkpoint has to change if checkpointing
     *       is moved before the force provider call in the MD-loop.
     */
    //void writeCheckpointData(MdModulesWriteCheckpointData checkpointWriting, const std::string& moduleName);


    //!From colvarproxy

    /*! \brief add energy to the total count of bias energy bias_energy
    * \param[in] energy the value of energy to add
    *
    */
    void add_energy (cvm::real energy) override;

    //!From colvars_io
    std::istream &input_stream(std::string const &input_name,
                               std::string const description,
                               bool error_on_fail) override;

};

} // namespace gmx

#endif // GMX_APPLIED_FORCES_COLVARSFORCEPROVIDER_H
