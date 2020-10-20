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

#include "gromacs/colvars/colvarmodule.h"
#include "gromacs/colvars/colvaratoms.h"
#include "gromacs/colvars/colvarproxy.h"

enum class PbcType : int;

#define COLVARPROXY_VERSION "MDMODULE-ALPHA"

namespace gmx
{


/*! \internal \brief
 * Implements IForceProvider for colvars forces.
 */
class ColvarsForceProvider final : public IForceProvider, public colvarproxy
{

private:

    PbcType pbcType;
    t_pbc gmx_pbc;

    t_atoms gmx_atoms;

    //! From colvarproxy
    //! The simulation time step
    cvm::real timestep;
    //! The temperature of the system
    //TODO: get value
    cvm::real thermostat_temperature;

    //!Other variables
    bool total_force_requested;
    double bias_energy;
    size_t restart_frequency_s;

    bool gmx_bNS; // Is this a neighbor-search step? Eventually will become unnecessary

    // GROMACS random number generation.
    DefaultRandomEngine           rng;   // gromacs random number generator
    TabulatedNormalDistribution<> normal_distribution;

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
    //! Last whole positions for all Colvars atoms.
    rvec      *xa_old = nullptr;
    //! Bias forces on all Colvars atoms
    rvec      *f_colvars = nullptr;


public:
    friend class cvm::atom;
    //! Construct force provider for colvars from its parameters
    ColvarsForceProvider(const std::string&    fileinput,
                         LocalAtomSetManager*  localAtomSetManager,
                         PbcType               pbcType,
                         double                simulationTimeStep,
                         t_atoms               atoms,
                         const t_commrec*      cr);
    ~ColvarsForceProvider();

    /*!\brief Calculate forces that maximise goodness-of-fit with a reference density map.
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



    void add_energy (cvm::real energy) override;

    // **************** SYSTEM-WIDE PHYSICAL QUANTITIES ****************
    cvm::real backend_angstrom_value() override;
    cvm::real boltzmann() override;
    cvm::real temperature() override;
    cvm::real dt() override;
    cvm::real rand_gaussian() override;

    // **************** INPUT/OUTPUT ****************
    /// Print a message to the main log
    void log (std::string const &message) override;
    /// Print a message to the main log and let the rest of the program handle the error
    void error (std::string const &message) override;
    /// Print a message to the main log and exit with error code
    void fatal_error (std::string const &message);
    /// Request to set the units used internally by Colvars
    int set_unit_system(std::string const &units_in, bool colvars_defined) override {}

    int init_atom(int atom_number) override;

    int check_atom_id(int atom_number) override;
    void update_atom_properties(int index);
    cvm::rvector position_distance(cvm::atom_pos const &pos1,
                                   cvm::atom_pos const &pos2) const;


};

} // namespace gmx

#endif // GMX_APPLIED_FORCES_COLVARSFORCEPROVIDER_H
