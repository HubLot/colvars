
#ifndef GMX_APPLIED_FORCES_COLVARSIMULATIONSPARAMETERS_H
#define GMX_APPLIED_FORCES_COLVARSIMULATIONSPARAMETERS_H

#include "gromacs/domdec/localatomset.h"
#include "gromacs/domdec/localatomsetmanager.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/pbcutil/pbc.h"

namespace gmx
{

/*! \internal
 * \brief Collect colvars parameters only available during simulation setup.
 *
 * To build the density fitting force provider during simulation setup,
 * one needs access to parameters that become available only during simulation setup.
 *
 * This class collects these parameters via MdModuleNotifications in the
 * simulation setup phase and provides a check if all necessary parameters have
 * been provided.
 */
class ColvarsSimulationsParameters
{
public:
    ColvarsSimulationsParameters() = default;

    //! Set the local atom set Manager for colvars.
    void setLocalAtomSetManager(LocalAtomSetManager* localAtomSetManager);
    //! Get the local atom set Manager for colvars.
    LocalAtomSetManager* localAtomSetManager() const;


    /*! \brief Construct the topology of the system.
     *
     * \param[in] mtop is the pointer to the global topology struct
     */
    void setTopology(gmx_mtop_t* mtop);

    //! Get the topology
    t_atoms topology() const;

    /*! \brief Set the periodic boundary condition via MdModuleNotifier.
     *
     * The pbc type is wrapped in PeriodicBoundaryConditionType to
     * allow the MdModuleNotifier to statically distinguish the callback
     * function type from other 'int' function callbacks.
     *
     * \param[in] pbcType enumerates the periodic boundary condition.
     */
    void setPeriodicBoundaryConditionType(const PbcType& pbcType);

    //! Get the periodic boundary conditions
    PbcType periodicBoundaryConditionType();

    //! Set the simulation time step
    void setSimulationTimeStep(double timeStep);
    //! Return the simulation time step
    double simulationTimeStep();

    //! Set the communicator
    void setComm(const t_commrec& cr);
    //! Return the communicator
    const t_commrec* comm() const;

private:

    //! The LocalAtomSetManager
    LocalAtomSetManager* localAtomSetManager_;
    //! The type of periodic boundary conditions in the simulation
    std::unique_ptr<PbcType> pbcType_;
    //! The simulation time step
    double simulationTimeStep_ = 1;
    //! The topology
    t_atoms gmx_atoms;
    //! The communicator
    const t_commrec* cr_;

    GMX_DISALLOW_COPY_AND_ASSIGN(ColvarsSimulationsParameters);
};

} // namespace gmx

#endif