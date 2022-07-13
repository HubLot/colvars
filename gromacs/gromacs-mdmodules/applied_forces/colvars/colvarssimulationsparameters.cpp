
#include "colvarssimulationsparameters.h"
#include "gromacs/utility/exceptions.h"

namespace gmx
{

void ColvarsSimulationsParameters::setLocalAtomSetManager(LocalAtomSetManager* localAtomSetManager)
{
    localAtomSetManager_ = localAtomSetManager;
}

LocalAtomSetManager* ColvarsSimulationsParameters::localAtomSetManager() const { return localAtomSetManager_; }

void ColvarsSimulationsParameters::setTopology(const gmx_mtop_t& mtop)
{
    gmx_atoms = gmx_mtop_global_atoms(mtop);
}

t_atoms ColvarsSimulationsParameters::topology() const { return gmx_atoms; }


void ColvarsSimulationsParameters::setPeriodicBoundaryConditionType(const PbcType& pbcType)
{
    pbcType_ = std::make_unique<PbcType>(pbcType);
}

PbcType ColvarsSimulationsParameters::periodicBoundaryConditionType()
{
    if (pbcType_ == nullptr)
    {
        GMX_THROW(InternalError(
                "Periodic boundary condition enum not set for colvars simulation."));
    }
    return *pbcType_;
}


void ColvarsSimulationsParameters::setSimulationTimeStep(double timeStep) { simulationTimeStep_ = timeStep; }

double ColvarsSimulationsParameters::simulationTimeStep() { return simulationTimeStep_; }


void ColvarsSimulationsParameters::setComm(const t_commrec& cr) { cr_ = &cr; }

const t_commrec* ColvarsSimulationsParameters::comm() const { return cr_; }


void ColvarsSimulationsParameters::setLogger(const MDLogger& logger) { logger_ = &logger; }


const MDLogger* ColvarsSimulationsParameters::logger() const
{
    if (logger_ == nullptr)
    {
        GMX_THROW(InternalError("Logger not set for Colvars simulation."));
    }
    return logger_;
}

} // namespace gmx
