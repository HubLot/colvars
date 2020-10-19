#include "gmxpre.h"

#include <memory>
#include <string>
#include <iostream>

#include "gromacs/domdec/localatomset.h"
#include "gromacs/domdec/localatomsetmanager.h"
#include "gromacs/fileio/checkpoint.h"
#include "gromacs/mdtypes/imdmodule.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/mdmodulenotification.h"
#include "gromacs/topology/mtop_util.h"

#include "colvarsMDModule.h"
#include "colvarsoptions.h"
#include "colvarssimulationsparameters.h"
#include "colvarsforceprovider.h"


namespace gmx
{

namespace
{

/*! \internal
 * \brief Colvars module
 *
 * Class that implements the colvars forces
 */
class ColvarsMDModule final : public IMDModule
{
public:
    ColvarsMDModule() {}

    //! From IMDModule
    IMdpOptionProvider* mdpOptionProvider() override { return &colvarsOptions_; }
    //! From IMDModule
    //! Colvars provide its own output
    IMDOutputProvider*  outputProvider() override { return nullptr; }

    //! From IMDModule
    //! Add this module to the force providers if active
    void initForceProviders(ForceProviders* forceProviders) override
    {
        std::cout << "ColvarsMDModule::initForceProviders()" << std::endl;

        if(colvarsOptions_.isActive())
        {

            colvarsForceProvider_ = std::make_unique<ColvarsForceProvider>(colvarsOptions_.colvarsFileName(),
                                                                           ColvarsSimulationsParameters_.localAtomSetManager(),
                                                                           ColvarsSimulationsParameters_.periodicBoundaryConditionType(),
                                                                           ColvarsSimulationsParameters_.simulationTimeStep(),
                                                                           ColvarsSimulationsParameters_.topology(),
                                                                           ColvarsSimulationsParameters_.comm());
            forceProviders->addForceProvider(colvarsForceProvider_.get());
        }
    }

    //! From IMDModule
    void subscribeToPreProcessingNotifications(MdModulesNotifier* notifier) override {}

    //! From IMDModule
    /*! \brief Request to be notified.
     * The colvars code subscribes to these notifications:
     *   - the LocalAtomSetManager sets in the simulation parameter setup
     *     by taking a LocalAtomSetManager * as parameter
     *   - the type of periodic boundary conditions that are used
     *     by taking a PeriodicBoundaryConditionType as parameter
     *   - the topology of the system
     *     by taking a gmx_mtop_t * as parameter
     *   - the communicator
     *     by taking a t_commrec as parameter
     */
    void subscribeToSimulationSetupNotifications(MdModulesNotifier* notifier) override
    {
        std::cout << "Inside subscribeToSimulationSetupNotifications()" << std::endl;
        if (!isActive())
        {
            return;
        }

        // Retrieve the LocalAtomSetManager during simulation setup
        const auto setLocalAtomManagerFunction = [this](LocalAtomSetManager* localAtomSetManager) {
            this->ColvarsSimulationsParameters_.setLocalAtomSetManager(localAtomSetManager);
        };
        notifier->simulationSetupNotifications_.subscribe(setLocalAtomManagerFunction);

        // constructing PBC during simulation setup
        const auto setPeriodicBoundaryContionsFunction = [this](const PbcType& pbc) {
            this->ColvarsSimulationsParameters_.setPeriodicBoundaryConditionType(pbc);
        };
        notifier->simulationSetupNotifications_.subscribe(setPeriodicBoundaryContionsFunction);

        // Retrieve the topology during simulation setup
        const auto setTopologyFunction = [this](gmx_mtop_t* mtop) {
            this->ColvarsSimulationsParameters_.setTopology(mtop);
        };
        notifier->simulationSetupNotifications_.subscribe(setTopologyFunction);

        // Retrieve the Communication Record during simulations setup
        const auto setCommFunction = [this](const t_commrec& cr) {
            this->ColvarsSimulationsParameters_.setComm(cr);
        };
        notifier->simulationSetupNotifications_.subscribe(setCommFunction);

        // setting the simulation time step
        const auto setSimulationTimeStepFunction = [this](const SimulationTimeStep& simulationTimeStep) {
            this->ColvarsSimulationsParameters_.setSimulationTimeStep(simulationTimeStep.delta_t);
        };
        notifier->simulationSetupNotifications_.subscribe(setSimulationTimeStepFunction);

    }


private:
    //! Return whether or not to apply a colvar biais
    bool isActive() const;

    //! The options provided for colvars
    ColvarsOptions colvarsOptions_;

    //! Parameters that become available at simulation setup time.
    ColvarsSimulationsParameters ColvarsSimulationsParameters_;

    //! Object that evaluates the forces
    std::unique_ptr<ColvarsForceProvider> colvarsForceProvider_;

};


bool ColvarsMDModule::isActive() const
{
    return colvarsOptions_.isActive();
}

} // namespace

std::unique_ptr<IMDModule> createColvarsMDModule()
{
    std::cout << "createColvarsMDModule (created in mpi process)" << std::endl;
    return std::make_unique<ColvarsMDModule>();
}

} // namespace gmx
