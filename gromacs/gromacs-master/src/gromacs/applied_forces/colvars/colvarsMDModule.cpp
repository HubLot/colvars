#include "gmxpre.h"

#include "colvarsMDModule.h"

#include <cmath>

#include <memory>
#include <string>
#include <iostream>

#include "gromacs/domdec/localatomset.h"
#include "gromacs/domdec/localatomsetmanager.h"
#include "gromacs/fileio/checkpoint.h"
#include "gromacs/mdtypes/imdmodule.h"
#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/mdmodulenotification.h"


#include "colvarsoptions.h"
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
    //! Colvars tools provide its own output
    IMDOutputProvider*  outputProvider() override { return nullptr; }
     //! Add this module to the force providers if active
    void                initForceProviders(ForceProviders* forceProviders) override
    {
        std::cout << "ColvarsMDModule::initForceProviders()" << std::endl;
        if(colvarsOptions_.isActive())
        {

            //Read colvars input here
            colvarsForceProvider_ = std::make_unique<ColvarsForceProvider>();
            forceProviders->addForceProvider(colvarsForceProvider_.get());
        }
    }

    void subscribeToPreProcessingNotifications(MdModulesNotifier* notifier) override
    {

        if (!isActive())
        {
            return;
        }
        std::cout << "Inside subscribeToPreProcessingNotifications()" << std::endl;
    }

    void subscribeToSimulationSetupNotifications(MdModulesNotifier* notifier) override
    {
        std::cout << "Inside subscribeToSimulationSetupNotifications()" << std::endl;
        if (!isActive())
        {
            return;
        }
    }

private:
    //! Return whether or not to apply a colvar biais
    bool isActive() const;

    //! The options provided for colvars
    ColvarsOptions colvarsOptions_;

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
    std::cout << "createColvarsMDModule (created for mpi process)" << std::endl;
    return std::make_unique<ColvarsMDModule>();
}

} // namespace gmx
