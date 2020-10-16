#ifndef GMX_APPLIED_FORCES_COLVARSFORCEPROVIDER_H
#define GMX_APPLIED_FORCES_COLVARSFORCEPROVIDER_H

#include <memory>

#include "gromacs/fileio/checkpoint.h"
#include "gromacs/mdspan/extensions.h"
#include "gromacs/mdtypes/iforceprovider.h"
#include "gromacs/utility/classhelpers.h"

enum class PbcType : int;

namespace gmx
{


/*! \internal \brief
 * Implements IForceProvider for colvars forces.
 */
class ColvarsForceProvider final : public IForceProvider
{
public:
    //! Construct force provider for colvars from its parameters
    // ColvarsForceProvider(const DensityFittingParameters&             parameters,
    //                             basic_mdspan<const float, dynamicExtents3D> referenceDensity,
    //                             const TranslateAndScale& transformationToDensityLattice,
    //                             const LocalAtomSet&      localAtomSet,
    //                             PbcType                  pbcType,
    //                             double                   simulationTimeStep,
    //                             const DensityFittingForceProviderState& state);
    ColvarsForceProvider();
    ~ColvarsForceProvider();

    /*!\brief Calculate forces that maximise goodness-of-fit with a reference density map.
     * \param[in] forceProviderInput input for force provider
     * \param[out] forceProviderOutput output for force provider
     */
    void calculateForces(const ForceProviderInput& forceProviderInput,
                         ForceProviderOutput*      forceProviderOutput) override;

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

private:

};

} // namespace gmx

#endif // GMX_APPLIED_FORCES_DENSITYFITTINGFORCEPROVIDER_H
