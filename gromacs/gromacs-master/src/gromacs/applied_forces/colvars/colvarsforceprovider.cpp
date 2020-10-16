#include <iostream>

#include "gmxpre.h"

#include "colvarsforceprovider.h"

#include <numeric>
#include <optional>

#include "gromacs/domdec/localatomset.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/math/coordinatetransformation.h"
#include "gromacs/math/densityfit.h"
#include "gromacs/math/densityfittingforce.h"
#include "gromacs/math/gausstransform.h"
#include "gromacs/mdlib/broadcaststructs.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/utility/strconvert.h"



namespace gmx
{


ColvarsForceProvider::~ColvarsForceProvider() = default;

ColvarsForceProvider::ColvarsForceProvider() {}

// ColvarsForceProvider::ColvarsForceProvider(const DensityFittingParameters& parameters,
//                                                          basic_mdspan<const float, dynamicExtents3D> referenceDensity,
//                                                          const TranslateAndScale& transformationToDensityLattice,
//                                                          const LocalAtomSet& localAtomSet,
//                                                          PbcType             pbcType,
//                                                          double              simulationTimeStep,
//                                                          const DensityFittingForceProviderState& state) :
//     impl_(new Impl(parameters, referenceDensity, transformationToDensityLattice, localAtomSet, pbcType, simulationTimeStep, state))
// {
// }

void ColvarsForceProvider::calculateForces(const ForceProviderInput& forceProviderInput,
                                                  ForceProviderOutput*      forceProviderOutput)
{
    std::cout << "Inside ColvarsForceProvider::calculateForces" << std::endl;
}


} // namespace gmx
