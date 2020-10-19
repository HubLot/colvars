#include <iostream>

#include "gmxpre.h"

#include "colvarsforceprovider.h"

#include <numeric>

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


ColvarsForceProvider::ColvarsForceProvider(const std::string&  fileinput,
                                            LocalAtomSetManager*  localAtomSetManager,
                                            PbcType               pbcType,
                                            double                simulationTimeStep,
                                            t_atoms               atoms,
                                            const t_commrec*      cr) :
    pbcType_(pbcType), gmx_atoms_(atoms), timestep(simulationTimeStep)
{
    colvars_atoms = std::make_unique<LocalAtomSet>(localAtomSetManager->add(std::vector<index>({337, 715})));


}

void ColvarsForceProvider::calculateForces(const ForceProviderInput& forceProviderInput,
                                                  ForceProviderOutput*      forceProviderOutput)
{
    std::cout << "Inside ColvarsForceProvider::calculateForces; step: " << timestep << std::endl;

    std::cout << "size localAtomset_ : " << colvars_atoms->numAtomsGlobal() << std::endl;
    std::cout << "list localAtomset_ : " << colvars_atoms->globalIndex()[0] << ", "
              <<  colvars_atoms->globalIndex()[1] << std::endl;
}


} // namespace gmx
