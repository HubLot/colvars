#include <iostream>
#include <string>

#include "gmxpre.h"

#include "colvarspreprocessor.h"

#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/stringutil.h"

#include "external/colvars/colvarmodule.h"
#include "external/colvars/colvaratoms.h"

namespace gmx
{

ColvarsPreProcessor::ColvarsPreProcessor(const std::string&     colvarsConfigString,
                                         t_atoms                atoms,
                                         PbcType               pbcType,
                                         const MDLogger*       logger,
                                         const matrix          box,
                                         ArrayRef<const RVec>   x) :
  ColvarProxyGromacs(colvarsConfigString, atoms, pbcType, logger, true,std::map<std::string, std::string>()),
  x_(x)
{

    // Initialize t_pbc struct
    set_pbc(&gmx_pbc, pbcType, box);

    cvm::log(cvm::line_marker);
    cvm::log("End colvars Initialization.\n\n");

}

std::vector<RVec> ColvarsPreProcessor::getColvarsCoords()
{

    std::vector<RVec> colvarsCoords;

    for (const auto& atom_id : atoms_ids)
    {
        colvarsCoords.push_back(x_[atom_id]);
    }
    return colvarsCoords;
}

bool ColvarsPreProcessor::inputStreamsToKVT(KeyValueTreeObjectBuilder treeBuilder, std::string tag)
{

    for ( const auto &p : cached_input_buffers_ )
    {
        treeBuilder.addValue<std::string>(tag+"-"+p.first, p.second);
    }
    return true;
}


} // namespace gmx
