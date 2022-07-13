#ifndef GMX_APPLIED_FORCES_COLVARSPREPROCESSOR_H
#define GMX_APPLIED_FORCES_COLVARSPREPROCESSOR_H


#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/keyvaluetreebuilder.h"

#include "colvarproxygromacs.h"



namespace gmx
{

/*! \internal \brief
 * Class that read a colvars configuration file during pre-processing and
 * retrieve the colvars atoms coordinates to be stored in tpr KVT.
 */
class ColvarsPreProcessor: public ColvarProxyGromacs
{
public:
    /*! \brief Construct ColvarsPreProcessor from its parameters
     *

     * \param[in] fileinput Content of the colvars input file.
     * \param[in] atoms Atoms topology
     * \param[in] pbcType Periodic boundary conditions
     * \param[in] logger GROMACS logger instance
     * \param[in] box Matrix with full box of the system
     * \param[in] x Coordinates of each atom in the system
     */
    ColvarsPreProcessor(const std::string&    colvarsConfigString,
                        t_atoms               atoms,
                        PbcType               pbcType,
                        const MDLogger*       logger,
                        const matrix          box,
                        ArrayRef<const RVec>  x);


    //! Return a vector of the colvars atoms coordinates
    std::vector<RVec> getColvarsCoords();

    //! Save all input files of colvars into the KVT
    bool inputStreamsToKVT(KeyValueTreeObjectBuilder treeBuilder, std::string tag);

private:

    //! Atoms coordinates of the whole system
    ArrayRef<const RVec> x_;

};


} // namespace gmx

#endif // GMX_APPLIED_FORCES_COLVARSPREPROCESSOR_H
