#ifndef GMX_APPLIED_FORCES_COLVARSOPTIONS_H
#define GMX_APPLIED_FORCES_COLVARSOPTIONS_H

#include <string>

#include "gromacs/mdtypes/imdpoptionprovider.h"


namespace gmx
{

class KeyValueTreeObject;
class KeyValueTreeBuilder;


/*! \internal
 * \brief Input data storage for colvars
 */
class ColvarsOptions final : public IMdpOptionProvider
{
public:
    //! From IMdpOptionProvider
    void initMdpTransform(IKeyValueTreeTransformRules* rules) override;

    /*! \brief
     * Build mdp parameters for colvars to be output after pre-processing.
     * \param[in, out] builder the builder for the mdp options output KV-tree.
     */
    void buildMdpOutput(KeyValueTreeObjectBuilder* builder) const override;

    /*! \brief
     * Connect option name and data.
     */
    void initMdpOptions(IOptionsContainerWithSections* options) override;

    //! Store the paramers that are not mdp options in the tpr file
    void writeInternalParametersToKvt(KeyValueTreeObjectBuilder treeBuilder);

    //! Set the internal parameters that are stored in the tpr file
    void readInternalParametersFromKvt(const KeyValueTreeObject& tree);

    //! Report if this set of options is active
    bool isActive() const;

    //! Return the file name of the colvars input
    const std::string& colvarsFileName() const;

    //! Return the file name of the colvars restart input
    const std::string& colvarsRestartFileName() const;

private:
    const std::string c_activeTag_ = "active";
    //! Indicate if colvars is active
    bool active_ = false;

    const std::string colvarsFileNameTag_ = "filename";
    //! Colvars input filename
    std::string colvarsFileName_ = "colvars.dat";

    const std::string colvarsRestartFileNameTag_ = "filename-restart";
    //! Colvars input restart filename
    std::string colvarsRestartFileName_ = "";
};

} // namespace gmx

#endif

