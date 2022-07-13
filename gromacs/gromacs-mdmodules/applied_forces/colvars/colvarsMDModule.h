#ifndef COLVARS_MDMODULE_H
#define COLVARS_MDMODULE_H

#include <memory>

namespace gmx
{

class IMDModule;

/*! \internal
    \brief Information about the colvars module.
 *
 * Provides name and method to create a colvars module.
 */
struct ColvarsModuleInfo
{
    /*! \brief
     * Creates a module for applying forces according to a colvar bias.
     *
     */
    static std::unique_ptr<IMDModule> create();
    //! The name of the module
    static const std::string name_;
};


} // namespace gmx

#endif
