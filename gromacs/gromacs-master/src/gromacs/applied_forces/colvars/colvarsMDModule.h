#ifndef COLVARS_MDMODULE_H
#define COLVARS_MDMODULE_H

#include <memory>

namespace gmx
{

class IMDModule;

std::unique_ptr<IMDModule> createColvarsMDModule();

} // namespace gmx

#endif
