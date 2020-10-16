#include <iostream>

#include "gmxpre.h"

#include "colvarsoptions.h"

#include "gromacs/options/basicoptions.h"
#include "gromacs/options/optionsection.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/keyvaluetreetransform.h"
#include "gromacs/utility/strconvert.h"


namespace gmx
{

namespace
{

/*! \brief Helper to declare mdp transform rules.
 *
 * Enforces uniform mdp options that are always prepended with the correct
 * string for the densityfitting mdp options.
 *
 * \tparam ToType type to be transformed to
 * \tparam TransformWithFunctionType type of transformation function to be used
 *
 * \param[in] rules KVT transformation rules
 * \param[in] transformationFunction the function to transform the flat kvt tree
 * \param[in] optionTag string tag that describes the mdp option, appended to the
 *                      default string for the density guided simulation
 */
template<class ToType, class TransformWithFunctionType>
void colvarsMdpTransformFromString(IKeyValueTreeTransformRules* rules,
                                          TransformWithFunctionType    transformationFunction,
                                          const std::string&           optionTag)
{
    rules->addRule()
            .from<std::string>("/colvars-" + optionTag)
            .to<ToType>("/colvars/" + optionTag)
            .transformWith(transformationFunction);
}

} //namespace


void ColvarsOptions::initMdpTransform(IKeyValueTreeTransformRules* rules)
{
    const auto& stringIdentityTransform = [](std::string s) { return s; };
    colvarsMdpTransformFromString<bool>(rules, &fromStdString<bool>, c_activeTag_);
    colvarsMdpTransformFromString<std::string>(rules, stringIdentityTransform,
                                                      colvarsFileNameTag_);
    colvarsMdpTransformFromString<std::string>(rules, stringIdentityTransform,
                                                      colvarsRestartFileNameTag_);
}


void ColvarsOptions::buildMdpOutput(KeyValueTreeObjectBuilder* builder) const
{
    std::cout << "ColvarsOptions::buildMdpOutput()" << std::endl;
    builder->addValue<std::string>("comment-colvars-module", "; Colvars biais");
    builder->addValue<bool>("colvars-"+c_activeTag_, active_);
    builder->addValue<std::string>("comment-colvars-"+colvarsFileNameTag_,
                                            "; colvars input file");
    builder->addValue<std::string>("colvars-"+colvarsFileNameTag_, colvarsFileName_);
    if(!colvarsRestartFileName_.empty())
    {
        builder->addValue<std::string>("comment-colvars-"+colvarsRestartFileNameTag_,
                                                "; colvars restart file");
        builder->addValue<std::string>("colvars-"+colvarsRestartFileNameTag_, colvarsRestartFileName_);
    }
}


void ColvarsOptions::initMdpOptions(IOptionsContainerWithSections* options)
{
    std::cout << "ColvarsOptions::initMdpOptions()" << std::endl;
    auto section = options->addSection(OptionSection("colvars"));
    section.addOption(BooleanOption(c_activeTag_.c_str()).store(&active_));
    section.addOption(StringOption(colvarsFileNameTag_.c_str()).store(&colvarsFileName_));
    section.addOption(StringOption(colvarsRestartFileNameTag_.c_str()).store(&colvarsRestartFileName_));
}


bool ColvarsOptions::isActive() const
{
    return active_;
}

const std::string& ColvarsOptions::colvarsFileName() const
{
    return colvarsFileName_;
}

const std::string& ColvarsOptions::colvarsRestartFileName() const
{
    return colvarsRestartFileName_;
}

} //namespace gmx
