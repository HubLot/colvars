#ifndef GMX_APPLIED_FORCES_COLVARPROXYGROMACS_H
#define GMX_APPLIED_FORCES_COLVARPROXYGROMACS_H

#include <memory>


#include "gromacs/fileio/checkpoint.h"
#include "gromacs/mdspan/extensions.h"
#include "gromacs/utility/classhelpers.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/random/tabulatednormaldistribution.h"
#include "gromacs/random/threefry.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/utility/logger.h"

#include "external/colvars/colvarmodule.h"
#include "external/colvars/colvaratoms.h"
#include "external/colvars/colvarproxy.h"


#define COLVARPROXY_VERSION "MDMODULE-ALPHA"


/**
 * TODO:
 * - créer méthode colvars pour lire un string (contenu d'un fichier input)
 * - voir comment récup les données manquantes nécessaires à colvars (temp, step, etc)
 * - gerer le checkpointing
 **/

namespace gmx
{


/*! \internal \brief
 * Implements a GROMACS version of colvarproxy.
 * This class hold for the communication between colvars and GROMACS.
 * 2 child class will inherit from this one: one during pre processing (ColvarsPreProcessor)
 * and one during the simulation (ColvarsForceProvider).
 * Most of the work needed for the communication will be implemented in this class.
 */
class ColvarProxyGromacs : public colvarproxy
{

protected:

    //! Atoms topology
    t_atoms gmx_atoms;

    //! Box infos
    PbcType pbcType_;
    t_pbc gmx_pbc;

    //! Activate or not the parsing of the Colvars config file
    bool doParsing_;

    //! The simulation time step
    cvm::real timestep;
    //! The temperature of the system
    cvm::real thermostat_temperature;

    size_t restart_frequency_s;


    // GROMACS random number generation.
    DefaultRandomEngine           rng;   // gromacs random number generator
    TabulatedNormalDistribution<> normal_distribution;

    //GROMACS logger instance
    const MDLogger* logger_ = nullptr;

public:
    friend class cvm::atom;

    /*! \brief Construct ColvarProxyGromacs from its parameters
     *
     * \param[in] fileinput Content of the colvars input file.
     * \param[in] atoms Atoms topology
     * \param[in] pbcType Periodic boundary conditions
     * \param[in] logger GROMACS logger instance
     * \param[in] doParsing Wether the input file should be parsed.
     */
    ColvarProxyGromacs(const std::string&    colvarsConfigString,
                       t_atoms               atoms,
                       PbcType               pbcType,
                       const MDLogger*       logger,
                       bool doParsing,
                       const std::map<std::string, std::string> &input_strings);
    ~ColvarProxyGromacs();

    /// Update colvars topology of one atom mass and charge from the GROMACS topology
    void update_atom_properties(int index);

    //!From colvarproxy

    // **************** SYSTEM-WIDE PHYSICAL QUANTITIES ****************
    cvm::real backend_angstrom_value() override;
    cvm::real boltzmann() override;
    cvm::real temperature() override;
    cvm::real dt() override;
    cvm::real rand_gaussian() override;

    // **************** INPUT/OUTPUT ****************
    /// Print a message to the main log
    void log (std::string const &message) override;
    /// Print a message to the main log and let the rest of the program handle the error
    void error (std::string const &message);
    /// Print a message to the main log and exit with error code
    ///possible suppression
    void fatal_error (std::string const &message);
    /// Request to set the units used internally by Colvars
    int set_unit_system(std::string const &units_in, bool colvars_defined) override;

    /// Initialize colvars atom from GROMACS topology
    int init_atom(int atom_number) override;

    /*! \brief Check if atom belongs to the global index of atoms
     *  \param[in] atom_number Colvars index of the atom to check
     */
    int check_atom_id(int atom_number) override;

    // **************** PERIODIC BOUNDARY CONDITIONS ****************
    cvm::rvector position_distance (cvm::atom_pos const &pos1,
                                    cvm::atom_pos const &pos2) const override;

};

} // namespace gmx

#endif // GMX_APPLIED_FORCES_COLVARPROXYGROMACS_H
