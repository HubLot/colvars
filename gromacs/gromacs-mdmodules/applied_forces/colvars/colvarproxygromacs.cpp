#include <iostream>

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
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/mdlib/groupcoord.h"
#include "gromacs/utility/exceptions.h"

#include "colvarproxygromacs.h"


namespace gmx
{

ColvarProxyGromacs::ColvarProxyGromacs(const std::string&  colvarsConfigString,
                                       t_atoms             atoms,
                                       PbcType             pbcType,
                                       const MDLogger*     logger,
                                       bool                doParsing,
                                       const std::map<std::string, std::string> &input_strings) :
   gmx_atoms(atoms), pbcType_(pbcType), logger_(logger), doParsing_(doParsing)
{

    restart_frequency_s = 0;

    //!From colvarproxy

    // Retrieve masses and charges from input file
    updated_masses_ = updated_charges_ = true;

    // User-scripted forces are not available in GROMACS
    have_scripts = false;

    angstrom_value = 0.1;

    // Get the thermostat temperature.
    // NEEDED for parsing. would be good to be set here
    //thermostat_temperature = ir->opts.ref_t[0];

    // GROMACS random number generation.
    rng.seed(makeRandomSeed());

    //set the initial step
    // would be good to be set here

     // GROMACS timestep
    //timestep = ir->delta_t;


    // Read configuration file and set up the proxy during Pre processing
    // and during simulation phase but only on the master node.
    if (doParsing)
    {

        cache_input_buffers = true;
        cached_input_buffers_ = input_strings;

        colvars = new colvarmodule(this);
        cvm::log(cvm::line_marker);
        cvm::log("Start colvars Initialization.");

        //version_int = get_version_from_string(COLVARPROXY_VERSION);


        colvars->cite_feature("GROMACS engine");
        colvars->cite_feature("Colvars-GROMACS interface");

        if (cvm::debug()) {
            cvm::log("Initializing the colvars proxy object.\n");
        }

        cvm::log("Using GROMACS interface, version "+
        cvm::to_str(COLVARPROXY_VERSION)+".\n");

        colvars->read_config_string(colvarsConfigString);

        colvars->setup();
        colvars->setup_input();


        // Citation Reporter
        cvm::log(std::string("\n")+colvars->feature_report(0)+std::string("\n"));

        //TODO: Retrieve step
        // if (step != 0) {
        //     cvm::log("Initializing step number to "+cvm::to_str(step)+".\n");
        // }

        // colvars->it = colvars->it_restart = step;
        colvars->it = 0;
    }
}


// GROMACS uses nanometers and kJ/mol internally
cvm::real ColvarProxyGromacs::backend_angstrom_value() { return 0.1; }

// From Gnu units
// $ units -ts 'k' 'kJ/mol/K/avogadro'
// 0.0083144621
cvm::real ColvarProxyGromacs::boltzmann() { return 0.0083144621; }

// Temperature of the simulation (K)
cvm::real ColvarProxyGromacs::temperature()
{
  return thermostat_temperature;
}

// Time step of the simulation (fs)
// GROMACS uses picoseconds.
cvm::real ColvarProxyGromacs::dt() { return 1000.0*timestep; }

cvm::real ColvarProxyGromacs::rand_gaussian()
{
  return  normal_distribution(rng);
}

void ColvarProxyGromacs::log (std::string const &message)
{
    if (logger_)
    {
        GMX_LOG(logger_->info).appendText(message);
    }
}

void ColvarProxyGromacs::error (std::string const &message)
{
  // In GROMACS, all errors are fatal.
  fatal_error (message);
}

void ColvarProxyGromacs::fatal_error (std::string const &message)
{
  log(message);
  GMX_THROW(InternalError("Error in collective variables module.\n"));
}


int ColvarProxyGromacs::set_unit_system(std::string const &units_in, bool /*colvars_defined*/)
{
  if (units_in != "gromacs") {
    cvm::error("Specified unit system \"" + units_in + "\" is unsupported in Gromacs. Supported units are \"gromacs\" (nm, kJ/mol).\n");
    return COLVARS_ERROR;
  }
  return COLVARS_OK;
}


// **************** ATOMS ****************

int ColvarProxyGromacs::check_atom_id(int atom_number)
{
  // GROMACS uses zero-based arrays.
  int const aid = (atom_number-1);

  if (cvm::debug())
    log("Adding atom "+cvm::to_str(atom_number)+
        " for collective variables calculation.\n");

  if ( (aid < 0) || (aid >= gmx_atoms.nr) ) {
    cvm::error("Error: invalid atom number specified, "+
               cvm::to_str(atom_number)+"\n", COLVARS_INPUT_ERROR);
    return COLVARS_INPUT_ERROR;
  }

  return aid;
}


int ColvarProxyGromacs::init_atom(int atom_number)
{
  // GROMACS uses zero-based arrays.
  int aid = atom_number-1;

  for (size_t i = 0; i < atoms_ids.size(); i++) {
    if (atoms_ids[i] == aid) {
      // this atom id was already recorded
      atoms_ncopies[i] += 1;
      return i;
    }
  }

  aid = check_atom_id(atom_number);

  if(aid < 0) {
    return COLVARS_INPUT_ERROR;
  }

  int const index = add_atom_slot(aid);
  update_atom_properties(index);
  return index;
}

void ColvarProxyGromacs::update_atom_properties(int index)
{

  // update mass
  double const mass = gmx_atoms.atom[atoms_ids[index]].m;
  if (mass <= 0.001) {
    this->log("Warning: near-zero mass for atom "+
              cvm::to_str(atoms_ids[index]+1)+
              "; expect unstable dynamics if you apply forces to it.\n");
  }
  atoms_masses[index] = mass;
  // update charge
  atoms_charges[index] = gmx_atoms.atom[atoms_ids[index]].q;
}

ColvarProxyGromacs::~ColvarProxyGromacs()
{
  if (colvars != NULL) {
    delete colvars;
    colvars = NULL;
  }
}


cvm::rvector ColvarProxyGromacs::position_distance (cvm::atom_pos const &pos1,
                                                    cvm::atom_pos const &pos2) const
{
  rvec r1, r2, dr;
  r1[0] = pos1.x;
  r1[1] = pos1.y;
  r1[2] = pos1.z;
  r2[0] = pos2.x;
  r2[1] = pos2.y;
  r2[2] = pos2.z;

  pbc_dx(&gmx_pbc, r2, r1, dr);
  return cvm::atom_pos( dr[0], dr[1], dr[2] );
}


} // namespace gmx
