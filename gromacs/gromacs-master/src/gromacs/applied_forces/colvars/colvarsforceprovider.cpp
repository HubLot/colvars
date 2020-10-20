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

#include "colvarsforceprovider.h"


namespace gmx
{

ColvarsForceProvider::ColvarsForceProvider(const std::string&  fileinput,
                                            LocalAtomSetManager*  localAtomSetManager,
                                            PbcType               pbcType,
                                            double                simulationTimeStep,
                                            t_atoms               atoms,
                                            const t_commrec*      cr) :
   pbcType(pbcType), gmx_atoms(atoms), timestep(simulationTimeStep)
{
    std::cout << "Inside ColvarsForceProvider::ColvarsForceProvider "<< std::endl;

    // Initialize colvars.
    total_force_requested = false;
    restart_frequency_s = 0;
    // Retrieve masses and charges from input file
    updated_masses_ = updated_charges_ = true;

    // User-scripted forces are not available in GROMACS
    force_script_defined = false;
    have_scripts = false;

    angstrom_value = 0.1;

    // Get the thermostat temperature.
    // TODO
    //thermostat_temperature = ir->opts.ref_t[0];

    // GROMACS random number generation.
    // TODO
    //rng.seed(ir->ld_seed);

    // Handle input filenames and prefix/suffix for colvars files.
    // TODO
    output_prefix_str = "output";
    restart_output_prefix_str = output_prefix_str + ".restart";

    //Get the initial step
    //TODO
    colvars->it = 0;

    // Read configuration file and set up the proxy only on the master node.
    if (MASTER(cr))
    {
        // initiate module: this object will be the communication proxy
        // colvarmodule pointer is only defined on the Master due to the static pointer to colvarproxy.
        colvars = new colvarmodule(this);

        //version_int = get_version_from_string(COLVARPROXY_VERSION);

        if (cvm::debug()) {
            log("Initializing the colvars proxy object.\n");
        }

        cvm::log("Using GROMACS interface, version "+
        cvm::to_str(COLVARPROXY_VERSION)+".\n");

        colvars->read_config_file(fileinput.c_str());

        colvars->setup();
        colvars->setup_input();
        colvars->setup_output();

        //TODO: REtrieve step
        // if (step != 0) {
        //     cvm::log("Initializing step number to "+cvm::to_str(step)+".\n");
        // }

        // colvars->it = colvars->it_restart = step;

    } // end master

  // MPI initialisation

    // Initialise attributs for the MPI communication
    if(MASTER(cr)) {
        // Retrieve the number of colvar atoms
        n_colvars_atoms = atoms_ids.size();
    }

    if(PAR(cr)) {
        // Let the other nodes know the number of colvar atoms and their ids to construct a gmx::LocalAtomSet
        block_bc(cr->mpi_comm_mygroup, n_colvars_atoms);
        atoms_ids.reserve(n_colvars_atoms);
        nblock_bc(cr->mpi_comm_mygroup, n_colvars_atoms, atoms_ids.data());

        // Initialise atoms_new_colvar_forces on non-master nodes
        if(!MASTER(cr)) {
        atoms_new_colvar_forces.reserve(n_colvars_atoms);
        }
    }

    colvars_atoms = std::make_unique<LocalAtomSet>(localAtomSetManager->add(atoms_ids));

    snew(x_colvars_unwrapped, n_colvars_atoms);
    snew(xa_shifts,           n_colvars_atoms);
    snew(xa_eshifts,          n_colvars_atoms);
    snew(f_colvars,           n_colvars_atoms);



    // // Communicate initial coordinates to all processes
    // if (PAR(cr))
    // {
    //     nblock_bc(cr->mpi_comm_mygroup,  n_colvars_atoms, xa_old);
    // }


    if (MASTER(cr) && cvm::debug()) {
        cvm::log ("atoms_ids = "+cvm::to_str (atoms_ids)+"\n");
        cvm::log ("atoms_ncopies = "+cvm::to_str (atoms_ncopies)+"\n");
        cvm::log ("positions = "+cvm::to_str (atoms_positions)+"\n");
        cvm::log ("total_forces = "+cvm::to_str (atoms_total_forces)+"\n");
        cvm::log ("atoms_new_colvar_forces = "+cvm::to_str (atoms_new_colvar_forces)+"\n");
        cvm::log (cvm::line_marker);
        log("done initializing the colvars proxy object.\n");
    }


}

void ColvarsForceProvider::calculateForces(const ForceProviderInput& forceProviderInput,
                                                  ForceProviderOutput*      forceProviderOutput)
{

    //Construct t_pbc struct
    set_pbc(&gmx_pbc, pbcType, forceProviderInput.box_);

    const t_commrec *cr           = &(forceProviderInput.cr_);
    // Local atom coords
    const gmx::ArrayRef<const gmx::RVec> x  = forceProviderInput.x_;
    // Local atom coords (coerced into into old gmx type)
    const rvec *x_pointer          = &(x.data()->as_vec());
    const auto& box = forceProviderInput.box_;

    //TODO: Get the real value from ForceProviderInput
    gmx_bNS = true;


    //!!!!!!!!!! DO NOT WORK ON MPI VERSION !!!!!!!!!!!!!!!!!!
    //Copy for the first time only (hence compare to nullptr) x coordinates of colvars atoms into xa_old
    //to determine shifts in communicate_group_positions
    //cant work on mpi version because we will need to reduce the coordinates first.
    //should be done in the class constructor but no coordinates avaliable.
    if(xa_old == nullptr) {
      snew(xa_old,              n_colvars_atoms);
      for (size_t i = 0; i < colvars_atoms->numAtomsGlobal(); i++)
      {
        copy_rvec(x[colvars_atoms->globalIndex()[i]], xa_old[i]);
      }
    }



    // Eventually there needs to be an interface to update local data upon neighbor search
    // We could check if by chance all atoms are in one node, and skip communication
    communicate_group_positions(cr, x_colvars_unwrapped, xa_shifts, xa_eshifts, gmx_bNS, x_pointer,
                                colvars_atoms->numAtomsGlobal(), colvars_atoms->numAtomsLocal(),
                                colvars_atoms->localIndex().data(), colvars_atoms->collectiveIndex().data(),
                                xa_old, box);

    // Communicate_group_positions takes care of removing shifts (unwrapping)
    // in single node jobs, communicate_group_positions() is efficient and adds no overhead

    if (MASTER(cr))
    {
        // On non-master nodes, jump directly to applying the forces

        // backup applied forces if necessary to calculate total forces (if available in future version of Gromacs)
        //if (total_force_requested)
        //  previous_atoms_new_colvar_forces = atoms_new_colvar_forces;

        // Zero the forces on the atoms, so that they can be accumulated by the colvars.
        for (size_t i = 0; i < atoms_new_colvar_forces.size(); i++) {
        atoms_new_colvar_forces[i].x = atoms_new_colvar_forces[i].y = atoms_new_colvar_forces[i].z = 0.0;
        }

        // Get the atom positions from the Gromacs array.
        for (size_t i = 0; i < atoms_ids.size(); i++) {
        atoms_positions[i] = cvm::rvector(x_colvars_unwrapped[i][0], x_colvars_unwrapped[i][1], x_colvars_unwrapped[i][2]);
        }

        // // Get total forces if required (if available in future version of Gromacs)
        // if (total_force_requested && cvm::step_relative() > 0) {
        //   for (size_t i = 0; i < atoms_ids.size(); i++) {
        //     size_t aid = atoms_ids[i];
        //     atoms_total_forces[i] = cvm::rvector(f[aid][0], f[aid][1], f[aid][2]);
        //   }
        // }

        bias_energy = 0.0;
        // Call the collective variable module to fill atoms_new_colvar_forces
        if (colvars->calc() != COLVARS_OK) {
        cvm::fatal_error("Error calling colvars->calc()\n");
        }

        // Copy the forces to C array for broadcasting
        for (int i = 0; i < n_colvars_atoms; i++)
        {
        f_colvars[i][0] = atoms_new_colvar_forces[i].x;
        f_colvars[i][1] = atoms_new_colvar_forces[i].y;
        f_colvars[i][2] = atoms_new_colvar_forces[i].z;
        }

        forceProviderOutput->enerd_.term[F_COM_PULL] += bias_energy;
    } // master node

    //Broadcast the forces to all the nodes
    if (PAR(cr))
    {
        nblock_bc(cr->mpi_comm_mygroup, n_colvars_atoms, f_colvars);
    }

    const gmx::ArrayRef<gmx::RVec> &f_out = forceProviderOutput->forceWithVirial_.force_;
    matrix local_colvars_virial = { { 0 } };
    const auto& localcolvarsIndex = colvars_atoms->localIndex();
    const auto& collectivecolvarsIndex = colvars_atoms->collectiveIndex();
    // Loop through local atoms to aply the colvars forces
    for (gmx::index l = 0; l < localcolvarsIndex.ssize(); l++)
    {
        /* Get the right index of the local colvars atoms */
        int i_local = localcolvarsIndex[l];
        /* Index of this local atom in the collective colvars atom arrays */
        int i_colvars = collectivecolvarsIndex[l];
        /* Add */
        rvec_inc(f_out[i_local], f_colvars[i_colvars]);
        add_virial_term(local_colvars_virial, f_colvars[i_colvars], x_colvars_unwrapped[i_colvars]);
    }

    forceProviderOutput->forceWithVirial_.addVirialContribution(local_colvars_virial);

    if(MASTER(cr))
    {
        colvars->write_output_files();
    }

    colvars->it++;

}

void ColvarsForceProvider::add_virial_term(matrix vir, rvec const f, gmx::RVec const x)
{
  for (int j = 0; j < DIM; j++) {
    for (int m = 0; m < DIM; m++) {
      vir[j][m] -= 0.5 * f[j] * x[m];
    }
  }
}

// Pass restraint energy value for current timestep to MD engine
void ColvarsForceProvider::add_energy (cvm::real energy)
{
  bias_energy += energy;
}

// GROMACS uses nanometers and kJ/mol internally
cvm::real ColvarsForceProvider::backend_angstrom_value() { return 0.1; }

// From Gnu units
// $ units -ts 'k' 'kJ/mol/K/avogadro'
// 0.0083144621
cvm::real ColvarsForceProvider::boltzmann() { return 0.0083144621; }

// Temperature of the simulation (K)
cvm::real ColvarsForceProvider::temperature()
{
  return thermostat_temperature;
}

// Time step of the simulation (fs)
// GROMACS uses picoseconds.
cvm::real ColvarsForceProvider::dt() { return 1000.0*timestep; }

cvm::real ColvarsForceProvider::rand_gaussian()
{
  return  normal_distribution(rng);
}

void ColvarsForceProvider::log (std::string const &message)
{
  // Gromacs prints messages on the stderr FILE.
  fprintf(stderr, "colvars: %s", message.c_str());
}

void ColvarsForceProvider::error (std::string const &message)
{
  // In GROMACS, all errors are fatal.
  fatal_error (message);
}

void ColvarsForceProvider::fatal_error (std::string const &message)
{
  log(message);
  if (!cvm::debug())
    log("If this error message is unclear, "
	"try recompiling with -DCOLVARS_DEBUG.\n");
  gmx_fatal(FARGS,"Error in collective variables module.\n");
}


// **************** ATOMS ****************

int ColvarsForceProvider::check_atom_id(int atom_number)
{
  // GROMACS uses zero-based arrays.
  int const aid = (atom_number-1);

  if (cvm::debug())
    log("Adding atom "+cvm::to_str(atom_number)+
        " for collective variables calculation.\n");

  if ( (aid < 0) || (aid >= gmx_atoms.nr) ) {
    cvm::error("Error: invalid atom number specified, "+
               cvm::to_str(atom_number)+"\n", INPUT_ERROR);
    return INPUT_ERROR;
  }

  return aid;
}


int ColvarsForceProvider::init_atom(int atom_number)
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
    return INPUT_ERROR;
  }

  int const index = add_atom_slot(aid);
  update_atom_properties(index);
  return index;
}

void ColvarsForceProvider::update_atom_properties(int index)
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

ColvarsForceProvider::~ColvarsForceProvider()
{
  if (colvars != NULL) {
    delete colvars;
    colvars = NULL;
  }
}

// **************** PERIODIC BOUNDARY CONDITIONS ****************
//  Get the PBC-aware distance vector between two positions
cvm::rvector ColvarsForceProvider::position_distance(cvm::atom_pos const &pos1,
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
