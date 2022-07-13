#include <iostream>
#include <string>

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

ColvarsForceProvider::ColvarsForceProvider(const std::string&  colvarsConfigString,
                                            LocalAtomSetManager*  localAtomSetManager,
                                            PbcType               pbcType,
                                            double                simulationTimeStep,
                                            t_atoms               atoms,
                                            const t_commrec*      cr,
                                            const MDLogger*       logger,
                                            const std::vector<RVec>&    colvarsCoords,
                                            const std::string&    outputPrefix,
                                            const std::map<std::string, std::string> &KVTInputs) :
   ColvarProxyGromacs(colvarsConfigString, atoms, pbcType, logger, MASTER(cr),KVTInputs)
{


    timestep = simulationTimeStep;

    output_prefix_str = outputPrefix;
    restart_output_prefix_str = output_prefix_str + ".restart";

    //TODO Handle restart file

    if (doParsing_)
    {
        colvars->setup_output();
    }


    // MPI initialisation

    // Initialise attributs for the MPI communication
    if(MASTER(cr)) {
        // Retrieve the number of colvar atoms
        n_colvars_atoms = atoms_ids.size();
    }

    if(PAR(cr)) {
        // Let the other nodes know the number of colvar atoms and their ids to construct a gmx::LocalAtomSet
        block_bc(cr->mpi_comm_mygroup, n_colvars_atoms);
        atoms_ids.resize(n_colvars_atoms);
        nblock_bc(cr->mpi_comm_mygroup, n_colvars_atoms, atoms_ids.data());

        // Initialise atoms_new_colvar_forces on non-master nodes
        if(!MASTER(cr)) {
            atoms_new_colvar_forces.resize(n_colvars_atoms);
        }
    }

    colvars_atoms = std::make_unique<LocalAtomSet>(localAtomSetManager->add(atoms_ids));


    snew(x_colvars_unwrapped, n_colvars_atoms);
    snew(xa_shifts,           n_colvars_atoms);
    snew(xa_eshifts,          n_colvars_atoms);
    snew(f_colvars,           n_colvars_atoms);
    snew(xa_old_whole,        n_colvars_atoms);


    for (int i = 0; i < n_colvars_atoms; i++)
    {
        copy_rvec(colvarsCoords[i], xa_old_whole[i]);
    }

    // // Communicate initial coordinates to all processes
    if (PAR(cr))
    {
        nblock_bc(cr->mpi_comm_mygroup,  n_colvars_atoms, xa_old_whole);
    }


    if (MASTER(cr) && cvm::debug()) {
        cvm::log ("atoms_ids = "+cvm::to_str (atoms_ids)+"\n");
        cvm::log ("atoms_ncopies = "+cvm::to_str (atoms_ncopies)+"\n");
        cvm::log ("positions = "+cvm::to_str (atoms_positions)+"\n");
        cvm::log ("total_forces = "+cvm::to_str (atoms_total_forces)+"\n");
        cvm::log ("atoms_new_colvar_forces = "+cvm::to_str (atoms_new_colvar_forces)+"\n");
        cvm::log (cvm::line_marker);
        log("done initializing the colvars proxy object.\n");
    }

    if(MASTER(cr))
    {
        cvm::log(cvm::line_marker);
        cvm::log("End colvars Initialization.\n\n");
    }

}

void ColvarsForceProvider::calculateForces(const ForceProviderInput& forceProviderInput,
                                                  ForceProviderOutput*      forceProviderOutput)
{

    //std::cout << "Inside ColvarsForceProvider::calculateForces step  "<< forceProviderInput.step_ << std::endl;

    //Construct t_pbc struct
    set_pbc(&gmx_pbc, pbcType_, forceProviderInput.box_);

    const t_commrec *cr           = &(forceProviderInput.cr_);
    // Local atom coords
    const gmx::ArrayRef<const gmx::RVec> x  = forceProviderInput.x_;
    // Local atom coords (coerced into into old gmx type)
    const rvec *x_pointer          = &(x.data()->as_vec());
    const auto& box = forceProviderInput.box_;

    colvars->it = forceProviderInput.step_;

    //TODO: Get the real value from ForceProviderInput
    gmx_bNS = true;

    // Eventually there needs to be an interface to update local data upon neighbor search
    // We could check if by chance all atoms are in one node, and skip communication
    communicate_group_positions(cr, x_colvars_unwrapped, xa_shifts, xa_eshifts, gmx_bNS, x_pointer,
                                colvars_atoms->numAtomsGlobal(), colvars_atoms->numAtomsLocal(),
                                colvars_atoms->localIndex().data(), colvars_atoms->collectiveIndex().data(),
                                xa_old_whole, box);


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
        cvm::error("Error calling colvars->calc()\n");
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

std::istream &ColvarsForceProvider::input_stream(std::string const &input_name,
                                                 std::string const description,
                                                 bool error_on_fail)
{

  if (!io_available()) {
    cvm::error("Error: trying to access an input file/channel "
               "from the wrong thread.\n", COLVARS_BUG_ERROR);
    return *input_stream_error_;
  }

  //Populate the istream map from the cached_input_buffers map
  if (cached_input_buffers_.count(input_name) > 0) {
      int ret = set_input_stream_from_string(input_name, cached_input_buffers_[input_name]);
      if (ret > 0 && error_on_fail) {
        cvm::error("Error: cannot open "+description+" \""+input_name+"\".\n",
               COLVARS_FILE_ERROR);
        return *input_stream_error_;
      }
  }
  else
  {
    cvm::error("Error: trying to read a cached input buffer, but caching is "
               "not enabled.\n", COLVARS_BUG_ERROR);
    return *input_stream_error_;
  }

    return *(input_streams_[input_name]);;
}


} // namespace gmx
