// clang-format off
// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author:  Axel Kohlmeyer (Temple U)
   Currently maintained by:  Giacomo Fiorin (NIH)
------------------------------------------------------------------------- */

#include "fix_colvars.h"

#include "atom.h"
#include "citeme.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "memory.h"
#include "modify.h"
#include "output.h"
#include "respa.h"
#include "universe.h"
#include "update.h"

#include <cstring>
#include <iostream>
#include <memory>
#include <vector>

#include "colvarproxy_lammps.h"
#include "colvarmodule.h"
#include "colvarscript.h"


/* struct for packed data communication of coordinates and forces. */
struct LAMMPS_NS::commdata {
  int tag,type;
  double x,y,z,m,q;
};

inline std::ostream & operator<< (std::ostream &out, const LAMMPS_NS::commdata &cd)
{
  out << " (" << cd.tag << "/" << cd.type << ": "
      << cd.x << ", " << cd.y << ", " << cd.z << ") ";
  return out;
}

/* re-usable integer hash table code with static linkage. */

/** hash table top level data structure */
typedef struct inthash_t {
  struct inthash_node_t **bucket;        /* array of hash nodes */
  int size;                           /* size of the array */
  int entries;                        /* number of entries in table */
  int downshift;                      /* shift cound, used in hash function */
  int mask;                           /* used to select bits for hashing */
} inthash_t;

/** hash table node data structure */
typedef struct inthash_node_t {
  int data;                           /* data in hash node */
  int key;                            /* key for hash lookup */
  struct inthash_node_t *next;        /* next node in hash chain */
} inthash_node_t;

#define HASH_FAIL  -1
#define HASH_LIMIT  0.5

/* initialize new hash table  */
static void inthash_init(inthash_t *tptr, int buckets);
/* lookup entry in hash table */
static int inthash_lookup(void *tptr, int key);
/* insert an entry into hash table. */
static int inthash_insert(inthash_t *tptr, int key, int data);
/* delete the hash table */
static void inthash_destroy(inthash_t *tptr);

/************************************************************************
 * integer hash code:
 ************************************************************************/

/* inthash() - Hash function returns a hash number for a given key.
 * tptr: Pointer to a hash table, key: The key to create a hash number for */
static int inthash(const inthash_t *tptr, int key) {
  int hashvalue;

  hashvalue = (((key*1103515249)>>tptr->downshift) & tptr->mask);
  if (hashvalue < 0) {
    hashvalue = 0;
  }

  return hashvalue;
}

/*
 *  rebuild_table_int() - Create new hash table when old one fills up.
 *
 *  tptr: Pointer to a hash table
 */
static void rebuild_table_int(inthash_t *tptr) {
  inthash_node_t **old_bucket, *old_hash, *tmp;
  int old_size, h, i;

  old_bucket=tptr->bucket;
  old_size=tptr->size;

  /* create a new table and rehash old buckets */
  inthash_init(tptr, old_size<<1);
  for (i=0; i<old_size; i++) {
    old_hash=old_bucket[i];
    while (old_hash) {
      tmp=old_hash;
      old_hash=old_hash->next;
      h=inthash(tptr, tmp->key);
      tmp->next=tptr->bucket[h];
      tptr->bucket[h]=tmp;
      tptr->entries++;
    } /* while */
  } /* for */

  /* free memory used by old table */
  free(old_bucket);
}

/*
 *  inthash_init() - Initialize a new hash table.
 *
 *  tptr: Pointer to the hash table to initialize
 *  buckets: The number of initial buckets to create
 */
void inthash_init(inthash_t *tptr, int buckets) {

  /* make sure we allocate something */
  if (buckets==0)
    buckets=16;

  /* initialize the table */
  tptr->entries=0;
  tptr->size=2;
  tptr->mask=1;
  tptr->downshift=29;

  /* ensure buckets is a power of 2 */
  while (tptr->size<buckets) {
    tptr->size<<=1;
    tptr->mask=(tptr->mask<<1)+1;
    tptr->downshift--;
  } /* while */

  /* allocate memory for table */
  tptr->bucket=(inthash_node_t **) calloc(tptr->size, sizeof(inthash_node_t *));
}

/*
 *  inthash_lookup() - Lookup an entry in the hash table and return a pointer to
 *    it or HASH_FAIL if it wasn't found.
 *
 *  tptr: Pointer to the hash table
 *  key: The key to lookup
 */
int inthash_lookup(void *ptr, int key) {
  const inthash_t *tptr = (const inthash_t *) ptr;
  int h;
  inthash_node_t *node;


  /* find the entry in the hash table */
  h=inthash(tptr, key);
  for (node=tptr->bucket[h]; node!=nullptr; node=node->next) {
    if (node->key == key)
      break;
  }

  /* return the entry if it exists, or HASH_FAIL */
  return(node ? node->data : HASH_FAIL);
}

/*
 *  inthash_insert() - Insert an entry into the hash table.  If the entry already
 *  exists return a pointer to it, otherwise return HASH_FAIL.
 *
 *  tptr: A pointer to the hash table
 *  key: The key to insert into the hash table
 *  data: A pointer to the data to insert into the hash table
 */
int inthash_insert(inthash_t *tptr, int key, int data) {
  int tmp;
  inthash_node_t *node;
  int h;

  /* check to see if the entry exists */
  if ((tmp=inthash_lookup(tptr, key)) != HASH_FAIL)
    return(tmp);

  /* expand the table if needed */
  while (tptr->entries>=HASH_LIMIT*tptr->size)
    rebuild_table_int(tptr);

  /* insert the new entry */
  h=inthash(tptr, key);
  node=(struct inthash_node_t *) malloc(sizeof(inthash_node_t));
  node->data=data;
  node->key=key;
  node->next=tptr->bucket[h];
  tptr->bucket[h]=node;
  tptr->entries++;

  return HASH_FAIL;
}

/*
 * inthash_destroy() - Delete the entire table, and all remaining entries.
 *
 */
void inthash_destroy(inthash_t *tptr) {
  inthash_node_t *node, *last;
  int i;

  for (i=0; i<tptr->size; i++) {
    node = tptr->bucket[i];
    while (node != nullptr) {
      last = node;
      node = node->next;
      free(last);
    }
  }

  /* free the entire array of buckets */
  if (tptr->bucket != nullptr) {
    free(tptr->bucket);
    memset(tptr, 0, sizeof(inthash_t));
  }
}

/***************************************************************/

using namespace LAMMPS_NS;
using namespace FixConst;

// initialize static class members
int FixColvars::instances=0;

/***************************************************************
 create class and parse arguments in LAMMPS script. Syntax:

 fix ID group-ID colvars <config_file> [optional flags...]

 optional keyword value pairs:

  input   <input prefix>    (for restarting/continuing, defaults to
                             nullptr, but set to <output prefix> at end)
  output  <output prefix>   (defaults to 'out')
  seed    <integer>         (seed for RNG, defaults to '1966')
  tstat   <fix label>       (label of thermostatting fix)

 ***************************************************************/

FixColvars::FixColvars(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 4)
    error->all(FLERR,"Illegal fix colvars command: too few arguments");

  if (instances > 0)
    error->all(FLERR,"Only one colvars fix can be active at a time");
  ++instances;

  scalar_flag = 1;
  global_freq = 1;
  nevery = 1;
  extscalar = 1;
  restart_global = 1;
  energy_global_flag = 1;

  me = comm->me;
  root2root = MPI_COMM_NULL;
  proxy = nullptr;

  if (strcmp(arg[3], "none") == 0) {
    conf_file = nullptr;
  } else {
    conf_file = utils::strdup(arg[3]);
  }

  rng_seed = 1966;
  unwrap_flag = 1;

  inp_name = nullptr;
  out_name = nullptr;
  tfix_name = nullptr;

  /* initialize various state variables. */
  energy = 0.0;
  nlevels_respa = 0;
  init_flag = 0;
  num_coords = 0;
  comm_buf = nullptr;
  taglist = nullptr;
  force_buf = nullptr;
  idmap = nullptr;

  script_args[0] = reinterpret_cast<unsigned char *>(strdup("fix_modify"));

  parse_fix_arguments(narg, arg, true);

  if (!out_name) out_name = utils::strdup("out");

  if (me == 0) {
#ifdef LAMMPS_BIGBIG
    utils::logmesg(lmp,
                   "colvars: Warning: cannot handle atom ids > 2147483647\n");
#endif
    proxy = new colvarproxy_lammps(lmp);
    proxy->init();
    if (conf_file) {
      proxy->add_config("configfile", conf_file);
    }
  }

  /* storage required to communicate a single coordinate or force. */
  size_one = sizeof(struct commdata);
}


int FixColvars::parse_fix_arguments(int narg, char **arg, bool fix_constructor)
{
  int const iarg_start = fix_constructor ? 4 : 0;
  int iarg = iarg_start;
  while (iarg < narg) {

    bool is_fix_keyword = false;

    if (0 == strcmp(arg[iarg], "input")) {
      inp_name = utils::strdup(arg[iarg+1]);
      // input prefix is set in FixColvars::setup()
      is_fix_keyword = true;
    } else if (0 == strcmp(arg[iarg], "output")) {
      out_name = utils::strdup(arg[iarg+1]);
      // output prefix is set in FixColvars::setup()
      is_fix_keyword = true;
    } else if (0 == strcmp(arg[iarg], "seed")) {
      rng_seed = utils::inumeric(FLERR, arg[iarg+1], false, lmp);
      if (me == 0) proxy->set_random_seed(rng_seed);
      is_fix_keyword = true;
    } else if (0 == strcmp(arg[iarg], "unwrap")) {
      unwrap_flag = utils::logical(FLERR, arg[iarg+1], false, lmp);
      is_fix_keyword = true;
    } else if (0 == strcmp(arg[iarg], "tstat")) {
      tfix_name = utils::strdup(arg[iarg+1]);
      if (me == 0) set_thermostat_temperature();
      is_fix_keyword = true;
    }

    if (is_fix_keyword) {

      // Valid LAMMPS fix keyword: raise error if it has no argument
      if (iarg + 1 == narg) {
        if (fix_constructor) {
          error->all(FLERR, ("Missing argument to keyword \""+
                             std::string(arg[iarg]) +"\""));
        } else {
          // Error code consistent with Fix::modify_param()
          return 0;
        }
      }

    } else {

      if (fix_constructor) {
        error->all(FLERR, "Unrecognized fix colvars argument: please note that "
                   "Colvars script commands are not allowed until after the "
                   "fix is created");
      } else {
        if (iarg > iarg_start) {
          error->all(FLERR,
                     "Unrecognized fix colvars argument: please note that "
                     "you cannot combine fix colvars keywords and Colvars "
                     "script commands in the same line");
        } else {
          // Return negative error code to try the Colvars script commands
          return -1;
        }
      }
    }

    iarg += 2;
  }

  return iarg;
}


FixColvars::~FixColvars()
{
  delete[] conf_file;
  delete[] inp_name;
  delete[] out_name;
  delete[] tfix_name;
  memory->sfree(comm_buf);

  if (proxy) {
    delete proxy;
    inthash_t *hashtable = (inthash_t *)idmap;
    inthash_destroy(hashtable);
    delete hashtable;
  }

  if (root2root != MPI_COMM_NULL)
    MPI_Comm_free(&root2root);

  --instances;
}

/* ---------------------------------------------------------------------- */

int FixColvars::setmask()
{
  int mask = 0;
  mask |= MIN_POST_FORCE;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= END_OF_STEP;
  mask |= POST_RUN;
  return mask;
}


void FixColvars::init()
{
  if (atom->tag_enable == 0)
    error->all(FLERR, "Cannot use fix colvars without atom IDs");

  if (atom->map_style == Atom::MAP_NONE)
    error->all(FLERR, "Fix colvars requires an atom map, see atom_modify");

  if ((me == 0) && (update->whichflag == 2))
    error->warning(FLERR, "Using fix colvars with minimization");

  if (utils::strmatch(update->integrate_style, "^respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;

  if (init_flag) return;
  init_flag = 1;

  if (universe->nworlds > 1) {
    // create inter root communicator
    int color = 1;
    if (me == 0) {
      color = 0;
    }
    MPI_Comm_split(universe->uworld, color, universe->iworld, &root2root);
    if (me == 0) {
      proxy->set_replicas_communicator(root2root);
    }
  }
}


void FixColvars::set_thermostat_temperature()
{
  if (me == 0) {
    // Try to determine thermostat target temperature
    double t_target = 0.0;
    if (tfix_name) {
      if (strcmp(tfix_name, "NULL") != 0) {
        Fix *tstat_fix = modify->get_fix_by_id(tfix_name);
        if (!tstat_fix) {
          error->one(FLERR, "Could not find thermostat fix ID {}", tfix_name);
        }
        int tmp = 0;
        double *tt = reinterpret_cast<double *>(tstat_fix->extract("t_target",
                                                                   tmp));
        if (tt) {
          t_target = *tt;
          proxy->set_target_temperature(t_target);
        } else {
          error->one(FLERR, "Fix ID {} is not a thermostat fix", tfix_name);
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixColvars::init_taglist()
{
  int new_taglist_size = -1;

  if (me == 0) {

    // Number of atoms requested by Colvars
    num_coords = static_cast<int>(proxy->modify_atom_positions()->size());

    if (proxy->modified_atom_list()) {
      new_taglist_size = num_coords;
      proxy->reset_modified_atom_list();
    } else {
      new_taglist_size = -1;
    }
  }

  // Broadcast number of colvar atoms; negative means no updates
  MPI_Bcast(&new_taglist_size, 1, MPI_INT, 0, world);

  if (new_taglist_size < 0) {
    return;
  }

  num_coords = new_taglist_size;

  if (taglist) {
    memory->destroy(taglist);
    memory->destroy(force_buf);
  }
  memory->create(taglist, num_coords, "colvars:taglist");
  memory->create(force_buf, 3*num_coords, "colvars:force_buf");

  if (me == 0) {

    // Initialize and build hashtable on MPI rank 0

    std::vector<int> const &tl = *(proxy->get_atom_ids());

    if (idmap) {
      delete reinterpret_cast<inthash_t *>(idmap);
      idmap = nullptr;
    }

    inthash_t *hashtable = new inthash_t;
    inthash_init(hashtable, num_coords);
    for (int i = 0; i < num_coords; ++i) {
      taglist[i] = tl[i];
      inthash_insert(hashtable, tl[i], i);
    }

    idmap = reinterpret_cast<void *>(hashtable);
  }

  // Broadcast colvar atom ID list
  MPI_Bcast(taglist, num_coords, MPI_LMP_TAGINT, 0, world);
}


int FixColvars::modify_param(int narg, char **arg)
{
  if (narg > 100) {
    error->one(FLERR, "Too many arguments for fix_modify command");
    return 2;
  }

  // Parse arguments to fix colvars
  int return_code = parse_fix_arguments(narg, arg, false);

  if (return_code >= 0) {
    // A fix colvars argument was detected, return directly
    return return_code;
  }

  // Any unknown arguments will go through the Colvars scripting interface
  if (me == 0) {
    int error_code = COLVARSCRIPT_OK;
    colvarscript *script = proxy->script;
    script->set_cmdline_main_cmd("fix_modify " + std::string(id));
    for (int i = 0; i < narg; i++) {
      script_args[i+1] = reinterpret_cast<unsigned char *>(arg[i]);
    }

    // Run the command through Colvars
    error_code |= script->run(narg+1, script_args);

    std::string const result = proxy->get_error_msgs() + script->str_result();
    if (result.size()) {
      std::istringstream is(result);
      std::string line;
      while (std::getline(is, line)) {
        if (lmp->screen) fprintf(lmp->screen, "%s\n", line.c_str());
        if (lmp->logfile) fprintf(lmp->logfile, "%s\n", line.c_str());
      }
    }

    return (error_code == COLVARSCRIPT_OK) ? narg : 0;

  } else {

    // Return without error, don't block Fix::modify_params()
    return narg;
  }

  return 0;
}


void FixColvars::setup_io()
{
  if (me == 0) {
    proxy->set_input_prefix(std::string(inp_name ? inp_name : ""));
    if (proxy->input_prefix().size() > 0) {
      proxy->log("Will read input state from file \""+
                 proxy->input_prefix()+".colvars.state\"");
    }

    proxy->set_output_prefix(std::string(out_name ? out_name : ""));

    // Try to extract a restart prefix from a potential restart command
    LAMMPS_NS::Output *outp = lmp->output;
    if ((outp->restart_every_single > 0) &&
        (outp->restart1 != nullptr)) {

      proxy->set_default_restart_frequency(outp->restart_every_single);
      proxy->set_restart_output_prefix(std::string(outp->restart1));

    } else if ((outp->restart_every_double > 0) &&
               (outp->restart2a != nullptr)) {

      proxy->set_default_restart_frequency(outp->restart_every_double);
      proxy->set_restart_output_prefix(std::string(outp->restart2a));
    }
  }
}


void FixColvars::setup(int vflag)
{
  const tagint * const tag  = atom->tag;
  const int * const type = atom->type;
  int i,nme,tmp,ndata;
  int nlocal = atom->nlocal;

  MPI_Status status;
  MPI_Request request;

  if (me == 0) {
    setup_io();
    proxy->parse_module_config();
    proxy->setup();
  }

  init_taglist();

  // determine size of comm buffer
  nme=0;
  for (i=0; i < num_coords; ++i) {
    const tagint k = atom->map(taglist[i]);
    if ((k >= 0) && (k < nlocal))
      ++nme;
  }

  MPI_Allreduce(&nme,&nmax,1,MPI_INT,MPI_MAX,world);
  memory->create(comm_buf,nmax,"colvars:comm_buf");

  const double * const * const x = atom->x;
  const imageint * const image = atom->image;

  const double xprd = domain->xprd;
  const double yprd = domain->yprd;
  const double zprd = domain->zprd;
  const double xy = domain->xy;
  const double xz = domain->xz;
  const double yz = domain->yz;

  if (me == 0) {

    std::vector<int>     const &id = *(proxy->get_atom_ids());
    std::vector<int>           &tp = *(proxy->modify_atom_types());
    std::vector<cvm::atom_pos> &cd = *(proxy->modify_atom_positions());
    std::vector<cvm::rvector>  &of = *(proxy->modify_atom_total_forces());
    std::vector<cvm::real>     &m  = *(proxy->modify_atom_masses());
    std::vector<cvm::real>     &q  = *(proxy->modify_atom_charges());

    // store coordinate data in holding array, clear old forces


    for (i=0; i<num_coords; ++i) {
      const tagint k = atom->map(taglist[i]);
      if ((k >= 0) && (k < nlocal)) {

        of[i].x = of[i].y = of[i].z = 0.0;

        if (unwrap_flag) {
          const int ix = (image[k] & IMGMASK) - IMGMAX;
          const int iy = (image[k] >> IMGBITS & IMGMASK) - IMGMAX;
          const int iz = (image[k] >> IMG2BITS) - IMGMAX;

          cd[i].x = x[k][0] + ix * xprd + iy * xy + iz * xz;
          cd[i].y = x[k][1] + iy * yprd + iz * yz;
          cd[i].z = x[k][2] + iz * zprd;
        } else {
          cd[i].x = x[k][0];
          cd[i].y = x[k][1];
          cd[i].z = x[k][2];
        }
        if (atom->rmass_flag) {
          m[i] = atom->rmass[k];
        } else {
          m[i] = atom->mass[type[k]];
        }
        if (atom->q_flag) {
          q[i] = atom->q[k];
        }
      }
    }

    // loop over procs to receive and apply remote data

    for (i=1; i < comm->nprocs; ++i) {
      int maxbuf = nmax*size_one;
      MPI_Irecv(comm_buf, maxbuf, MPI_BYTE, i, 0, world, &request);
      MPI_Send(&tmp, 0, MPI_INT, i, 0, world);
      MPI_Wait(&request, &status);
      MPI_Get_count(&status, MPI_BYTE, &ndata);
      ndata /= size_one;

      for (int k=0; k<ndata; ++k) {

        const int j = inthash_lookup(idmap, comm_buf[k].tag);

        if (j != HASH_FAIL) {

          tp[j] = comm_buf[k].type;

          cd[j].x = comm_buf[k].x;
          cd[j].y = comm_buf[k].y;
          cd[j].z = comm_buf[k].z;

          m[j] = comm_buf[k].m;
          q[j] = comm_buf[k].q;

          of[j].x = of[j].y = of[j].z = 0.0;
        }
      }
    }
  } else { // me != 0

    // copy coordinate data into communication buffer

    nme = 0;
    for (i=0; i<num_coords; ++i) {
      const tagint k = atom->map(taglist[i]);
      if ((k >= 0) && (k < nlocal)) {

        comm_buf[nme].tag  = tag[k];
        comm_buf[nme].type = type[k];

        if (unwrap_flag) {
          const int ix = (image[k] & IMGMASK) - IMGMAX;
          const int iy = (image[k] >> IMGBITS & IMGMASK) - IMGMAX;
          const int iz = (image[k] >> IMG2BITS) - IMGMAX;

          comm_buf[nme].x = x[k][0] + ix * xprd + iy * xy + iz * xz;
          comm_buf[nme].y = x[k][1] + iy * yprd + iz * yz;
          comm_buf[nme].z = x[k][2] + iz * zprd;
        } else {
          comm_buf[nme].x = x[k][0];
          comm_buf[nme].y = x[k][1];
          comm_buf[nme].z = x[k][2];
        }

        if (atom->rmass_flag) {
          comm_buf[nme].m = atom->rmass[k];
        } else {
          comm_buf[nme].m = atom->mass[type[k]];
        }

        if (atom->q_flag) {
          comm_buf[nme].q = atom->q[k];
        }

        ++nme;
      }
    }
    /* blocking receive to wait until it is our turn to send data. */
    MPI_Recv(&tmp, 0, MPI_INT, 0, 0, world, MPI_STATUS_IGNORE);
    MPI_Rsend(comm_buf, nme*size_one, MPI_BYTE, 0, 0, world);
  }

  // run pre-run setup in colvarproxy
  if (me == 0)
    proxy->setup();

  // initialize forces
  if (utils::strmatch(update->integrate_style,"^verlet") || (update->whichflag == 2))
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
    post_force_respa(vflag,nlevels_respa-1,0);
    ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
  }
}

/* ---------------------------------------------------------------------- */
/* Main colvars handler:
 * Send coodinates and add colvar forces to atoms. */
void FixColvars::post_force(int /*vflag*/)
{
  // some housekeeping: update status of the proxy as needed.
  if (me == 0) {
    if (proxy->want_exit())
      error->one(FLERR,"Run aborted on request from colvars module.\n");
  }

  const tagint * const tag = atom->tag;
  const double * const * const x = atom->x;
  double * const * const f = atom->f;
  const imageint * const image = atom->image;

  const double xprd = domain->xprd;
  const double yprd = domain->yprd;
  const double zprd = domain->zprd;
  const double xy = domain->xy;
  const double xz = domain->xz;
  const double yz = domain->yz;
  const int nlocal = atom->nlocal;

  /* check and potentially grow local communication buffers. */
  int i,nmax_new,nme=0;
  for (i=0; i < num_coords; ++i) {
    const tagint k = atom->map(taglist[i]);
    if ((k >= 0) && (k < nlocal))
      ++nme;
  }

  MPI_Allreduce(&nme,&nmax_new,1,MPI_INT,MPI_MAX,world);
  if (nmax_new > nmax) {
    nmax = nmax_new;
    memory->grow(comm_buf,nmax,"colvars:comm_buf");
  }

  MPI_Status status;
  MPI_Request request;
  int tmp, ndata;

  if (me == 0) {

    std::vector<cvm::atom_pos> &cd = *(proxy->modify_atom_positions());

    // store coordinate data

    for (i=0; i<num_coords; ++i) {
      const tagint k = atom->map(taglist[i]);
      if ((k >= 0) && (k < nlocal)) {

        if (unwrap_flag) {
          const int ix = (image[k] & IMGMASK) - IMGMAX;
          const int iy = (image[k] >> IMGBITS & IMGMASK) - IMGMAX;
          const int iz = (image[k] >> IMG2BITS) - IMGMAX;

          cd[i].x = x[k][0] + ix * xprd + iy * xy + iz * xz;
          cd[i].y = x[k][1] + iy * yprd + iz * yz;
          cd[i].z = x[k][2] + iz * zprd;
        } else {
          cd[i].x = x[k][0];
          cd[i].y = x[k][1];
          cd[i].z = x[k][2];
        }
      }
    }

    /* loop over procs to receive remote data */
    for (i=1; i < comm->nprocs; ++i) {
      int maxbuf = nmax*size_one;
      MPI_Irecv(comm_buf, maxbuf, MPI_BYTE, i, 0, world, &request);
      MPI_Send(&tmp, 0, MPI_INT, i, 0, world);
      MPI_Wait(&request, &status);
      MPI_Get_count(&status, MPI_BYTE, &ndata);
      ndata /= size_one;

      for (int k=0; k<ndata; ++k) {
        const int j = inthash_lookup(idmap, comm_buf[k].tag);
        if (j != HASH_FAIL) {
          cd[j].x = comm_buf[k].x;
          cd[j].y = comm_buf[k].y;
          cd[j].z = comm_buf[k].z;
        }
      }
    }

  } else { // me != 0
    /* copy coordinate data into communication buffer */
    nme = 0;
    for (i=0; i<num_coords; ++i) {
      const tagint k = atom->map(taglist[i]);
      if ((k >= 0) && (k < nlocal)) {
        comm_buf[nme].tag = tag[k];

        if (unwrap_flag) {
          const int ix = (image[k] & IMGMASK) - IMGMAX;
          const int iy = (image[k] >> IMGBITS & IMGMASK) - IMGMAX;
          const int iz = (image[k] >> IMG2BITS) - IMGMAX;

          comm_buf[nme].x = x[k][0] + ix * xprd + iy * xy + iz * xz;
          comm_buf[nme].y = x[k][1] + iy * yprd + iz * yz;
          comm_buf[nme].z = x[k][2] + iz * zprd;
        } else {
          comm_buf[nme].x = x[k][0];
          comm_buf[nme].y = x[k][1];
          comm_buf[nme].z = x[k][2];
        }

        ++nme;
      }
    }
    /* blocking receive to wait until it is our turn to send data. */
    MPI_Recv(&tmp, 0, MPI_INT, 0, 0, world, MPI_STATUS_IGNORE);
    MPI_Rsend(comm_buf, nme*size_one, MPI_BYTE, 0, 0, world);
  }

  ////////////////////////////////////////////////////////////////////////
  // call our workhorse and retrieve additional information.
  if (me == 0) {
    energy = proxy->compute();
    store_forces = proxy->total_forces_enabled();
  }
  ////////////////////////////////////////////////////////////////////////

  // broadcast store_forces flag and energy data to all processors
  MPI_Bcast(&energy, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&store_forces, 1, MPI_INT, 0, world);

  // broadcast and apply biasing forces

  if (me == 0) {

    std::vector<cvm::rvector> &fo = *(proxy->modify_atom_applied_forces());

    double *fbuf = force_buf;
    for (int j=0; j < num_coords; ++j) {
      *fbuf++ = fo[j].x;
      *fbuf++ = fo[j].y;
      *fbuf++ = fo[j].z;
    }
  }
  MPI_Bcast(force_buf, 3*num_coords, MPI_DOUBLE, 0, world);

  for (int i=0; i < num_coords; ++i) {
    const tagint k = atom->map(taglist[i]);
    if ((k >= 0) && (k < nlocal)) {
      f[k][0] += force_buf[3*i+0];
      f[k][1] += force_buf[3*i+1];
      f[k][2] += force_buf[3*i+2];
    }
  }
}

/* ---------------------------------------------------------------------- */
void FixColvars::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */
void FixColvars::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  /* only process colvar forces on the outmost RESPA level. */
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */
void FixColvars::end_of_step()
{
  if (store_forces) {

    const tagint * const tag = atom->tag;
    double * const * const f = atom->f;
    const int nlocal = atom->nlocal;

    /* check and potentially grow local communication buffers. */
    int i,nmax_new,nme=0;
    for (i=0; i < num_coords; ++i) {
      const tagint k = atom->map(taglist[i]);
      if ((k >= 0) && (k < nlocal))
        ++nme;
    }

    MPI_Allreduce(&nme,&nmax_new,1,MPI_INT,MPI_MAX,world);
    if (nmax_new > nmax) {
      nmax = nmax_new;
      memory->grow(comm_buf,nmax,"colvars:comm_buf");
    }

    MPI_Status status;
    MPI_Request request;
    int tmp, ndata;

    if (me == 0) {

      // store old force data
      std::vector<cvm::rvector> &of = *(proxy->modify_atom_total_forces());

      for (i=0; i<num_coords; ++i) {
        const tagint k = atom->map(taglist[i]);
        if ((k >= 0) && (k < nlocal)) {

          const int j = inthash_lookup(idmap, tag[k]);
          if (j != HASH_FAIL) {
            of[j].x = f[k][0];
            of[j].y = f[k][1];
            of[j].z = f[k][2];
          }
        }
      }

      /* loop over procs to receive remote data */
      for (i=1; i < comm->nprocs; ++i) {
        int maxbuf = nmax*size_one;
        MPI_Irecv(comm_buf, maxbuf, MPI_BYTE, i, 0, world, &request);
        MPI_Send(&tmp, 0, MPI_INT, i, 0, world);
        MPI_Wait(&request, &status);
        MPI_Get_count(&status, MPI_BYTE, &ndata);
        ndata /= size_one;

        for (int k=0; k<ndata; ++k) {
          const int j = inthash_lookup(idmap, comm_buf[k].tag);
          if (j != HASH_FAIL) {
            of[j].x = comm_buf[k].x;
            of[j].y = comm_buf[k].y;
            of[j].z = comm_buf[k].z;
          }
        }
      }

    } else { // me != 0
      /* copy total force data into communication buffer */
      nme = 0;
      for (i=0; i<num_coords; ++i) {
        const tagint k = atom->map(taglist[i]);
        if ((k >= 0) && (k < nlocal)) {
          comm_buf[nme].tag  = tag[k];
          comm_buf[nme].x    = f[k][0];
          comm_buf[nme].y    = f[k][1];
          comm_buf[nme].z    = f[k][2];
          ++nme;
        }
      }
      /* blocking receive to wait until it is our turn to send data. */
      MPI_Recv(&tmp, 0, MPI_INT, 0, 0, world, MPI_STATUS_IGNORE);
      MPI_Rsend(comm_buf, nme*size_one, MPI_BYTE, 0, 0, world);
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixColvars::write_restart(FILE *fp)
{
  if (me == 0) {
    std::string rest_text;
    proxy->serialize_status(rest_text);
    // TODO call write_output_files()
    const char *cvm_state = rest_text.c_str();
    int len = strlen(cvm_state) + 1; // need to include terminating null byte.
    fwrite(&len,sizeof(int),1,fp);
    fwrite(cvm_state,1,len,fp);
  }
}

/* ---------------------------------------------------------------------- */

void FixColvars::restart(char *buf)
{
  if (me == 0) {
    std::string rest_text(buf);
    proxy->deserialize_status(rest_text);
  }
}

/* ---------------------------------------------------------------------- */

void FixColvars::post_run()
{
  if (me == 0) {
    proxy->post_run();
    if (lmp->citeme) {
      lmp->citeme->add(proxy->colvars->feature_report(1));
    }
  }
}

/* ---------------------------------------------------------------------- */

double FixColvars::compute_scalar()
{
  return energy;
}

/* ---------------------------------------------------------------------- */
/* local memory usage. approximately. */
double FixColvars::memory_usage()
{
  double bytes = (double) (num_coords * (2*sizeof(int)+3*sizeof(double)));
  bytes += (double)(double) (nmax*size_one) + sizeof(this);
  return bytes;
}
