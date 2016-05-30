#ifndef MC_H_INCLUDED
#define MC_H_INCLUDED

#include <stdio.h>
#include <vector>
#if defined(PARALLEL) || defined(EXCHANGE)
#include <mpi.h>
#endif
//#include "tables.h"
#include "ffield.h"
#include "topology.h"
#include "go_model.h"
#include "util.h"
//#include "fragments.h"
//#include "covalent_tables.h"
//#include "nblist.h"

#define MOVE_BACKBONE      1
#define MOVE_SIDECHAIN     2
#define MOVE_BACKRUB       3
#define MOVE_LIGAND_TRANS  4
#define MOVE_LIGAND_ROT    5
#define MOVE_HEAVY_TRANS   6
#define MOVE_HEAVY_ROT     7
#define NUM_MOVES          7


//For some reason the above include files are not prividing these declarations.
//class table;
//class fragmenttype;
using std::vector;

static const char * mc_move_names[NUM_MOVES+1] = {"","Backbone","Sidechain","Backrub","Ligand-trans","Ligand-rot","Heavy-atom-trans","Heavy-atom-rot"};


/*File names: coordinate output, quaternion output, starting restart file (if needed), ending restart file*/
class simulation {
private:
    //table * * tables; //[MAX_FRAGMENT_TYPES][MAX_FRAGMENT_TYPES];
    //covalent_table * * covalent_tables;
    //use_std_tables, use_cov_tables,
    bool pbc, interp, enwrite, rdie;//,do_mc,do_energy,do_dock;
    double eps, cutoff, cutoff2, boxsize, halfboxsize, rmargin, tables_lambda;
    char forcefieldfname[255],deffname[255];
    subset aaregion_res, aaregion_atoms;
    bool aaregion_specified, initialized;
    forcefield * ffield;
    go_model_info * go_model;
    go_model_params go_params;
    //long int entablecount, enexactcount, enevalcount;
    FILE * energy_output;
    FILE * pairs_output;
    char structfname[255],struct2fname[255],xyzfname[255],quatfname[255],restartfname[255];
    FILE * xyzoutput;
    FILE * quatoutput;
    double * initcoords;
    /*double * newcenter;
    double * neworient;*/
    double * newcoords;
    /*double * oldcenter;
    double * oldorient;*/
    double * oldcoords;
    /*other data*/
    double beta; /*Monte Carlo temperature, expressed as beta*/

    long int nmcstep, nsave, nprint,ncheck,nprevstep; /*number of MC steps, frequency of saving conformations, frequency of printing*/
    double prob[NUM_MOVES+1],cumprob[NUM_MOVES+1],movesize[NUM_MOVES+1];
    vector<mc_move> backbone_moves;
    vector<mc_move> sidechain_moves;
    vector<mc_move> backrub_moves;
    int startoption;
    topology * top;
    char * sequence;
    char ligand_resname[6];
    /*bool use_nb_list;
    fragment_nblist * frag_nblist;*/
    float * xwrite; //the DCD file format requires single precision
    float * ywrite;
    float * zwrite;
    unsigned long int seed;

    std::vector<atom_nb_entry> non_tab_list, overlap_list;
    //double * en_by_table;
    int nres; //nres is for temporary use until the topology object can be fully initialized
#if defined(PARALLEL) || defined(EXCHANGE)
    int mynod, numnod;
    void parallel_start(void);
#endif
#ifdef EXCHANGE
    int myrep, nrep, exchfreq;
    double * betas;
    MPI_File rexlog;
    void exchange_init(FILE * input, char * fragfmt);
    void exchange(int icycle, double * cum_energies, double * cum_energy, double * center, double * orient, double * coords);
#endif
    //std::vector<atom_nb_entry> nb_atom_list;
    //double evdw_internal, eelec_internal; //Total internal VDW and electrostatic energies over all the fragments.
    void recenter(void);
    void create_lists(void);
    //void write_frame_xyz(FILE * output, long int istep, double * coords);
    //void write_frame_pdb(FILE * output, long int istep, double * coords);
    void write_pair_pdb(FILE * output, int ifrag, int jfrag, double * coords);
    void read_frame_quat(FILE * input, long int * istep, double * center, double * orient);
    void write_frame_quat(FILE * output, long int istep, double * center, double * orient);
    void copy_frag(int ifrag, double * center1, double * orient1, double * coords1, double * center2, double * orient2, double * coords2);
    double interaction_energy(int ifrag, int jfrag, double * center, double * orient,double * coords);
    void moved_energy(int movetype,subset& movedatoms, double * coords, double * energies, double * etot);
    void total_energy(double * coords, double * energies, double * etot);
    void read_restart(char * fname);
    void write_restart(long int istep, char * fname);
    void create_nonbond_list(double * center);
    bool check_nb_list(bool * moved, double * center);
    //MC move support -- these are in mcmoves.cpp
    void mcmove(int * movetype, subset * movedatoms, double * coords);
    void rotate_atoms_by_axis(mc_move * move, const double angle, double * coords);
    void rotate_atoms_by_point(subset atoms, const double * quat, const double * point, double * coords);
    void do_ligand_trans(double movesize, double * coords);
    void do_ligand_rot(double movesize, double * coords);
    void heavy_atom_trans(subset * movedatoms, double movesize, double * coords);
    void heavy_atom_rot(subset * movedatoms, double movesize, double * coords);
    void prepare_docking(double trans_size, double rot_size, double bond_rot_size, int nsearch, double * coords);
    //i/o related stuff -- these are in io.cpp
    void write_dcd_header(FILE * dcdfile);
    void write_dcd_frame(FILE * dcdfile, double * coords);
    void print_energies(FILE * output,int hdr, const char * title, long int istep, double * energies, double etot);
    void print_energies_by_table(void);
    //void load_tables(char * fmt);
public:
    void comparison_test(void);
#if defined(PARALLEL) || defined(EXCHANGE)
    simulation(int _mynod, int _numnod);
#else
    simulation();
#endif
    ~simulation();
    void process_commands(char * infname);
    void prepare_docking(double trans_size, double rot_size);
    void finish_initialization(void);
    void mcloop(void);
    void fakeloop(void);
    void do_energies(char * type, char * fname,  char * enfname);
};



#endif // MC_H_INCLUDED
