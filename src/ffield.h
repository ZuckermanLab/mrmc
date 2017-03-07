
#ifndef __MAIN_H__
#define __MAIN_H__



#define _CRT_SECURE_NO_WARNINGS


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "util.h" //for subset

//#include "fragments.h"


//#define USE_GPU_AFTER_EQUIL 0// 1 for GPU solvation calculation, 0 for CPU solvation calculation

#define MAX_NUM_OF_ATOM_CLASSES     120//should be the same as the number of classes
#define MAX_NUM_OF_ATOM_TYPES       3000//amber force field requires more types
#define MAX_LIST_SIZE               10000
#define KBOLTZ                      1.987e-3

#define RAD_TO_DEG                  57.2957795//converts radians to degrees
#define DEG_TO_RAD                  1.74532925e-2//converts degrees to radians


#ifdef AMBER
#define COUL_CONST                  332.0522173 // I got this from consta.fcm in CHARMM
#define ELEC14_SCALE                0.833333333333
#elif CHARMM19
#define COUL_CONST                  332.0716//conversion from electron^2/Ang to kcal/mole
#define ELEC14_SCALE                0.4
#else
#error Must specify a force field (CHARMM19 or AMBER)
#endif



#define FALSE 0
#define TRUE 1

#define INVALID_ENERGY -1.0E20
#define DUMMY_ENERGY 1.0E20

#define INTERP_NONE    0
#define INTERP_1D      1
#define INTERP_3D      2
#define INTERP_GRAD_6D 3
#define INTERP_6D      4

/*energy term indices*/
//#define EN_INTERACTION    0
//#define EN_COVALENT_TABLE 1
#define EN_BOND           0
#define EN_ANGLE          1
#define EN_DIHEDRAL       2
#define EN_IMPROPER       3
#define EN_VDW_EXACT      4
#define EN_ELEC_EXACT     5
#define EN_GO             6
#define EN_TERMS          7 //adjust this if more terms are added!


typedef struct _atoms_ {
  int num;
  char name[6];
  char resName[6];
  int resNum;
  int chainNum;
  //int fragType;//this is fragment type from 0 to 26
  int type;//this is tinker atom type
  int classx;//this is tinker atom class
  int atomicNum;//this is element's atomic number as in the Mendeleev's table
  //The "evaluable" bonds do not include bonds involving the sidechain in the CG region.
  int numOfBondedAtoms;//this is from tinker xyz file, regular bonded1-2
  int bondedAtomList[6];
  bool evaluableBond[6]; //This is a set of flags determining whether the bond is "evaluable."
  int bondedParamType[6];
  //The "total" bonds include all bonds, regardless of where they are.
  int total_numOfBondedAtoms;//this is from tinker xyz file, regular bonded1-2
  int total_bondedAtomList[6];
  int numOfAngles;//these are regular bonded1-3 interactions
  int angleAtomList[3*20];
  int angleParamType[20];
  /*int numOfBonded12Atoms;//these are lbmc specific bonded1-2, some bonds are excluded
  int bonded12AtomList[6];
  int bonded12ParamType[6];
  int numOfBonded13Atoms;//this are lbmc specific bonded1-3, some angles are excluded
  int bonded13AtomList[3*20];
  int bonded13ParamType[20];*/
  int numOfImprops;
  int impropAtomList[4*20];
  int impropParamType[20];
  int numOfBonded14Atoms;
  int bonded14AtomList[4*30];
  int bonded14DihedParamType[30];
  double bonded14Flag[30];
  //int numOfBonded14AddAtoms;
  //int bonded14AddAtomList[4*50];
  //int bonded14AddDihedParamType[50];
  /*int numOfNonBondedAtoms;
  int nonBondedAtomList[MAX_LIST_SIZE];*/
  /*double epsilon;//This is taken care of through the fragment.
  double sigma;
  double q;*/
  //These fields added by justin.
  //int frag_within_res; //which fragment is this atom a part of? Initially within residue.
  int fragment;
  //int fragtype; //which fragment type is this atom.
  int fragatom; //which atom number within the fragment corresponds to this atom?
  double mass;
  double radius; //needed for solvation
  bool is_backbone; //is this atom part of the side chain?
  bool is_in_aa_region; //is this atom in the all-atom region?
}ATOMS;
//ATOMS atoms[MAX_NUM_OF_ATOMS];



typedef struct _atomTypeLookUp_ {//looks up different parameters based on tinker atom type
  int classx;//tinker atom class
  int atomicNum;//element's atomic number as in the Mendeleev's table
  double mass; //mass of the atom. Needed for center of mass calculation.
}ATOMTYPELOOKUP;


typedef struct _vdwParams_{//based on atom class as in tinker param file
  double sigma;
  double epsilon;
  double sigma14;   //The CHARMM 19 force field has special parameters for "1-4" van der Waals and electrostatic interacitons for some atoms.
  double epsilon14;
}VDWPARAMS;


typedef struct _bondParams_{
  int atomClass[2];
  double K;
  double r0;//equilbrium bond lenght
}BONDPARAMS;


typedef struct _angleParams_{
  int atomClass[3];
  double K;
  double theta0;//equilibrium bond angle
}ANGLEPARAMS;


typedef struct _impropParams_{
  int atomClass[4];
  double K;
  double chi0; //equilibrium angle --the CHARMM 19 ff needs this extra field
}IMPROPPARAMS;


typedef struct _dihedParams_{
  int atomClass[4];
  bool phase[4]; //need to store phases, true if phase == 180
  double V[4];
}DIHEDPARAMS;


struct atom_nb_entry {
    int iatom;
    int jatom;
    bool is14;
};

#ifdef SEDDD
struct seddd_params {
    double c, eps0, eps1, delta_eps, frac_vol_tol;
    double hydration_volume[MAX_NUM_OF_ATOM_CLASSES];
    double hydration_shell_thickness[MAX_NUM_OF_ATOM_CLASSES];
    bool read[MAX_NUM_OF_ATOM_CLASSES];
};
#endif

class forcefield {
    friend class topology;
    int numOfBondParams;
    int numOfAngleParams;
    int numOfImpropParams;
    int numOfDihedParams;

    VDWPARAMS vdwParams[MAX_NUM_OF_ATOM_CLASSES];
    BONDPARAMS bondParams[MAX_LIST_SIZE];
    ANGLEPARAMS angleParams[MAX_LIST_SIZE];
    IMPROPPARAMS impropParams[MAX_LIST_SIZE];
    DIHEDPARAMS dihedParams[MAX_LIST_SIZE];
    double vdwAFact[MAX_NUM_OF_ATOM_CLASSES][MAX_NUM_OF_ATOM_CLASSES];
    double vdwBFact[MAX_NUM_OF_ATOM_CLASSES][MAX_NUM_OF_ATOM_CLASSES];
    double vdwAFact14[MAX_NUM_OF_ATOM_CLASSES][MAX_NUM_OF_ATOM_CLASSES];
    double vdwBFact14[MAX_NUM_OF_ATOM_CLASSES][MAX_NUM_OF_ATOM_CLASSES];
    //void nonbond_energy(int rdie, int type1,  int type2, int is14, double dx, double dy, double dz, double * evdw, double * eelec);
    double MDE(double dihed, int type);
    double bond_energy(int type, int iatom, int jatom, double * coords);
    double angle_energy(int type, int a, int b, int c, double * coords);
    double dihedral_energy(int type, int a, int b, int c, int d, double * coords);
    double improper_energy(int type, int a, int b, int c, int d, double * coords);
public:
    ATOMTYPELOOKUP atomTypeLookUp[MAX_NUM_OF_ATOM_TYPES]; //fragment's constructor needs the atomic masses
    double chargeParams[MAX_NUM_OF_ATOM_TYPES];//based on atom types, fragment's constructor needs this to calculate dipole moments
    forcefield(char * fname);
    void nonbond_energy(int rdie, int type1,  int type2, int is14, double r2, double * evdw, double * eelec);
    double exact_interaction_energy(int pbc, double halfboxsize, double boxsize,double eps,  int rdie, int natom1, int * types1, double * coords1, int natom2, int * types2, double * coords2);
#ifdef SEDDD
    void moved_non_tabulated_energy(seddd_params * params, double cutoff2, int numOfAtoms, ATOMS * atoms, subset& movedatoms, subset& changedvol, bool do_bonds, int nb_atom_list_size, atom_nb_entry * nb_atom_list, double * coords, double * frac_volumes, double * energies);
    void non_tabulated_energy(seddd_params * params, double cutoff2, int numOfAtoms, ATOMS * atoms, int nb_atom_list_size, atom_nb_entry * nb_atom_list, double * coords, double * frac_volumes, double * energies);
    void subset_energy(seddd_params * params, double cutoff2, int numOfAtoms, ATOMS * atoms, subset& atomset, int nb_atom_list_size, atom_nb_entry * nb_atom_list,  double * coords, double * frac_volumes, double * internal_energies, double * intxn_energies);
#else
    void moved_non_tabulated_energy(double eps, int rdie, double cutoff2, int numOfAtoms, ATOMS * atoms, subset& movedatoms, bool do_bonds, int nb_atom_list_size, atom_nb_entry * nb_atom_list, double * coords, double * energies);
    void non_tabulated_energy(double eps, int rdie, double cutoff2, int numOfAtoms, ATOMS * atoms, int nb_atom_list_size, atom_nb_entry * nb_atom_list, double * coords, double * energies);
    void subset_energy(double eps, int rdie, double cutoff2, int numOfAtoms, ATOMS * atoms, subset& atomset, int nb_atom_list_size, atom_nb_entry * nb_atom_list,  double * coords, double * internal_energies, double * intxn_energies);
#endif
    void link_fragments(void);
    void find_parameters(int numOfAtoms, ATOMS * atoms);
    void build_coords(double * coords, int numOfAtoms, ATOMS * atoms, subset& valid_coords);
    //solvation stuff
    void nonbond_energy_gb(int type1, int type2, bool is14, double r2, double a1a2, double * evdw, double * eelec, double * egb);
    //void calculate_born_radii(gb_param_info * params, int natom, ATOMS * atoms, double * coords, double * reff, double * energies);
    //double calculate_sasa(gb_param_info * params, int natom, ATOMS * atoms, int * ineighbor, double * coords);
    //void create_virtual_site_list(int natom, ATOMS * atoms, int * nsite, virtual_site * * sites);
};

double angle(double a[], double b[]);
//justin altered this to return the cosine of the dihedral
double dihedral(double a[],double b[],double c[]);


#endif// __MAIN_H__
