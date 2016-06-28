#ifndef TOPOLOGY_H_INCLUDED
#define TOPOLOGY_H_INCLUDED

#include <vector>
//#include "fragments.h"
//#include "solvation.h"
#include "ffield.h"
#include "topology.h"
#include "util.h"


#define MAX_FRAGMENTS_PER_RESIDUE 4
#define MAX_SEGMENTS              10
#define MAX_ATOMS_PER_RESIDUE     100
#define MAX_BONDS_PER_RESIDUE     MAX_ATOMS_PER_RESIDUE*4
#define UNKNOWN_COORD             9999.0
#define RT_NOT_ROTATABLE          0
#define RT_BACKBONE               1
#define RT_SIDECHAIN              2
#define MAX_FRAGMENT_NAME         32

using namespace std;
/*&struct fragment {
    int type; //Fragment type index.
    int is_side_chain; //Whether or not the fragment is part of the side chain.
    //Each of the following residues is -1 if
    int main_chain_prev; //Previous residue along main chain.
    int main_chain_next; //Next residue along main chain.
    int side_chain_prev; //Previous residue along side chain (towards main chain)
    int side_chain_next; //Next residue along side chain (away from main chain)
    int atoms[MAX_ATOMS_PER_FRAGMENT]; //Atoms in the fragment.
    int atn, atca, atc; //N, CA, and C atoms, for covalent table lookup
};
//Information on fragments within residues.  The arrays in this structure are indexed by the index of an atom within fragment.
struct resfraginfo {
    int fragtype;
    int is_side_chain; //is this fragment part of the side chain or the backbone?
    //Names of atoms within residue that are part of the fragment.
    char atomnames[MAX_ATOMS_PER_FRAGMENT][6];
    //int overlap; //Does this fragment overlap with a fragment in the previous residue?
    //int atomsprevious[MAX_ATOMS_PER_FRAGMENT]; //Atoms that overlap in previous residue.
    int offset[MAX_ATOMS_PER_FRAGMENT]; //Residue offset.  If the residue overlaps, this designates the atoms that are actually in the previous residue.
};*/

/*struct covtablelookup {
    int ifrag,jfrag;
    int itype,jtype;
    int */
//Definition of a residue.

struct mc_move {
    //atoms forming axis around which
    int iaxis, jaxis;
    subset movedatoms;
    //void perform(double * coords);
};

struct residuedef {
    char name[4];
    int natom; //Number of atoms per residue.
    char atomnames[MAX_ATOMS_PER_RESIDUE][6];
    int atomtypes[MAX_ATOMS_PER_RESIDUE];
    int nbond; //Number of bonds.
    char iname[MAX_BONDS_PER_RESIDUE][6];//Atoms pair for each bond.
    char jname[MAX_BONDS_PER_RESIDUE][6];
    int rottype[MAX_BONDS_PER_RESIDUE]; //Is this bond rotatable? If so, MC moves are generated based on rotation about this bond.
    int joffset[MAX_BONDS_PER_RESIDUE]; //Residue offset on "j" atom.
    //int nfrag; //Number of fragments intersecting this residue.
    //resfraginfo frags[MAX_FRAGMENTS_PER_RESIDUE];
    int branchatom; //Name of the branch atom.
};

struct reslookup {
    int restype; //Index of the definition.
    int whichseg; //Which segment is it a part of?
    int branchatom; //Alpha carbon, or other branching atom.
    int atomstart;
    int atomend;
    //int nscrot; //Number of backbone and sidechain rotatable bonds.
    int nbbrot;
    /*int iscrot[MAX_BONDS_PER_RESIDUE];
    int jscrot[MAX_BONDS_PER_RESIDUE];*/
    int ibbrot[MAX_BONDS_PER_RESIDUE];
    int jbbrot[MAX_BONDS_PER_RESIDUE]; //Is the bond part of the side chain?
    //int peptidebond; //Number of the actual peptide bond fragment.
    //int sidechainstart; //Starting fragment of the side chain.
    //int sidechainend; //Ending fragment of the side chain.
};

struct topology {
    //int nfragtypes;
    //fragmenttype * * fragtypes; //[MAX_FRAGMENT_TYPES];
    /*number of actual fragments, total number of atoms*/
    int natom,nseg,nres,nscrot;
    /*which type of fragments each fragment is*/
    //fragment * frags; //[MAX_FRAGMENTS];
    //virtual_site * virtual_sites;
    /*starting atom number for each fragment*/
    //int * fragstart;  //[MAX_FRAGMENTS];
    /*Atom information*/
    ATOMS * atoms;  //[MAX_ATOMS];
    int * iscrot;
    int * jscrot;
    int nresdef;
    residuedef * resdef;
    reslookup * resinfo;
    //int * sequence; //Sequence of residue index numbers;
    //Start and end of each segment, segment index for each residue.
    int * segstart;
    int * segend;
    int * first_main_chain_frag;
    //Chain code for each segment. To match PDB file.
    char * chaincodes;
    bool * closefragments;
    double qsystem;
    subset ligand;
    int ligand_res;
    //int * whichseg;
    topology(const char * commandfile, forcefield * ffield);
    ~topology();
    int frag_type_by_name(const char * name);
    int frag_type_by_file(const char * fname);
    int resdefbyname(const char * name);

    void read_residue_definition(FILE * f, residuedef * def);
    int find_atom(int actual_res, const char * aname);
    int find_atom(char chain, int res, const char * aname);
    //int is_bonded(int ifrag, int jfrag);
    //void addseg(int nsegres, char * seq);
    void print_detailed_info(subset aaregion_res);
    void print_summary_info(void);
    void insert_residue(const char * res, subset aaregion_res);
    void link_fragments(void);
    void create_angle_dihedral_lists(bool using_cov_tables);
    void create_improper_dihedral_lists(bool using_cov_tables, forcefield * ffield);
    void create_non_tab_list(bool using_cov_tables,std::vector<atom_nb_entry> * atom_nb_list);
    bool use_covalent_table(int itype, int jtype);
/*    bool term_in_covalent_tables(int iatom, int jatom);
    bool term_in_covalent_tables(int iatom, int jatom, int katom);
    bool term_in_covalent_tables(int iatom, int jatom, int katom, int latom);*/
    //void create_nb_atom_exact_list(int exact, int nb_list_per_frag, int * nb_list_count, int * nonbond_list, std::vector<atom_nb_entry> * atom_nb_list);
    void add_segment(char chain, const char * sequence, subset aaregion_res);
    void assemble_fragments(double * orig_coords, double * center, double * orient, double * new_coords);
    void fit_all_fragments(double * orig_coords, double * center, double * orient, double * new_coords, double * rmsds);
    //void load_tables(const char * fmt, const char * fragfmt, table * * tables);
    //void load_covalent_tables(const char * covtablefmt, covalent_table * * cov_tables);
    void update_coords(int ifrag, double * center, double * orient, double * coords);
    void copy_frag(int ifrag, double * center1, double * orient1, double * coords1, double * center2, double * orient2, double * coords2);
    /*double exact_interaction_energy(forcefield * ffield, int pbc, double halfboxsize, double boxsize, int gb_mode, gb_param_info * gb_params,
    int frag1, int frag2, double * coords, double * born_radii);*/
    //double covalent_table_energy(double * coords, bool * moved, covalent_table * * covalent_tables);
    //void total_internal_interaction_energy(forcefield * ffield, double eps, int rdie, double * evdw_internal, double * eelec_internal);
    //I/O routines.
    void read_pdb_stream(FILE * input, double * coords);
    void read_pdb_file(char * fname, double * coords);
    void write_pdb_stream(FILE * output, double * coords);
    void write_pdb_file(char * fname, double * coords);
    void write_pqr_file(char * fname, double * coords, int ichargedfrag, forcefield * ffield);
    void write_psf_file(char * fname, forcefield * ffield);
    //Monte Carlo move generation.
    void create_backbone_move(bool * moved, int * atom1, int * atom2);
    void create_sidechain_move(bool * moved, int * atom1, int * atom2);
    void create_backrub_move(bool * moved, int * atom1, int * atom2);
    void get_moved_atoms(bool * movedfrag, bool * movedatoms);
    void print_atom_subset(subset atomset);
    subset follow_bonds(int iatom, const subset exclude);
    void generate_backbone_moves(vector<mc_move> * backbone_moves);
    void generate_sidechain_moves(vector<mc_move> * sidechain_moves);
    void generate_backrub_moves(vector<mc_move> * backrub_moves);
    //SASA stuff.
    //void add_to_overlap_list(int ifrag, int jfrag, double * coords, std::vector<atom_nb_entry> * overlap_list);
    //void create_overlap_list(int nb_list_per_frag, int * nb_list_count, int * nonbond_list, double * coords, std::vector<atom_nb_entry> * overlap_list);
    //born radii stuff
    /*void partial_get_born_radii(bool per_atom, facts_param_info * facts_params, int n, const double * coords,  double * __restrict__ numer, double * __restrict__ denom, double * __restrict__ a, double * __restrict__ b);
    void finish_born_radii(bool per_atom, facts_param_info * facts_params, int n,  forcefield * ffield, const double tau, const double * a, const double * b, double * __restrict__ born_radii, double * __restrict__ dg_self, double * __restrict__ sasa);void convert_born_radii(forcefield * ffield, const double * per_atom_born_radii, double * per_fragment_born_radii);
    double born_solvation_energy1(forcefield * ffield, int ifrag, int jfrag, double tau, double * coords, double * per_atom_born_radii);
    double born_solvation_energy2(forcefield * ffield, int ifrag, int jfrag, double tau, double rkl2, double * coords, double xx);
    double born_solvation_energy4(forcefield * ffield, int ifrag, int jfrag, double tau, double rkl2, double * coords, double akal);*/
};

#endif // TOPOLOGY_H_INCLUDED
