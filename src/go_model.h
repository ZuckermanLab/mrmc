#ifndef go_model_H_INCLUDED
#define go_model_H_INCLUDED

//#include "fragments.h"
#include "topology.h"
#define INVALID_ENERGY -1.0E20
#define DUMMY_ENERGY 1.0E20

//U_native = -e_native * (x^-m - (m/n) (x^-n))/(1-m/n) where x=r_ij / r_native(ij)
//U_nonnative = e_nonnative * (x/r_ij)^-m where x=r_ij/r_hardcore
//m and n must be even (to avoid square root) with m>n
//cutoff -- distance in A
struct go_model_params {
    int m,n;
    double hardcore,cutoff,subcutoff,native_energy,nonnative_energy;
    //auxiliary -- scaled_native_en = -native_energy/(1-m/n)), ratio=m/n
    double scaled_native_en,scaled_nonnative_en,ratio, hardcore2,cutoff2,rsubcutoff;
    int hm,hn; //half of m, half of n
};

/*struct go_model_entry {
    //int ifragatom, jfragatom; //atoms IN FRAGMENT
    int ires, jres;
    bool native; //is it native?
    double native_distance2;
    //double low_distance2; //either ((sigma_i+sigma_j)/2)*(1-delta) or r_ij(1-delta)
    //double hi_distance2; //either R_cut or r_ij(1+delta)
};*/

class go_model_info {
    //go_model_entry * entries;
    //int find_entry(int ires, int jres);
    double * distances;//[ires][jres]
    //int get_index(fragmenttype * ifragtype, int ifragatom, fragmenttype * jfragtype, int jfragatom);
public:
    go_model_info();
    ~go_model_info();
    int nentries,nnative;

    char ifragname[MAX_FRAGMENT_NAME],jfragname[MAX_FRAGMENT_NAME];
    //int itype, jtype;//,inatom,jnatom;
    void create_contact_map(int nres, reslookup * resinfo, int natom, double * nativecoords, go_model_params * params, subset aaregion_res);
    //void set_parameters(go_model_params * params);
    double energy(int pbc, double halfboxsize, double boxsize, go_model_params * params, int nres, reslookup * resinfo, subset aaregion_res, int natom, double * coords); //need actual parameters
    //int count_native_contacts(int pbc, double halfboxsize, double boxsize, double cutoff, double ratio, fragmenttype * ifragtype, double * icoords, fragmenttype * jfragtype, double * jcoords);
    //void read_file(char * fname);
    //void write_file(char * fname);
};

//void get_alpha_carbon_list(int itype, fragmenttype * fragtype, int natom, fragmentinfo * fraginfo, int * nlist, int * * list);
void read_go_params(char * line, go_model_params * params);
void print_go_params(go_model_params params);
#endif // go_model_H_INCLUDED
