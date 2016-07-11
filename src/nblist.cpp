#include <math.h>
#include <stdio.h>
#include "util.h"
#include "mc.h"
#include "ffield.h"
#include "rotations.h"

//"by groups" pair list generator using residues as groups.
//1. Compute the COM and radius of each residue.
//2. Calculate the distance between COM for each pair of residues, and determine which residues can be close enough for their atoms to be on the list.
//(The minimum distance between atoms in residues i and j = dist(i,j) - radius(i) - radius(j).  If this is > listcutoff, can skip atom pairs.
//3. Create short pair list based on this.
void topology::create_pair_list(bool pbc, double halfboxsize, double boxsize, double listcutoff, std::vector<atom_nb_entry> * pair_list, double * coords)
{
    int iatom, jatom, ires, jres, k, ientry;
    double d;
    double * res_com;
    double * res_mass;
    double * res_radius;
    //bool * close_residues;
    std::vector<atom_nb_entry> * piece_of_list;
#ifdef TIMERS
    switch_timer(TIMER_NB_LIST);
#endif
    res_com = (double *) checkalloc(3*nres, sizeof(double));
    res_mass = (double *) checkalloc(nres, sizeof(double));
    res_radius = (double *) checkalloc(nres, sizeof(double));
    //1. compute COM and radius for each residue
    for (ires=0; ires<nres; ires++) {
        res_mass[ires]=0.0;
        for (k=0; k<3; k++) res_com[3*ires+k]=0.0;
        res_radius[ires]=0.0;
    }
    for (iatom=0; iatom<natom; iatom++) {
        ires=atoms[iatom].resNum;
        res_mass[ires]+=atoms[iatom].mass;
        for (k=0; k<3; k++) res_com[3*ires+k]+=atoms[iatom].mass*coords[3*iatom+k];
    }
    for (ires=0; ires<nres; ires++) for (k=0; k<3; k++) {
        res_com[3*ires+k]/=res_mass[ires];
    }
    for (iatom=0; iatom<natom; iatom++) {
        ires=atoms[iatom].resNum;
        d=pbc_distance2(pbc,halfboxsize,boxsize,&coords[3*iatom],&res_com[3*ires]);
        if (d>res_radius[ires]) res_radius[ires]=d;
    }
    for (ires=0; ires<nres; ires++) res_radius[ires]=sqrt(res_radius[ires]);
    //2. determine close residues and assemble the short list (the original list is kept by residue to make this more efficient
    pair_list->clear();
    for (ires=0; ires<nres; ires++) for (jres=ires; jres<nres; jres++) {
        if (jres==ires) {
            d=0;
        } else { //jres>ires
            d=sqrt(pbc_distance2(pbc,halfboxsize,boxsize,&res_com[3*ires],&res_com[3*jres]));
            d=d-res_radius[ires]-res_radius[jres];
        }
        if (d<listcutoff) {
            piece_of_list=&orig_pair_list_by_res[ires*nres+jres];
            pair_list->insert(pair_list->end(),piece_of_list->begin(),piece_of_list->end());
        }
    }
    //printf("Non bond pair list original entries = %ld, short list = %ld\n",non_tab_list.size(),pair_list->size());
    //clean up
    free(res_com);
    free(res_mass);
    free(res_radius);
#ifdef TIMERS
    switch_timer(TIMER_OTHER);
#endif
}
//(bool pbc, double halfboxsize, double boxsize, double cutoff, int nfrag, bool * moved, double * center)

bool simulation::update_pair_list_if_needed(long int istep, double * coords)
{
    bool redo_nb_list;
    int iatom,k;
    double rij[3],r2,margin,margin2;
    margin=(listcutoff-cutoff)/2.0;
    margin2=margin*margin;
    redo_nb_list=FALSE;
    if (last_nb_list_coords==NULL) {
        last_nb_list_coords = (double *) checkalloc(3*top->natom,sizeof(double));
        redo_nb_list=true;
    } else for (iatom=0; iatom<top->natom; iatom++) {
        r2=pbc_distance2(pbc,halfboxsize,boxsize,&coords[3*iatom],&last_nb_list_coords[3*iatom]);
        if (r2>margin2) {
            redo_nb_list=true;
            break;
        }
    }
    if (redo_nb_list) {
        //printf("Updating pair list at step %ld\n",istep);
        top->create_pair_list(pbc,halfboxsize,boxsize,listcutoff,&pair_list,coords);
    }
    return redo_nb_list;
}
