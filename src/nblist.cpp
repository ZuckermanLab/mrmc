#include <math.h>
#include <stdio.h>
#include "util.h"
#include "mc.h"
#include "ffield.h"
#include "rotations.h"

//This suborutine is called during initialization to set up pair lists by residue pair.  Each pair list includes all pairs of atoms meeting certain conditions belonging to that residue pair.
//The nonbond list generator called during simulations works by combining these pair lists to create the complete pair list, depending on the distance between COMs of residues.
void topology::create_non_tab_list(void)
{
    int ifrag,jfrag,j,ii,jj,k,l,m,iatom,jatom,katom,latom,temp,ires,jres;
    bool is12, is13, is14;
    atom_nb_entry newentry;
    //non_tab_list.clear();
    orig_pair_list_by_res = (std::vector<atom_nb_entry> * ) checkalloc(nres*nres,sizeof(std::vector<atom_nb_entry>));
    orig_solv_list_by_res = (std::vector<atom_nb_entry> * ) checkalloc(nres*nres,sizeof(std::vector<atom_nb_entry>));
    for (ires=0; ires<nres; ires++) for (jres=ires; jres<nres; jres++) {
	orig_pair_list_by_res[ires*nres+jres].clear();
	orig_solv_list_by_res[ires*nres+jres].clear();
    }
    for (iatom=0; iatom<natom; iatom++)
        for (jatom=iatom+1; jatom<natom; jatom++) {
                    //include in list only if both atoms don't belong to the CG region
                    //if (!atoms[iatom].is_in_aa_region && !atoms[jatom].is_in_aa_region) continue;
                    //iatom<jatom and check to see if 1-2, 1-3, or 1-4.  We do not use the angle or dihedral lists as they do not include all possibilities in the CG region.
                    is12=false;
                    is13=false;
                    is14=false;
                    for (k=0; k<atoms[iatom].numOfBondedAtoms; k++) {
                        katom=atoms[iatom].bondedAtomList[k];
                        if (katom==jatom) {
                            is12=true; //this pair is a 1-2 pair
                        } else {
                            //iatom - katom -...
                            for (l=0; l<atoms[katom].numOfBondedAtoms; l++) {
                                latom=atoms[katom].bondedAtomList[l];
                                if (latom==jatom) {
                                    is13=true; //this pair is a 1-3 pair (iatom - katom - jatom)
                                } else {
                                    for (m=0; m<atoms[latom].numOfBondedAtoms; m++)
                                        if (atoms[latom].bondedAtomList[m]==jatom) is14=true; //a 1-4 pair (iatom - katom - latom - jatom), need to flag it.
                                }
                            }
                        }
                    }
                    if (is12 || is13) continue; //skip it
                    newentry.iatom=iatom;
                    newentry.jatom=jatom;
                    //newentry.is12_or_13=is12 || is13;
                    newentry.is14=is14;
                    //non_tab_list.push_back(newentry);
                    ires=atoms[iatom].resNum;
                    jres=atoms[jatom].resNum;
                    if (ires>jres) {
                       temp=ires;
                       ires=jres;
                       jres=temp;
                    }
                    orig_pair_list_by_res[ires*nres+jres].push_back(newentry);
                    //special "solvation" list contains only heavy atom pairs -- 1-4 pairs are excluded (see eef1.src in CHARMM)
                    if ((atoms[iatom].atomicNum>1) && (atoms[jatom].atomicNum>1) && !newentry.is14) {
			orig_solv_list_by_res[ires*nres+jres].push_back(newentry);
		    }
        }
    //printf("Number of atom pairs to be evaluated exactly: %ld\n",non_tab_list.size());
}


//"by groups" pair list generator using residues as groups.
//1. Compute the COM and radius of each residue.
//2. Calculate the distance between COM for each pair of residues, and determine which residues can be close enough for their atoms to be on the list.
//(The minimum distance between atoms in residues i and j = dist(i,j) - radius(i) - radius(j).  If this is > listcutoff, can skip atom pairs.
//3. Create short pair list based on this.
//There are two pair lists generated: the regular pair list and the solvation list.
//The regular pair list contains all pairs of atoms (including hydrogens) that (a) are not in a 1-2 or 1-3 relationship, (b) at least one atom belongs to the all-atom region, (c) potentially within cutoffs.
//The solvation pair list contains all pairs of atoms that (a) are both non-hydrogen, (b) are not in a 1-2, 1-3, or 1-4 relationship, (c) potentially within cutoffs. 
//The differences are: 
//(a) the solvation pair list contains pairs in which both atoms belong to the CG region, the regular pair list does not
//(b) the regular pair list contains pairs in 1-4 relationship (flagged by the "is14" flag), the solvation pair list does not
//(c) the solvation pair list contains only pairs of heavy atoms, the regular pair list contains all pairs including hydrogens
void topology::create_pair_list(bool pbc, double halfboxsize, double boxsize, double listcutoff, std::vector<atom_nb_entry> * pair_list,  std::vector<atom_nb_entry> * solv_list,double * coords)
{
    int iatom, jatom, ires, jres, k, ientry;
    long int oldsize;
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
    //atom pairs are not to be on the pair list unless at least one residue is in the all-atom region
    oldsize=pair_list->capacity();
    pair_list->clear();
    pair_list->reserve(oldsize);//restore old capacity to avoid having to extend.
    oldsize=solv_list->capacity();
    solv_list->clear();
    solv_list->reserve(oldsize);
    for (ires=0; ires<nres; ires++) for (jres=ires; jres<nres; jres++) {
        if (jres==ires) {
            d=0;
        } else { //jres>ires
            d=sqrt(pbc_distance2(pbc,halfboxsize,boxsize,&res_com[3*ires],&res_com[3*jres]));
            d=d-res_radius[ires]-res_radius[jres];
        }
        if (d<listcutoff) {
            piece_of_list=&orig_pair_list_by_res[ires*nres+jres];
            //the main pair list contains every atom pair that is within cutoffs and in which at least one atom is in the AA region
            if (aaregion_res[ires] || aaregion_res[jres]) pair_list->insert(pair_list->end(),piece_of_list->begin(),piece_of_list->end());
            //the solvation list contains all pairs of non-hydrogen atoms that are within cutoffs, regardless of the AA region
            piece_of_list=&orig_solv_list_by_res[ires*nres+jres];
            solv_list->insert(solv_list->end(),piece_of_list->begin(),piece_of_list->end());
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

/*bool simulation::update_pair_list_if_needed(long int istep, double * coords)
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
}*/
