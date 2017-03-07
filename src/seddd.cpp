#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <vector>
#include "ffield.h"
#include "util.h"
#include "rotations.h"
#include "topology.h"
#include "seddd.h"

#define SCALE_FACTOR 2.0/(4.0*M_PI*sqrt(M_PI))
void topology::calculate_solvation_volumes(seddd_params * params, double cutoff2, std::vector<atom_nb_entry> * solv_list, double * coords, double * frac_volumes, forcefield * ffield)
{
    double * group_com;
    double * group_mass;
    double mass;
    double r,r2, voli, volj, ri, rj, lambdai, lambdaj,aux,dvi, dvj;
    int i,iatom, jatom, iclass, jclass, igroup, jgroup,k;
#ifdef TIMERS
    switch_timer(TIMER_CALC_VOLUMES);
#endif
    //compute the COM of each group -- not sure if this is how it should be done
    /*group_com=(double *) checkalloc(3*ngroup,sizeof(double));
    group_mass=(double *) checkalloc(ngroup,sizeof(double));
    for (igroup=0; igroup<ngroup; igroup++) {
        for (k=0; k<3; k++) group_com[3*igroup+k]=0.0;
    for (iatom=0; iatom<natom; iatom++) {
        igroup=atoms[iatom].group;
        mass=atoms[iatom].mass;
        for (k=0; k<3; k++) group_com[]
        }
    }*/
    //
    for (iatom=0; iatom<natom; iatom++) frac_volumes[iatom]=0.0;
    //We've constructed the solvation list so it includes only heavy atom pairs potentially within cutoff.
    for (i=0; i<solv_list->size(); i++) {
        iatom=(*solv_list)[i].iatom;
        jatom=(*solv_list)[i].jatom;
        //both atoms are guaranteed to be heavy atoms by construction
        //if ((atoms[iatom].mass<2.0) || (atoms[jatom].mass<2.0)) continue; //only for pairs of heavy atoms
        r2=0.0;
        for (k=0; k<3; k++) r2+=(coords[3*jatom+k]-coords[3*iatom+k])*(coords[3*jatom+k]-coords[3*iatom+k]);
        if (r2>cutoff2) continue;
        iclass=atoms[iatom].classx;
        jclass=atoms[jatom].classx;
        ri=ffield->vdwParams[iclass].sigma;
        rj=ffield->vdwParams[jclass].sigma;
        voli=params->hydration_volume[iclass];
        volj=params->hydration_volume[jclass];
        lambdai=params->hydration_shell_thickness[iclass];
        lambdaj=params->hydration_shell_thickness[jclass];
        r=sqrt(r2);
        //apply eq. 3 twice, once for i and once for j
        aux=(r-ri)/lambdai;
        aux=exp(-aux*aux);
        dvi=SCALE_FACTOR*volj*aux/(lambdai*r2);
        frac_volumes[iatom]+=dvi;
        /*if (frac_volumes[iatom]>1.0) {
            printf("warning: too high fractional solvation for atom %s %d %s = %.10f\n",atoms[iatom].resName,atoms[iatom].resNum+1,atoms[iatom].name,frac_volumes[iatom]);
        }*/
        aux=(r-rj)/lambdaj;
        aux=exp(-aux*aux);
        dvj=SCALE_FACTOR*voli*aux/(lambdaj*r2);
        frac_volumes[jatom]+=dvj;
        /*if (frac_volumes[jatom]>1.0) {
                printf("warning: too high fractional solvation for atom %s %d %s = %.10f\n",atoms[jatom].resName,atoms[jatom].resNum+1,atoms[jatom].name,frac_volumes[jatom]);
        }*/
    }
    //The fractional volume can be above 1.0 if there are very closely spaced pairs of atoms (e.g. <0.5 A).  While this is implausible we do have to deal with it during docking preparation.
    //The most sensible thing to do is to restrict each fractional volume to be less than 1
    for (iatom=0; iatom<natom; iatom++) if (frac_volumes[iatom]>1.0) frac_volumes[iatom]=1.0;
    //now fill in for the hydrogen atoms
    for (iatom=0; iatom<natom; iatom++) if (atoms[iatom].mass<2.0) {
        jatom=atoms[iatom].bondedAtomList[0];
        frac_volumes[iatom]=frac_volumes[jatom];
    }
#ifdef TIMERS
    switch_timer(TIMER_OTHER);
#endif
}

void topology::check_solvation_params(seddd_params * params)
{
    int iatom, iclass;
    for (iatom=0; iatom<natom; iatom++) if (atoms[iatom].atomicNum>1) {
        iclass=atoms[iatom].classx;
        if (!params->read[iclass]) {
           printf("Error: solvation parameters for atom %s %d %s class %d have not been read.\n",atoms[iatom].resName,atoms[iatom].resNum+1,atoms[iatom].name,iclass);
           die();
        }
    }
}

void read_solvation_params(char * line, seddd_params * params)
{
    char fname[255],buf[255],clname[255];
    char * p;
    FILE * input;
    int iclass;
    double vol, thickness;
    for (iclass=0; iclass<MAX_NUM_OF_ATOM_CLASSES; iclass++) params->read[iclass]=false;
    sscanf(line,"%lg %lg %lg %lg %s",&params->eps0,&params->eps1,&params->c, &params->frac_vol_tol, fname);
    params->delta_eps=params->eps1-params->eps0;
    printf("Garden and Zhorov SEDDD method enabled.  Dielectric constant given by epsilon = r * (eps_0 + (1 - s_kl) * (eps_1 - eps_0))\n");
    printf("where s_kl = c * (v_k + v_l)\n");
    printf("eps_0, eps_1                        = %.1f %.1f\n",params->eps0,params->eps1);
    printf("c                                   = %.3f\n",params->c);
    printf("tolerance for changed frac. volumes = %g\n",params->frac_vol_tol); 
    input=fopen(fname,"r");
    if (input==NULL) {
        printf("Could not open file %s\n",fname);
        die();
    }
    while (true) {
        fgets(buf,sizeof(buf),input);
        if (feof(input)) break;
	//terminate the string at the first '#' in order to ignore comments
	p=strchr(buf,'#');
	if (p!=NULL) *p='\0';
        //skip over white space at the beginning, skip if empty
	p=&buf[0];
	while (isspace(*p)) p++;
	if (*p=='\0') continue;
        //format: atom class number, hydration volume, hydration shell thickness
        sscanf(p,"%d %s %lg %lg\n",&iclass,clname,&vol,&thickness);
        params->read[iclass]=true;
        params->hydration_volume[iclass]=vol;
        params->hydration_shell_thickness[iclass]=thickness;
    }
}
