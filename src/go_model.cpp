#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "go_model.h"
#include "util.h"
//#include "fragments.h"
#include "topology.h"
#include "rotations.h" //for pbc_distance2

//#define ALPHA_CARBON "CA"
#define CLASH_ENERGY 1e20

//i/o support routines
//format:  go_hardcore go_cutoff m n native_energy nonnative_energy
//m and n must be even and positive, with m>n
void read_go_params(char * line, go_model_params * params)
{
    sscanf(line,"%lg %lg %lg %d %d %lg %lg",&params->hardcore,&params->cutoff,&params->subcutoff,&params->m,&params->n,&params->native_energy,&params->nonnative_energy);
    if ((params->m%2!=0) || (params->m<params->n) || (params->n<0) || (params->n%2!=0)) {
        printf("Invalid set of Go parameters.\n");
        die();
    }
    finish_go_params(params);
}

/*bool read_go_params2(char * token1, char * token2)
{
    result=true;
    if (strncmp(*/
bool read_go_parameter(char * word, char * rest_of_line, go_model_params * params)
{
    if (strcasecmp("GO_HARDCORE",word)==0) {
        sscanf(rest_of_line,"%lg",&params->hardcore);
    } else if (strcasecmp("GO_NATIVE",word)==0) {
        sscanf(rest_of_line,"%lg",&params->cutoff);
    } else if (strcasecmp("GO_EXPONENTS",word)==0) {
        sscanf(rest_of_line,"%d %d",&params->m, &params->n);
    } else if (strcasecmp("GO_WELLDEPTH",word)==0) {
        sscanf(rest_of_line,"%lg",&params->native_energy);
        params->nonnative_energy=params->native_energy;
    } else {
        return false;
    }
    return true;
}

void finish_go_params(go_model_params * params)
{
    params->subcutoff=NAN;
    if (params->nonnative_energy==0) params->nonnative_energy=params->native_energy;
    params->hardcore2=4*params->hardcore*params->hardcore; //hardcore distance = 2 * hardcore radius
    params->cutoff2=params->cutoff*params->cutoff;
    if (params->n>0) params->ratio=((double) params->m)/(params->n); else params->ratio=0;
    params->scaled_native_en=-params->native_energy/(1-params->ratio);
    params->scaled_nonnative_en=-params->nonnative_energy/(1-params->ratio);
    params->hn=params->n/2;
    params->hm=params->m/2;
    params->rsubcutoff=0; //eliminate subcutoff
}



//this can be called from the setup or from table::print_header_info
void print_go_params(go_model_params params)
{
    printf("Go model hardcore and cutoff distances:           %.2f %.2f A\n",params.hardcore,params.cutoff);
    printf("Go model subcutoff (units of native distance):    %.2f\n",params.subcutoff);
    //printf("Go model delta:                                   %.2f A\n",hdr.go_delta);
    printf("Go model exponents:                               %d %d\n",params.m,params.n);
    printf("Go model native and nonnative energies:           %.2f %.2f kcal/mol\n",params.native_energy,params.nonnative_energy);
    if ((params.m%2!=0) || (params.m<params.n) || (params.n<0) || (params.n%2!=0)) {
        printf("Invalid set of Go parameters.\n");
        die();
    }
}



go_model_info::go_model_info()
{
    nentries=0;
    nnative=0;
    //entries=NULL;
    distances=NULL;
    ifragname[0]='\0';
    jfragname[0]='\0';
}



go_model_info::~go_model_info()
{
    //if (entries!=NULL) free(entries);
    if (distances!=NULL) free(distances);
}

/*int go_model_info::find_entry(int ires, int jres) {
    go_model_entry entry;
    go_model_entry * p;
    //int low, high, mid,comp;
    entry.ires=ires;
    entry.jres=jres;
    if (entries==NULL) return -1; //no entries? not there
    p=(go_model_entry *) bsearch(&entry,entries,nentries,sizeof(go_model_entry),compare_entries);
    if (p==NULL) return -1; else return p-entries; //index
}*/

//utility routine to convert fragment atoms to residues and find the appropriate position in the distance table
/*int go_model_info::get_index(int nres, ires, jres)
{
    int ires,jres,swap,index;
    ires=ifragtype->atoms[ifragatom].res;
    jres=jfragtype->atoms[jfragatom].res;

    index=(inres*(jres-jfragtype->startres))+(ires-ifragtype->startres);
    return index;
}*/




//to prepare lists of atoms for making Go models


void go_model_info::create_contact_map(int nres, reslookup * resinfo, int natom, double * nativecoords, go_model_params * params, subset aaregion_res)
{
    int iatom,jatom,ires,jres, index,k;
    int * fraglist;
    double distance2;
    nentries=0;

    distances=(double *) checkalloc(nres*nres,sizeof(double));
    nnative=0;
    //printf("Setting up Go model with hard-core radius %g A, cutoff %g A \n",params->hardcore,params->cutoff);

    //for (index=0; index<nres*nres; index++) distances[index]=1e20;
    for (ires=0; ires<nres; ires++)
        for (jres=ires+2; jres<nres; jres++)
            //Go energy only if both atoms are in the CG region.
            if (!aaregion_res[ires] && !aaregion_res[jres]) {
                    //Determine native distance.
                    iatom=resinfo[ires].branchatom;
                    jatom=resinfo[jres].branchatom;

                    distance2=0;
                    for (k=0; k<3; k++) distance2+=(nativecoords[3*iatom+k]-nativecoords[3*jatom+k])*(nativecoords[3*iatom+k]-nativecoords[3*jatom+k]);
                    //if ((((ifragatom==24) && (jfragatom==127)) || ((ifragatom==127) && (jfragatom==24))) && (distance2<2000)) printf("***: %d %d %d %d %.16f\n",ifrag,ifragatom,jfrag,jfragatom,sqrt(distance2));
                    /*counter++;
                    if (counter%10000==0) printf("Calculated %d distances\n",counter);*/
#ifdef DEBUG_NON_TABULATED
                    printf("Adding Go model term for residues %d %d atoms %d %d distance %.2f A\n",ires,jres,iatom,jatom,sqrt(distance2));
#endif
                    index=nres*ires+jres;
                    distances[index]=distance2;
                    index=nres*jres+ires; //just in case
                    distances[index]=distance2;
                    if (distance2<params->cutoff2) nnative++;
                }

}

double go_model_info::moved_energy(int pbc, double halfboxsize, double boxsize, go_model_params * params, int nres, reslookup * resinfo, subset& aaregion_res, subset& movedatoms, int natom, double * coords)
{
    double en,entot,dx[3],r2;
    int ientry,k,ifragatom,jfragatom,iexp,ires,jres,iatom,jatom;
    double a,am,an,native_distance2;
    entot=0.0;

     for (ires=0; ires<nres; ires++)
        for (jres=ires+2; jres<nres; jres++)
        //Go energy only if both atoms are in the CG region.
            if (!aaregion_res[ires] && !aaregion_res[jres]) {

                    //Determine native distance.
                iatom=resinfo[ires].branchatom;
                jatom=resinfo[jres].branchatom;
                //only do energy if one atom moved and not the other
                if (!(movedatoms[iatom]^movedatoms[jatom])) continue;
                r2=pbc_distance2(pbc,halfboxsize,boxsize,&coords[3*iatom],&coords[3*jatom]);
                if (r2<=0) {
                    printf("go_model_info::energy: error, ifragatom=%d jfragatom=%d r2=%.2f\n",ifragatom,jfragatom,r2);
                    die();
                }
            /*if (r2<entries[ientry].low_distance2) return CLASH_ENERGY;
            if (r2<entries[ientry].hi_distance2) { //intermediate range in equations
                if (entries[ientry].native) en+=native_energy; else en+=nonnative_energy;
            }*/
                native_distance2=distances[nres*ires+jres];
                if (native_distance2<params->cutoff2) a=native_distance2/r2; else a=params->hardcore2/r2;
                //if (a<params->rsubcutoff)  continue; //<1% of native energy
                an=1;
                for (iexp=1; iexp<=params->hn; iexp++) an*=a;
                //an now equals (r/r_nat)^n or (r/r_hc)^n
                am=an;
                if (params->m==2*params->n) am=an*an; else for (; iexp<=params->hm; iexp++) am*=a;
                //am now equals (r/r_nat)^m or (r/r_hc)^m
                if (native_distance2<params->cutoff2) {
                    en=params->scaled_native_en*(am-params->ratio*an);
                } else {
                    en=params->scaled_nonnative_en*am;
                }
#ifdef DEBUG
                if (en!=0) printf("Go model: %d %d %d %d %.16f %.16f %.16f %.16f\n",ires,jres,iatom,jatom,sqrt(r2),sqrt(native_distance2),a,en);
#endif
                entot+=en;
        }
    return entot;
}

double go_model_info::energy(int pbc, double halfboxsize, double boxsize, go_model_params * params, int nres, reslookup * resinfo, subset& aaregion_res, int natom, double * coords)
{
    double en,entot,dx[3],r2;
    int ientry,k,ifragatom,jfragatom,iexp,ires,jres,iatom,jatom;
    double a,am,an,native_distance2;
    entot=0.0;

     for (ires=0; ires<nres; ires++)
        for (jres=ires+2; jres<nres; jres++)
        //Go energy only if both atoms are in the CG region.
            if (!aaregion_res[ires] && !aaregion_res[jres]) {

                    //Determine native distance.
                iatom=resinfo[ires].branchatom;
                jatom=resinfo[jres].branchatom;

                r2=pbc_distance2(pbc,halfboxsize,boxsize,&coords[3*iatom],&coords[3*jatom]);
                if (r2<=0) {
                    printf("go_model_info::energy: error, ifragatom=%d jfragatom=%d r2=%.2f\n",ifragatom,jfragatom,r2);
                    die();
                }
            /*if (r2<entries[ientry].low_distance2) return CLASH_ENERGY;
            if (r2<entries[ientry].hi_distance2) { //intermediate range in equations
                if (entries[ientry].native) en+=native_energy; else en+=nonnative_energy;
            }*/
                native_distance2=distances[nres*ires+jres];
                if (native_distance2<params->cutoff2) a=native_distance2/r2; else a=params->hardcore2/r2;
                if (a<params->rsubcutoff)  continue; //<1% of native energy
                an=1;
                for (iexp=1; iexp<=params->hn; iexp++) an*=a;
                //an now equals (r/r_nat)^n or (r/r_hc)^n
                am=an;
                if (params->m==2*params->n) am=an*an; else for (; iexp<=params->hm; iexp++) am*=a;
                //am now equals (r/r_nat)^m or (r/r_hc)^m
                if (native_distance2<params->cutoff2) {
                    en=params->scaled_native_en*(am-params->ratio*an);
                } else {
                    en=params->scaled_nonnative_en*am;
                }
#ifdef DEBUG
                if (en!=0) printf("Go model: %d %d %d %d %.16f %.16f %.16f %.16f\n",ires,jres,iatom,jatom,sqrt(r2),sqrt(native_distance2),a,en);
#endif
                entot+=en;
        }
    return entot;
}

/*i
/*int go_model_info::count_native_contacts(int pbc, double halfboxsize, double boxsize, double cutoff, double ratio, fragmenttype * ifragtype, double * icoords, fragmenttype * jfragtype, double * jcoords)
{
    double dx[3],r2,cutoff2,ratio2;
    int ientry,k,ifragatom,jfragatom,iexp;
    double a,am,an,native_distance2;
    int ncontact;
    ncontact=0;
    cutoff2=cutoff*cutoff;
    ratio2=ratio*ratio;

    for (ifragatom=0; ifragatom<ifragtype->natom; ifragatom++)
        for (jfragatom=0; jfragatom<jfragtype->natom; jfragatom++) {
            native_distance2=distances[get_index(ifragtype,ifragatom,jfragtype,jfragatom)];
            if (native_distance2<=cutoff2) {
                r2=pbc_distance2(pbc,halfboxsize,boxsize,&icoords[3*ifragatom],&jcoords[3*jfragatom]);
                if (r2<=0) {
                    printf("go_model_info::energy: error, ifragatom=%d jfragatom=%d r2=%.2f\n",ifragatom,jfragatom,r2);
                    die();
                }
                a=r2/native_distance2;
                if (a<=ratio2) ncontact++;
            }
    }
    return ncontact;
}*/

//this restores to the same state as after create_contact_map.  need to call set_parameters to finish the go model before use.
/*void go_model_info::read_file(char * fname)
{
    FILE * input;
    int ires,jres,index;
    double r;
    input=fopen(fname,"r");
    if (input==NULL) {
        printf("Could not open file %s.\n",fname);
        die();
    }
    printf("Reading Go model information from file %s.\n",fname);
    fscanf(input,"%s %s %d %d\n",ifragname,jfragname,&inres,&jnres);
    distances=(double *) checkalloc(inres*jnres,sizeof(double));
    for (index=0; index<inres*jnres; index++) distances[index]=1.0e20;
    //fscanf(input,"%d\n",&nentries);
    while (!feof(input)) {
        fscanf(input,"%d %d %lg\n",&ires,&jres,&r);
        index=inres*(jres-1)+(ires-1);
        distances[index]=r*r;
    }
    fclose(input);
}

void go_model_info::write_file(char * fname)
{
    FILE * output;
    int ires,jres,index;
    output=fopen(fname,"w");
    if (output==NULL) {
        printf("Could not open file %s.\n",fname);
        die();
    }
    fprintf(output,"%s %s %d %d\n",ifragname,jfragname,inres,jnres);
    //fprintf(output,"%d\n",nentries);
    for (ires=1; ires<=inres; ires++) for (jres=1; jres<=jnres; jres++) {
        index=inres*(jres-1)+(ires-1);
        fprintf(output,"%d %d %.16f\n",ires,jres,sqrt(distances[index]));
    }
    fclose(output);
    printf("Go mod*
}*/
