#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
//#include "tables.h"
#include "mc.h"
#include "mt.h"
#include "rotations.h"
#include "ffield.h"
#include "go_model.h"
#include "util.h"
#include "topology.h"

void simulation::read_lambda_schedule(char * fname)
{
    FILE * input;
    char buffer[255];
    int isched;
    input=fopen(fname,"r");
    if (input==NULL) {
        printf("Could not open file %s for reading\n",fname);
        die();
    }
    printf("Reading lambda schedule from file %s.\n",fname);
    n_lambda_schedule = 0;
    while (!feof(input)) {
        fgets(buffer,sizeof(buffer),input);
        if (strlen(buffer)==0) continue;
        lambda_schedule = (lambda_sched_info *) checkrealloc(lambda_schedule,n_lambda_schedule+1,sizeof(lambda_sched_info));
        sscanf(buffer,"%ld %lg %lg",&lambda_schedule[n_lambda_schedule].step,
            &lambda_schedule[n_lambda_schedule].lambda_vdw,&lambda_schedule[n_lambda_schedule].lambda_elec);
        lambda_schedule[n_lambda_schedule].step %= nsteps_temper_move;
        n_lambda_schedule++;
    }
    n_lambda_schedule--;
    printf("Lambda schedule:\n");
    printf("Step number  VDW lambda  Elec lambda\n");
    for (isched=0; isched<n_lambda_schedule; isched++) printf("%d %.4f %.4f\n",
        lambda_schedule[isched].step,lambda_schedule[isched].lambda_vdw,lambda_schedule[isched].lambda_elec);
}

void simulation::adjust_lambdas_and_accumulate_work(long int istep, bool * lambda_changed)
{
    int isched, found_sched;
    long int reduced_step;
    double new_lambda_vdw, new_lambda_elec, delta_work;
    double internal_energies[EN_TERMS],intxn_energies[EN_TERMS];
    reduced_step = (istep % nsteps_temper_move);
    //search the lambda schedule for a corresponding entry
    found_sched=-1;
    for (isched=0; isched<n_lambda_schedule; isched++) if (lambda_schedule[isched].step==reduced_step) {
        found_sched=isched;
        break;
    }
    //if not on the schedule, we don't need to do anything
    if (found_sched<0) {
        *lambda_changed=false;
        return;
    }
    new_lambda_vdw=lambda_schedule[found_sched].lambda_vdw;
    new_lambda_elec=lambda_schedule[found_sched].lambda_elec;
    //compute interaction energy between ligand and protein
    //avoid recomputing solvation volumes or pair lists
    ffield->subset_energy(&solvation_params,cutoff2,top->natom,top->atoms,top->ligand,new_pair_list.size(),&new_pair_list[0],newcoords,new_frac_volumes,internal_energies,intxn_energies);
    delta_work=(new_lambda_vdw-current_lambda_vdw)*intxn_energies[EN_VDW_EXACT]+(new_lambda_elec-current_lambda_elec)*intxn_energies[EN_ELEC_EXACT];
    current_noneq_work+=delta_work;
    printf("Step %ld: Changing lambdas to VDW %.4f and elec %.4f Delta work %.6f\n",istep, new_lambda_vdw, new_lambda_elec, delta_work);
    current_lambda_vdw=new_lambda_vdw;
    current_lambda_elec=new_lambda_elec;
    *lambda_changed=true;
}

void simulation::perform_ncmc_move(void)
{
    double r,p;
    bool accepted;
    int i;
    natt_ncmc++;
    if (current_noneq_work<0) p=1.0; else p=exp(-beta*current_noneq_work);
    sumprob_ncmc+=p;
    r=genrand_real3();
    accepted=(r<p);
    if (accepted) {
        //NCMC move accepted.
        nacc_ncmc++;
        for (i=0; i<3*top->natom; i++) oldcoords_ncmc[i]=newcoords[i];
    } else {
        //NCMC move rejected.
        for (i=0; i<3*top->natom; i++) newcoords[i]=oldcoords_ncmc[i];
        //ensure that the pair list, solvation list, and fractional volumes are updated to reflect the new (actually old) coordiantes
        if (use_nb_list && !go_only) top->create_pair_list(pbc,halfboxsize,boxsize,listcutoff,&new_pair_list,&new_solv_list,newcoords);
#ifdef SEDDD
        if (!go_only) top->calculate_solvation_volumes(&solvation_params,cutoff2,&new_solv_list,newcoords,new_frac_volumes,ffield);
#endif
        //ensure that "old" and "new" coordinates for regular MC match.
        for (i=0; i<3*top->natom; i++) oldcoords[i]=newcoords[i];
        old_pair_list=new_pair_list;
        old_solv_list=new_solv_list;
        for (i=0; i<top->natom; i++) old_frac_volumes[i]=new_frac_volumes[i];
    }
    printf("NCMC moves   Attempted %ld   Accepted %ld   Avg. probability = %g\n",natt_ncmc,nacc_ncmc,sumprob_ncmc/natt_ncmc);
    //mcloop will reset current_noneq_work to zero.
}
