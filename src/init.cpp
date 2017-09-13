#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
//#include "tables.h"
#include "mc.h"
#include "mt.h"
#include <math.h>
#include "rotations.h"
#include "ffield.h"
#include "go_model.h"
#include "util.h"
#include "topology.h"
//#include "fragments.h"
#if defined(__unix__) || defined(__unix) || \
        (defined(__APPLE__) && defined(__MACH__))
#define UNIX
#include <sys/types.h>
#include <unistd.h>
#endif

#if defined(PARALLEL) || defined(EXCHANGE)
simulation::simulation(int _mynod, int _numnod)
#else
simulation::simulation(void)
#endif
{
    FILE * f;
    unsigned long seed;
    double listcutoff;
    double ptot,mctemp,p,size;
    double total_table_size;
    char buffer[255],word[255],tablefmt[255],covtablefmt[255],fragfmt[255],structf3name[255],psffname[255];
    bool start,restart;
    int itype,ifrag,i,jtype,move;
    double energies[EN_TERMS],etot;
    memset(this,0,sizeof(*this));
#if defined(PARALLEL) || defined(EXCHANGE)
    mynod=_mynod;
    numnod=_numnod;
#endif
    top=NULL;
    ffield=NULL;
    sequence=NULL;
    memset(ligand_resname,0,sizeof(ligand_resname));
    initcoords=NULL;
    oldcoords=NULL;
    pbc=false;
    aaregion_specified=false;
    initialized=false;
    go_only=false;
    seed = 0;
    do_ncmc = false;
    n_lambda_schedule = 0;
    lambda_schedule = NULL;
    current_lambda_vdw=1.0;
    current_lambda_elec=1.0;
    ncmc_write_frames = false;
}




simulation::~simulation()
{
    int i;
    /*if (tables!=NULL) {
        for (i=0; i<top->nfragtypes*top->nfragtypes; i++)
        	if (tables[i]!=NULL) delete tables[i];
    }
    if (en_by_table!=NULL) free(en_by_table);
    free(tables);*/
    free(sequence);
    free(initcoords);
    free(oldcoords);
    free(newcoords);
#ifdef SEDDD
    free(old_frac_volumes);
    free(new_frac_volumes);
#endif
    /*free(oldcenter);
    free(oldorient);
    free(newcenter);
    free(neworient);*/
    //if (frag_nblist!=NULL) delete frag_nblist;
    if (acc_by_atom!=NULL) free(acc_by_atom);
    if (att_by_atom!=NULL) free(att_by_atom);
    delete ffield;
    delete go_model;
    delete top;
}

void simulation::process_commands(char * infname)
{
    char command[255],command2[255],fname[255],fmt[255],word[255],word2[255],fname2[255],unit[255];
    char * token;
    const char * delim = " \t\n";
    char chain;
    FILE * input;
    FILE * output;
    double mctemp, trans_size, rot_size, bond_rot_size, desired_ligand_com[3],temp;
    int i,iatom,move,nsearch;

    double energies[EN_TERMS],etot;
    double intxn_energies[EN_TERMS],internal_energies[EN_TERMS],total_internal_energy,total_intxn_energy;
    mctemp=300.0;
    nres=0;
    initialized = false;
    input=fopen(infname,"r");
    if (input==NULL) {
        printf("Could not open command file %s.\n",infname);
        die();
    }
    while (true) {
        fgets(command,sizeof(command),input);
        if (feof(input)) break; //end of file
        for (i=0; i<strlen(command); i++) if (command[i]=='#') command[i]='\0';
        printf("Processing command: %s\n",command);
        //for (i=0; i<strlen(command); i++) command[i]=toupper(command[i]);
        strncpy(command2,command,sizeof(command2));
        token=strtok(command2,delim);
        if (token==NULL) continue; //blank line
        //else if (*token=='#') continue; //comment
        else if (strcasecmp("END",token)==0) return; //end of commands
        else if (strcasecmp("READ",token)==0) { //read a PDB file
            token=strtok(NULL,delim);
            strncpy(fmt,token,sizeof(fmt));
            token=strtok(NULL,delim);
            strncpy(fname,token,sizeof(fname));
            if (strcasecmp("DEFS",fmt)==0) {
                if (top!=NULL) delete top;
                if (ffield==NULL) {
                    printf("You need to read a forcefield first.\n");
                    die();
                }
                printf("Loading definitions file %s\n",fname);
                top = new topology(fname,ffield);
            } else if (strcasecmp("FFIELD",fmt)==0) {
                if (ffield!=NULL) delete ffield;
                printf("Loading forcefield file %s\n",fname);
                ffield = new forcefield(fname);
            } else if (strcasecmp("PDB",fmt)==0) {
                strncpy(structfname,fname,sizeof(structfname));
            } else {
                printf("Unrecognized format.\n");
                die();
            }
        } else if (strcasecmp("WRITE",token)==0) { //Fill in coordinates and write a pdb file.
            token=strtok(NULL,delim);
            strncpy(fmt,token,sizeof(fmt));
            token=strtok(NULL,delim);
            strncpy(fname,token,sizeof(fname));
            if (strcasecmp("PDB",fmt)==0) {
                top->calculate_solvation_volumes(&solvation_params,cutoff2,&old_solv_list,initcoords,old_frac_volumes,ffield);
                top->write_pdb_file(fname,initcoords,old_frac_volumes); //if old_frac_volumes is not allocated, will write "0.0" in column.
            } else if (strcasecmp("PSF",fmt)==0) {
                top->write_psf_file(fname,ffield);
            /*} else if (strcasecmp("REST",fmt)==0) {
#ifdef EXCHANGE
                if (strstr(fname,"%d")!=NULL) { //fill in replica number
                    strncpy(fmt,fname,sizeof(fmt));
                    snprintf(fname,sizeof(fname),fmt,myrep+1);
                }
#endif // EXCHANGE
                write_restart(0,fname);*/
            } else if (strcasecmp("MCBYATOM",fmt)==0) {
                //write info on the counts of moves attempted and accepted by atom
                output=fopen(fname,"w");
                if (output==NULL) {
                   printf("Could not open file %s.\n",fname);
                   die();
                }
                for (iatom=0; iatom<top->natom; iatom++) fprintf(output,"%d %s %d %s %ld %ld\n",iatom+1,
                   top->atoms[iatom].resName,top->atoms[iatom].resNum+1,top->atoms[iatom].name,att_by_atom[iatom],acc_by_atom[iatom]);
                fclose(output);
                printf("MC by-atom counts written to file %s.\n",fname);
            } else {
                printf("Unrecognized format.\n");
                die();
            }
        } else if (strcasecmp("INSERT",token)==0) {
            token=strtok(NULL,delim);
            strncpy(word,token,sizeof(word));
            if (strcasecmp("SEQUENCE",word)==0) {
                sequence=read_multiline(input);
                strupper(sequence,strlen(sequence));
                nres=count_words(sequence);
                if (nres==0) {
                    printf("You need to specify a sequence first.\n");
                    die();
                }
            } else if (strcasecmp("LIGAND",word)==0) {
                token=strtok(NULL,delim);
                strncpy(ligand_resname,token,sizeof(ligand_resname));
                strupper(ligand_resname,sizeof(ligand_resname));
                top->ligand_res=nres;
                nres++;
            } else {
                printf("Unrecognized command.\n");
                die();
            }
        } else if (strcasecmp("AAREGION",token)==0) {
            if ((top==NULL) || (nres==0)) {
                printf("You must load a definitions file and insert a sequence first.\n");
                die();
            }
            top->aaregion_res.init(nres);
            token=strtok(NULL,delim);
            strncpy(word,token,sizeof(word));
            top->aaregion_res.parse_int_list(word);
            if (top->ligand_res>=0) top->aaregion_res+=top->ligand_res; //force inclusion of ligand in AA region for docking
            aaregion_specified = true;
            create_lists();
            go_only=true;
            for (i=0; i<top->nres; i++) if (top->aaregion_res[i]) {
                go_only=false;
                break;
            }
        } else if (strcasecmp("GO_MODEL",token)==0) {
            read_go_model_info(input);
        } else if (strcasecmp("MOVES",token)==0) {
            //in mcmoves.cpp
            read_move_info(input);
        } else if (strcasecmp("BOXSIZE",token)==0) {
            pbc=true;
            token=strtok(NULL,delim);
            strncpy(word,token,sizeof(word));
            sscanf(word,"%lg",&boxsize);
            halfboxsize=0.5*boxsize;
        } else if (strncasecmp("CUTOFF",token,6)==0) {
            token+=strlen(token)+1;
            strncpy(word,token,sizeof(word));
            sscanf(word,"%lg %lg",&cutoff,&listcutoff);
            cutoff2=cutoff*cutoff;
#ifdef SEDDD
        } else if (strcasecmp("SEDDD",token)==0) {
            token+=strlen(token)+1;
            strncpy(word,token,sizeof(word));
            read_solvation_params(word,&solvation_params);
#else
        } else if (strcasecmp("EPS",token)==0) {
            token=strtok(NULL,delim);
            strncpy(word,token,sizeof(word));
            sscanf(word,"%lg",&eps);
#endif
        //This is temporary until the Mobley hamiltonian tempering method is fully implemented.
        } else if (strcasecmp("LAMBDA",token)==0) {
            token+=strlen(token)+1;
            strncpy(word,token,sizeof(word));
            sscanf(word,"%lg %lg",&current_lambda_vdw,&current_lambda_elec);
        } else if (strcasecmp("NCMC",token)==0) {
            //NCMC move_length block_length lambda_schedule_file [WRITE_FRAMES_IN_MOVE]
            do_ncmc=true;
            printf("Nonequilibriumm Candidate Monte Carlo enabled.\n");
            token+=strlen(token)+1;
            strncpy(word,token,sizeof(word));
            ncmc_write_frames=((strstr(word,"WRITE_FRAMES_IN_MOVE")!=NULL) || (strstr(word,"write_frames_in_move")!=NULL));
            sscanf(word,"%ld %ld %s",&nsteps_temper_move,&nsteps_block,fname);
            if (nsteps_block<nsteps_temper_move) {
               printf("Number of steps in block must exceed number in NCMC move.\n");
               die();
            }
            printf("Each block of %ld moves will include one NCMC move of %ld moves and %ld additional regular MC moves.\n",
                nsteps_block,nsteps_temper_move,nsteps_block-nsteps_temper_move);
            if (ncmc_write_frames) {
                printf("Frames within NCMC moves will be written to the trajectory.\n");
            } else {
                printf("Frames within NCMC moves will not be written to the trajectory.\n");
            }
            read_lambda_schedule(fname);
        } else if (strcasecmp("SEED",token)==0) {
            token=strtok(NULL,delim);
            strncpy(word,token,sizeof(word));
            sscanf(word,"%lu",&seed);
        } else if (strncasecmp("TEMP",token,4)==0) {
            token=strtok(NULL,delim);
            strncpy(word,token,sizeof(word));
            sscanf(word,"%lg",&mctemp);
            beta=1/(KBOLTZ*mctemp);
        } else if (strcasecmp("SAVEFREQ",token)==0) {
            token=strtok(NULL,delim);
            strncpy(word,token,sizeof(word));
            sscanf(word,"%d",&nsave);
        } else if (strcasecmp("PRINTFREQ",token)==0) {
            token=strtok(NULL,delim);
            strncpy(word,token,sizeof(word));
            sscanf(word,"%d",&nprint);
        } else if (strncasecmp("TRAJ",token,4)==0) {
            token=strtok(NULL,delim);
            strncpy(trajfmt,token,sizeof(trajfmt));
            token=strtok(NULL,delim);
            strncpy(xyzfname,token,sizeof(xyzfname));
        } else if (strncasecmp("MCLOG",token,4)==0) {
            token=strtok(NULL,delim);
            strncpy(mclogfname,token,sizeof(mclogfname));
        } else if (strncasecmp("INIT",token,4)==0) {
            if (!initialized) finish_initialization();
        } else if (strcasecmp("ENERGY",token)==0) {
            if (!initialized) finish_initialization();
            if (oldcoords==NULL) {
                oldcoords=(double *) checkalloc(3*top->natom,sizeof(double));
                for (i=0; i<3*top->natom; i++) oldcoords[i]=initcoords[i];
            }
            token=strtok(NULL,delim);
            strncpy(word,token,sizeof(word));
            if (strcasecmp(word,"TOTAL")==0) {
                top->create_pair_list(pbc,halfboxsize,boxsize,listcutoff,&old_pair_list,&old_solv_list,initcoords);
#ifdef SEDDD
                top->calculate_solvation_volumes(&solvation_params,cutoff2,&old_solv_list,initcoords,old_frac_volumes,ffield);
                total_energy(initcoords,&old_pair_list,old_frac_volumes,energies,&etot);
#else
                total_energy(initcoords,&old_pair_list,energies,&etot);
#endif
                print_energies(stdout,true,"Energy:",0,energies,etot);
            } else if (strcasecmp(word,"LIGAND")==0) {
                ligand_energies(oldcoords,&total_internal_energy,internal_energies,&total_intxn_energy,intxn_energies);
                print_energies(stdout,true,"Internal:",0,internal_energies,total_internal_energy);
                print_energies(stdout,true,"Intxn:",0,intxn_energies,total_intxn_energy);
            }
        } else if (strcasecmp("ANALYZE",token)==0) {
            if (!initialized) finish_initialization();
            token=strtok(NULL,delim);
            strncpy(word,token,sizeof(word));
            sscanf(word,"%s %s %s",word,fname,fname2);
            energy_analysis(word,fname,fname2);
        } else if (strcasecmp("SETLIGANDCOM",token)==0) {
            //possibilities SETLIGANDCOM real real real or SETLIGANDCOM
            token+=strlen(token)+1;
            strncpy(word,token,sizeof(word)); //entire rest of line
            if (!initialized) finish_initialization();
            if (strstr(word,"AAREGION")!=NULL) {
                printf("Aligning center of mass of ligand with center of mass of AA region.\n");
                align_ligand_with_aa_region(initcoords);
            } else {
                sscanf(word,"%lg %lg %lg",&desired_ligand_com[0],&desired_ligand_com[1],&desired_ligand_com[2]);
                set_ligand_com(desired_ligand_com,initcoords);
            }
        } else if (strcasecmp("DOCKPREP",token)==0) {
            token+=strlen(token)+1;
            strncpy(word,token,sizeof(word)); //entire rest of line
            sscanf(word,"%lg %lg %d",&trans_size, &rot_size, &nsearch);
            rot_size*=DEG_TO_RAD;
            if (!initialized) finish_initialization();
            prepare_docking(trans_size,rot_size,nsearch,initcoords);
        } else if ((strcasecmp("RUN",token)==0) || ((strcasecmp("MC",token)==0))) {
            token=strtok(NULL,delim);
            strncpy(word,token,sizeof(word));
            sscanf(word,"%d",&nmcstep);
            if (!initialized) finish_initialization();
            mcloop();
        } else {
            //check for Go parameters -- see go_model.cpp
            /*strncpy(word,token,sizeof(word)); //first word
            token+=strlen(token)+1;
            if (read_go_parameter(word,token,&go_params)) {
            } else {*/
                printf("Unrecognized command.\n");
                die();
           //}
        } //end of large if..else
    } //while(true)
    fclose(input);
}

void simulation::read_go_model_info(FILE * input) // go_model_params * params)
{
    char command[255],command2[255],  model_fname[255];
    char * token;
    char * rest_of_line;
    const char * delim = " \t\n";
    double * model_coords;
    subset valid_coords;
    memset(model_fname,'\0',sizeof(model_fname));
    for (;;) {
        fgets(command,sizeof(command),input);
        strncpy(command2,command,sizeof(command2));//This prevents parsing of the original command.
        token=strtok(command2,delim);
        if (token==NULL) continue; //blank line
        if (strcasecmp("END",token)==0) break;
        rest_of_line=token+strlen(token)+1; //point after the token
        if (strcasecmp("HARDCORE",token)==0) {
            sscanf(rest_of_line,"%lg",&go_params.hardcore);
        } else if (strcasecmp("NATIVE_CUTOFF",token)==0) {
            sscanf(rest_of_line,"%lg",&go_params.cutoff);
        } else if (strcasecmp("EXPONENTS",token)==0) {
            sscanf(rest_of_line,"%d %d",&go_params.m, &go_params.n);
        } else if (strcasecmp("WELLDEPTH",token)==0) {
            sscanf(rest_of_line,"%lg",&go_params.native_energy);
            go_params.nonnative_energy=go_params.native_energy;
        } else if (strcasecmp("NATIVE_STATE",token)==0) {
	    token=strtok(NULL,delim);
            strncpy(model_fname,token,sizeof(model_fname));
        } else {
            printf("Unrecognized Go model command.\n");
            die();
        }
    }
    finish_go_params(&go_params);
    print_go_params(go_params);
    go_model = new go_model_info();
    model_coords=(double *) checkalloc(3*top->natom,sizeof(double));
    printf("Creating native distances for Go model from PDB file %s.\n",model_fname);
    top->read_pdb_file(model_fname,model_coords,valid_coords);
    go_model->create_contact_map(top->nres,top->resinfo,top->natom,model_coords,&go_params,top->aaregion_res);
    free(model_coords);

}

void simulation::create_lists(void)
{
    top->aaregion_res.print("All atom region residues ");
    if ((top==NULL) || (ffield==NULL)) {
        printf("You need to load the definitions and parameter files first.\n");
        die();
    }
    if (nres<=0) {
        printf("You need to specify a sequence.\n");
        die();
    }
    /*if (!aaregion_specified) {
        printf("All atom region not explicitly specified.  By default assuming that all atoms are to be in the all-atom region.\n");
        aaregion_res.init(nres);
        for (i=0; i<nres; i++) aaregion_res+=i;
    }*/
    //Do we do the below stuff here or in finish_initialization?
    top->add_segment(' ',sequence,ffield);
    if (top->ligand_res>=0) top->add_segment(' ',ligand_resname,ffield);
    top->create_angle_dihedral_lists(false);
    top->create_non_tab_list();
    ffield->find_parameters(top->natom,top->atoms);
    top->create_improper_dihedral_lists(false,ffield);
}

void simulation::finish_initialization(void)
{
    int i,iatom;
    if (!aaregion_specified) {
        printf("All atom region not explicitly specified.  By default assuming that all atoms are to be in the all-atom region.\n");
        top->aaregion_res.init(nres);
        for (i=0; i<nres; i++) top->aaregion_res+=i;
        create_lists();
    }
    top->aaregion_res.print("All atom region residues ");
    //this really should be done in topology's scope
    top->ligand.init(top->natom);
    for (i=0; i<top->natom; i++) if (top->atoms[i].resNum==top->ligand_res) top->ligand+=i;
#ifdef DEBUG_NON_TABULATED
    top->print_detailed_info(aaregion_res);
#else
    top->print_summary_info();
#endif
    printf("Will read PDB file %s.\n",structfname);
    initcoords=(double *) checkalloc(3*top->natom,sizeof(double));
    top->read_pdb_file(structfname,initcoords,valid_coords);
    for (iatom=0; iatom<top->natom; iatom++) if (!valid_coords[iatom]) {
        printf("Failed to find coordinates for atom %s %d %s\n",top->atoms[iatom].resName,top->atoms[iatom].resNum+1,top->atoms[iatom].name);
        die();
    }
    top->check_solvation_params(&solvation_params);
    //ffield->build_coords(initcoords,top->natom,top->atoms,valid_coords);
    //we need to have set up go parameters by now
    /*finish_go_params(&go_params);
    go_model = new go_model_info();
    go_model->create_contact_map(top->nres,top->resinfo,top->natom,initcoords,&go_params,top->aaregion_res);*/

    //Create nonbond list object.  Set listcutoff=cutoff or less to disable the nb list.
    /*use_nb_list=(listcutoff>cutoff);
    if (use_nb_list) frag_nblist=new fragment_nblist(top->nfrag,listcutoff);*/

    cutoff2=cutoff*cutoff;
    //use_nb_list=(listcutoff>cutoff);
    use_nb_list=true;
    if (listcutoff<cutoff) listcutoff=cutoff;
    if (pbc) {
        printf("PBC is on. Box size =      %.2f A\n",boxsize);
    } else {
        printf("PBC is off.\n");
    }
    if (rdie) {
        printf("Distance dependent dielectric will be used.\n");
    }
    //if (use_nb_list) {
        printf("Nonbond list will be used.\n");
        printf("Spherical/list cutoff         %.2f %.2f\n",cutoff,listcutoff);
    /*} else {
        printf("Nonbond list will not be used.\n");
        printf("Spherical cutoff              %.2f\n",cutoff);
    }*/
    printf("Dielectric constant:          %.2f\n",eps);
    //print_go_params(go_params);
    if (seed==0) {
       	   	//seed=time(NULL);
        seed=read_random_seed();
        printf("Initialized seed based on /dev/urandom.  Seed used is %ld\n",seed);
    } else {
        printf("Seed specified in file.  Seed used is %ld\n",seed);
    }
    init_genrand(seed);
    top->generate_backbone_moves(&backbone_moves);
    top->generate_sidechain_moves(&sidechain_moves,&ligand_bond_rotation_moves);
    top->generate_backrub_moves(&backrub_moves);
    if ((prob[MOVE_SIDECHAIN]>0) && (sidechain_moves.size()<=0)) {
        //this usually
        printf("You cannot have sidechain moves with an empty AA region.\n");
        die();
    }
    printf("Number of monte carlo steps: %ld\n",nmcstep);
    printf("Save and print frequencies:  %ld %ld\n",nsave,nprint);
    //need these arrays for energy calculations
#ifdef SEDDD
    old_frac_volumes=(double *) checkalloc(top->natom,sizeof(double));
    new_frac_volumes=(double *) checkalloc(top->natom,sizeof(double));
#endif
    initialized=true;
}

void simulation::get_ligand_com(double * coords, double * com_ligand)
{
    int iatom,k;
    double mass_ligand;
    for (k=0; k<3; k++) com_ligand[k]=0.0;
    mass_ligand=0.0;
    for (iatom=0; iatom<top->natom; iatom++) if (top->ligand[iatom]) {
        mass_ligand+=top->atoms[iatom].mass;
        for (k=0; k<3; k++) com_ligand[k]+=top->atoms[iatom].mass*coords[3*iatom+k];
    }
    for (k=0; k<3; k++) com_ligand[k]/=mass_ligand;
}

//position the ligand's center of mass
void simulation::set_ligand_com(double * desired_com, double * coords)
{
    double com_ligand[3], mass_ligand, disp[3];
    int iatom, k;
    //Move the ligand COM to superimpose over the COM of the all-atom region.
    get_ligand_com(coords,com_ligand);
    for (k=0; k<3; k++) disp[k]=desired_com[k]-com_ligand[k];
    for (iatom=0; iatom<top->natom; iatom++) if (top->ligand[iatom]) {
        for (k=0; k<3; k++) coords[3*iatom+k]+=disp[k];
    }
    printf("Set ligand COM to coordinates %.3f %.3f %.3f\n",desired_com[0],desired_com[1],desired_com[2]);
}


void simulation::align_ligand_with_aa_region(double * coords)
{
    double com_aa_region[3], mass_aa;
    int iatom, k;
    //Move the ligand COM to superimpose over the COM of the all-atom region.
    for (k=0; k<3; k++) com_aa_region[k]=0.0;
    mass_aa=0.0;
    for (iatom=0; iatom<top->natom; iatom++) if ((top->atoms[iatom].is_in_aa_region) && (!top->ligand[iatom])) {
        mass_aa+=top->atoms[iatom].mass;
        for (k=0; k<3; k++) com_aa_region[k]+=top->atoms[iatom].mass*coords[3*iatom+k];
    }
    for (k=0; k<3; k++) com_aa_region[k]/=mass_aa;
    set_ligand_com(com_aa_region,coords);
}

void simulation::prepare_docking(double trans_size, double rot_size, int nsearch, double * coords)
{
    double com_aa_region[3],com_ligand[3],disp[3],mass_aa,mass_ligand,x[3],q[4],rmsd,junk;
    int iatom,k,imove;
    double angle;
    double * private_coords;
    double * weight;
    double * private_coords2;
    double * best_coords;
    double etot, best_etot, energies[EN_TERMS], best_energies[EN_TERMS];
    int isearch,ibest;
#ifdef MONITORPREDOCK
    FILE * predock_poses; //pdb file containing predock poses
#endif
    private_coords = (double *) checkalloc(3*top->natom,sizeof(double));
    weight = (double *) checkalloc(top->natom,sizeof(double));
    printf("Preparing for docking -- random displacement %.2f A, random rotation %.2f deg\n",
        trans_size, rot_size*RAD_TO_DEG);
    for (iatom=0; iatom<top->natom; iatom++) {
        //weight is to be used for RMSD analysis of heavy-atom RMSD in ligand.
        if (top->ligand[iatom] && top->atoms[iatom].atomicNum>1) weight[iatom]=1.0; else weight[iatom]=0.0;
        for (k=0; k<3; k++) private_coords[3*iatom+k]=coords[3*iatom+k];
    }
    printf("Randomizing bonds in ligand.\n");
    for (imove=0; imove<ligand_bond_rotation_moves.size(); imove++) 
	if (!ligand_bond_rotation_moves[imove].is_stiff) {
            //in theory there should be no bonds connecting the ligand to anything, but for safety we check both
            //if (top->ligand[ligand_bond_rotation_moves[imove].iaxis] && top->ligand[ligand_bond_rotation_moves[imove].jaxis]) {
                angle=(2.0*genrand_real3()-1.0)*M_PI;
                rotate_atoms_by_axis(&ligand_bond_rotation_moves[imove],angle,private_coords);
            //}
        }
    //Begin search for a low-energy ligand conformation
    printf("Generating %d random ligand conformations to search for a low-energy initial pose.\n",nsearch);
    private_coords2 = (double *) checkalloc(3*top->natom,sizeof(double));
    best_coords = (double *) checkalloc(3*top->natom,sizeof(double));

    best_etot=DUMMY_ENERGY; //very high
    ibest=-1;
#ifdef MONITORPREDOCK
    predock_poses=fopen("predock_poses.pdb","w");
#endif
    for (isearch=1; isearch<=nsearch; isearch++) {
        //restore coordinates from original
        for (iatom=0; iatom<top->natom; iatom++) for (k=0; k<3; k++) private_coords2[3*iatom+k]=private_coords[3*iatom+k];
        //give a random displacement and orientation
        do_ligand_trans(trans_size,0.0,trans_size,&junk,private_coords2);
        do_ligand_rot(rot_size,0.0,rot_size,&junk,private_coords2);
        //rotate by random extents around rotatable bonds -- do not include  "stiff" bonds
        for (imove=0; imove<ligand_bond_rotation_moves.size(); imove++)
            if (!ligand_bond_rotation_moves[imove].is_stiff) {
            //in theory there should be no bonds connecting the ligand to anything, but for safety we check both
            //if (top->ligand[ligand_bond_rotation_moves[imove].iaxis] && top->ligand[ligand_bond_rotation_moves[imove].jaxis]) {
                angle=(2.0*genrand_real3()-1.0)*M_PI;
                rotate_atoms_by_axis(&ligand_bond_rotation_moves[imove],angle,private_coords2);
            //}
        }
        //check teh energy
        top->create_pair_list(pbc,halfboxsize,boxsize,listcutoff,&old_pair_list,&old_solv_list,private_coords2);
#ifdef SEDDD
        top->calculate_solvation_volumes(&solvation_params,cutoff2,&old_solv_list,private_coords2,old_frac_volumes,ffield);
        total_energy(private_coords2,&old_pair_list,old_frac_volumes,energies,&etot);
#else
        total_energy(private_coords2,&old_pair_list,energies,&etot);
#endif
        if (etot<best_etot) {
            //if lower energy than best so far, save it
            best_etot=etot;
            ibest=isearch;
            for (k=0; k<EN_TERMS; k++) best_energies[k]=energies[k];
            for (iatom=0; iatom<top->natom; iatom++) for (k=0; k<3; k++) best_coords[3*iatom+k]=private_coords2[3*iatom+k];
        }
        if (isearch%nprint==0) print_energies(stdout,(isearch==nprint),"Search:",isearch,best_energies,best_etot);
#ifdef MONITORPREDOCK
        fprintf(predock_poses,"REMARK random poses %d, total energy %g\n",isearch,etot);
        fprintf(predock_poses,"MODEL     %4d\n",isearch);
        top->write_pdb_stream(predock_poses,private_coords2,old_frac_volumes);
        fprintf(predock_poses,"ENDMDL\n");
        fflush(predock_poses);
#endif
    }
#ifdef MONITORPREDOCK
    fclose(predock_poses);
    printf("Best pose is number %d\n",ibest);
#endif
    rmsd_fit(top->natom,weight,best_coords,coords,&x[0],&q[0],&rmsd);
    printf("Initial ligand RMSD: %.3f\n",rmsd);
    for (iatom=0; iatom<top->natom; iatom++) {
        for (k=0; k<3; k++) coords[3*iatom+k]=best_coords[3*iatom+k];
    }
    free(private_coords);
    free(private_coords2);
    free(best_coords);
    free(weight);
}
//Supervises the loading of tables. fmt is a format string with the names of fragments.
//void simulation::load_tables(char * fmt)
/*void topology::load_tables(const char * fmt, const char * fragfmt, table * * tables)
{
    int num_tables, count, frags_in_use, ifrag,jfrag,itype,jtype,i;
    double total_size;
    FILE * f;
    bool * need_table;
    char fname[255];
    need_table=(bool *) malloc(nfragtypes*nfragtypes*sizeof(bool));
    num_tables=0;
    for (i=0; i<nfragtypes*nfragtypes; i++) need_table[i]=false;
    for (ifrag=0; ifrag<nfrag; ifrag++)
        for (jfrag=0; jfrag<nfrag; jfrag++)
            if ((ifrag!=jfrag) && (!closefragments[ifrag*nfrag+jfrag])) {
                itype=frags[ifrag].type;
                jtype=frags[jfrag].type;
                need_table[itype*nfragtypes+jtype]=true;
            }
    for (itype=0; itype<nfragtypes; itype++)
        for (jtype=itype; jtype<nfragtypes; jtype++)
            if (need_table[itype*nfragtypes+jtype]) num_tables++;
    //total_size=0.0;
    frags_in_use=0;
    for (ifrag=0; ifrag<nfragtypes; ifrag++) if (fragtypes[ifrag]->n_used>0) frags_in_use++;
    //num_tables=frags_in_use*(frags_in_use+1)/2;
    printf("Fragment types in use: %d\n",frags_in_use);
    printf("Need to load a total of %d tables.\n",num_tables);
    count=0;
    for (itype=0; itype<nfragtypes; itype++)
        for (jtype=itype; jtype<nfragtypes; jtype++)
            if (need_table[itype*nfragtypes+jtype]) {
                //Maybe this code for finding the file name should be moved to table's constructor.
                count++;
                printf("\n");
                printf("--------------------------------------------------------------------------------\n");
                printf("Loading interaction table file %d of %d.\n",count,num_tables);
                //ifrag is the lesser fragment type.  Load the table! (Constructor will try both possible names.
                tables[itype*nfragtypes+jtype]=new table(fmt,fragfmt,fragtypes[itype]->fragname,fragtypes[jtype]->fragname,this);
                //tables[jtype*nfragtypes+itype]=tables[itype*nfragtypes+jtype];
                //total_size+=tables[itype*nfragtypes+jtype]->getsize();
		fflush(stdout);
            }
    free(need_table);
    printf("Total interaction tables loaded: %d\n",count);

    //printf("Total standard table size:     %.2f MB\n",total_size);
}

void topology::load_covalent_tables(const char * covtablefmt, covalent_table * * cov_tables)
{
    int num_cov_tables, count, frags_in_use, seg, ifrag,jfrag,itype,jtype,i;
    double total_size;
    FILE * f;
    bool * need_table;
    char fname[255];
    need_table=(bool *) malloc(nfragtypes*nfragtypes*sizeof(bool));
    num_cov_tables=0;
    for (i=0; i<nfragtypes*nfragtypes; i++) need_table[i]=false;
    for (seg=0; seg<nseg; seg++) {
        ifrag=first_main_chain_frag[seg];
        jfrag=frags[ifrag].main_chain_next;
        while (jfrag>0) {
            itype=frags[ifrag].type;
            jtype=frags[jfrag].type;
            //we rely on short-circuit evaluaton here
            if (use_covalent_table(itype,jtype)) need_table[itype*nfragtypes+jtype]=true;
            ifrag=jfrag;
            jfrag=frags[jfrag].main_chain_next;
        }
    }
    for (itype=0; itype<nfragtypes; itype++)
        for (jtype=0; jtype<nfragtypes; jtype++)
            if (need_table[itype*nfragtypes+jtype]) num_cov_tables++;
    //total_size=0.0;
    printf("Need to load a total of %d covalent tables.\n",num_cov_tables);
    count=0;
    for (itype=0; itype<nfragtypes; itype++)
        for (jtype=0; jtype<nfragtypes; jtype++)
            if (need_table[itype*nfragtypes+jtype]) {
                count++;
                printf("\n");
                printf("--------------------------------------------------------------------------------\n");
                printf("Loading covalent table file %d of %d.\n",count,num_cov_tables);
                snprintf(fname,sizeof(fname),covtablefmt,fragtypes[itype]->fragname,fragtypes[jtype]->fragname);
                strlower(fname,sizeof(fname));
                trim_string(fname);
                //ifrag is the lesser fragment type.  Load the table! (Constructor will try both possible names.
                cov_tables[itype*nfragtypes+jtype]=new covalent_table(fname,false);
                //tables[jtype*nfragtypes+itype]=tables[itype*nfragtypes+jtype];
                //total_size+=cov_tables[itype*nfragtypes+jtype]->getsize();
		fflush(stdout);
            }
    free(need_table);
    printf("Total covalent tables loaded: %d\n",count);
    //printf("Total covalent table size:     %.2f MB\n",total_size);
}*/
//this does not include the go term!
void simulation::ligand_energies(double * coords, double * total_internal_energy, double * internal_energies, double * total_intxn_energy, double * intxn_energies)
{
    int i;
    double * frac_volumes;
    std::vector<atom_nb_entry> pair_list, solv_list;
    for (i=0; i<EN_TERMS; i++) {
        internal_energies[i]=0.0;
        intxn_energies[i]=0.0;
    }
    top->create_pair_list(pbc,halfboxsize,boxsize,listcutoff,&pair_list,&solv_list,coords);
#ifdef SEDDD
    frac_volumes = (double *) checkalloc(top->natom,sizeof(double));
    top->calculate_solvation_volumes(&solvation_params,cutoff2,&solv_list,coords,frac_volumes,ffield);
    ffield->subset_energy(&solvation_params,cutoff2,top->natom,top->atoms,top->ligand,pair_list.size(),&pair_list[0],coords,frac_volumes,internal_energies,intxn_energies);
    free(frac_volumes);
#else
    ffield->subset_energy(eps,rdie,cutoff2,top->natom,top->atoms,top->ligand,pair_list.size(),&pair_list[0],coords,internal_energies,intxn_energies);
#endif
    *total_internal_energy=0.0;
    *total_intxn_energy=0.0;
    for (i=0; i<EN_TERMS; i++) {
        *total_internal_energy+=internal_energies[i];
        *total_intxn_energy+=intxn_energies[i];
    }
}

void simulation::energy_analysis(char * type, char * fname,  char * enfname)
{
    long int frame,istep;
    int ifrag,i;
    FILE * input;
    FILE * enoutput;
    subset valid_coords;
    double intxn_energies[EN_TERMS],internal_energies[EN_TERMS],total_internal_energy,total_intxn_energy,energies[EN_TERMS],etot;
    top->chaincodes[0]='P'; //hack to make compatible with pdb files derived from catdcd
    if ((enfname==NULL) || (strlen(enfname)==0)) enoutput=stdout; else enoutput=fopen(enfname,"w");
    if (enoutput==NULL) {
        printf("Could not open energies file %s\n",enfname);
        die();
    }
    /*if (strcasecmp(fname,"now")==0) {
        //if (strcasecmp(type,"total",)==0) {
        ligand_energies(oldcoords,&total_internal_energy,internal_energies,&total_intxn_energy,intxn_energies);
        print_energies(enoutput,true,"Internal:",0,internal_energies,total_internal_energy);
        print_energies(enoutput,true,"Intxn:",0,intxn_energies,total_intxn_energy);
    } else {*/
        input=fopen(fname,"r");
        if (input==NULL) {
            printf("Could not open trajectory file %s\n",fname);
            die();
        }
        frame=1;
        while (!feof(input)) {
            printf("Reading frame %ld\n",frame);
            top->read_pdb_stream(input,oldcoords,valid_coords);
            if (feof(input)) break;
            if (strcasecmp(type,"TOTAL")==0) {
                top->create_pair_list(pbc,halfboxsize,boxsize,listcutoff,&old_pair_list,&old_solv_list,initcoords);
#ifdef SEDDD
                top->calculate_solvation_volumes(&solvation_params,cutoff2,&old_solv_list,initcoords,old_frac_volumes,ffield);
                total_energy(initcoords,&old_pair_list,old_frac_volumes,energies,&etot);
#else
                total_energy(initcoords,&old_pair_list,energies,&etot);
#endif
                print_energies(enoutput,true,"Energy:",0,energies,etot);
            } else if (strcasecmp(type,"LIGAND")==0) {
                ligand_energies(oldcoords,&total_internal_energy,internal_energies,&total_intxn_energy,intxn_energies);
                print_energies(enoutput,true,"Internal:",0,internal_energies,total_internal_energy);
                print_energies(enoutput,true,"Intxn:",0,intxn_energies,total_intxn_energy);
            }
            frame++;
        }
    //}
    fclose(input);
    fclose(enoutput);
}

