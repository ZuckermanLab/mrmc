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
    seed = 0;
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
    free(initcoords);
    /*free(oldcoords);
    free(oldcenter);
    free(oldorient);
    free(newcoords);
    free(newcenter);
    free(neworient);*/
    //if (frag_nblist!=NULL) delete frag_nblist;
    delete ffield;
    delete go_model;
    delete top;
}

void simulation::process_commands(char * infname)
{
    char command[255],command2[255],fname[255],fmt[255],word[255],word2[255];
    char * token;
    const char * delim = " \t\n";
    char chain;
    FILE * input;
    FILE * output;
    double p,size,ptot,mctemp, trans_size, rot_size, bond_rot_size;
    int i,move;

    double energies[EN_TERMS],etot;
    mctemp=300.0;
    initialized = false;
    input=fopen(infname,"r");
    if (input==NULL) {
        printf("Could not open command file %s.\n",infname);
        die();
    }
    while (true) {
        fgets(command,sizeof(command),input);
        if (feof(input)) break; //end of file
        //for (i=0; i<strlen(command); i++) command[i]=toupper(command[i]);
        printf("Processing command: %s\n",command);
        //for (i=0; i<strlen(command); i++) command[i]=toupper(command[i]);
        strncpy(command2,command,sizeof(command2));
        token=strtok(command2,delim);
        if (token==NULL) continue; //blank line
        else if (*token=='#') continue; //comment
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
                top->write_pdb_file(fname,oldcoords);
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
            } else {
                printf("Unrecognized format.\n");
                die();
            }
        } else if (strcasecmp("INSERT",token)==0) {
            token=strtok(NULL,delim);
            strncpy(word,token,sizeof(word));
            if (strcasecmp("SEQUENCE",word)==0) {
                sequence=read_multiline(input);
                nres=count_words(sequence);
                if (nres==0) {
                    printf("You need to specify a sequence first.\n");
                    die();
                }
            } else if (strcasecmp("LIGAND",word)==0) {
                token=strtok(NULL,delim);
                strncpy(ligand_resname,token,sizeof(ligand_resname));
                top->ligand_res=nres;
                nres++;
            } else {
                printf("Unrecognized command.\n");
                die();
            }
        } else if (strcasecmp("AAREGION",token)==0) {
            aaregion_res.init(nres);
            token=strtok(NULL,delim);
            strncpy(word,token,sizeof(word));
            aaregion_res.parse_int_list(word);
            if (top->ligand_res>=0) aaregion_res+=top->ligand_res; //force inclusion of ligand in AA region for docking
            aaregion_specified = true;
            create_lists();
        } else if (strcasecmp("MOVES",token)==0) {
            for (i=1; i<=NUM_MOVES; i++) prob[i]=0.0;
            for(;;) {
                fgets(command,sizeof(command),input);
                sscanf(command,"%s %lg %lg\n",word,&p,&size);
                if (strcasecmp("END",word)==0) break;
                move=-1;
                for (i=1; i<=NUM_MOVES; i++) if (strcasecmp(mc_move_names[i],word)==0) move=i;
                if (move<0) break;
                prob[move]=p;
                movesize[move]=size*DEG_TO_RAD;
            }
            ptot=0.0;
            for(i=1;i<=NUM_MOVES;i++)ptot+=prob[i];
            for(i=1;i<=NUM_MOVES;i++)prob[i]/=ptot;
            cumprob[1]=prob[1];
            for(i=2;i<=NUM_MOVES;i++)cumprob[i]=cumprob[i-1]+prob[i];
            for(i=1;i<=NUM_MOVES;i++)printf("%.10s moves:  Maximum size %.2f degrees  Fraction %.2f%%\n",
                mc_move_names[i],movesize[i]*RAD_TO_DEG,prob[i]*100.0);
        } else if (strcasecmp("BOXSIZE",token)==0) {
            pbc=true;
            token=strtok(NULL,delim);
            strncpy(word,token,sizeof(word));
            sscanf(word,"%lg",&boxsize);
            halfboxsize=0.5*boxsize;
        } else if (strncasecmp("CUTOFF",token,6)==0) {
            token=strtok(NULL,delim);
            strncpy(word,token,sizeof(word));
            sscanf(word,"%lg",&cutoff);
            cutoff2=cutoff*cutoff;
        } else if (strcasecmp("EPS",token)==0) {
            token=strtok(NULL,delim);
            strncpy(word,token,sizeof(word));
            sscanf(word,"%lg",&eps);
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
            strncpy(xyzfname,token,sizeof(xyzfname));
        } else if (strcasecmp("ENERGY",token)==0) {
            finish_initialization();
            total_energy(initcoords,energies,&etot);
            print_energies(stdout,true,"Energy:",0,energies,etot);
        } else if (strcasecmp("DOCKPREP",token)==0) {
	    token+=strlen(token)+1; 
            strncpy(word,token,sizeof(word)); //entire rest of line		
            sscanf(word,"%lg %lg %lg",&trans_size, &rot_size, &bond_rot_size);
            rot_size*=DEG_TO_RAD;
            bond_rot_size*=DEG_TO_RAD;
            if (!initialized) finish_initialization();
            prepare_docking(trans_size,rot_size,bond_rot_size,initcoords);
        } else if ((strcasecmp("RUN",token)==0) || ((strcasecmp("MC",token)==0))) {
            token=strtok(NULL,delim);
            strncpy(word,token,sizeof(word));
            sscanf(word,"%d",&nmcstep);
            if (!initialized) finish_initialization();
            mcloop();
        } else {
            //check for Go parameters -- see go_model.cpp
            strncpy(word,token,sizeof(word)); //first word
            token+=strlen(token)+1;
            if (read_go_parameter(word,token,&go_params)) {
            } else {
                printf("Unrecognized command.\n");
                die();
            }
        } //end of large if..else
    } //while(true)
    fclose(input);
}

void simulation::create_lists(void)
{
    aaregion_res.print("All atom region residues ");
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
    top->add_segment(' ',sequence,aaregion_res);
    if (top->ligand_res>=0) top->add_segment(' ',ligand_resname,aaregion_res);
    top->create_angle_dihedral_lists(false);
    top->create_non_tab_list(false,&non_tab_list);
    ffield->find_parameters(top->natom,top->atoms);
    top->create_improper_dihedral_lists(false,ffield);
}

void simulation::finish_initialization(void)
{
    int i;
    if (!aaregion_specified) {
        printf("All atom region not explicitly specified.  By default assuming that all atoms are to be in the all-atom region.\n");
        aaregion_res.init(nres);
        for (i=0; i<nres; i++) aaregion_res+=i;
        create_lists();
    }

    aaregion_res.print("All atom region residues ");
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
    top->read_pdb_file(structfname,initcoords);
    //we need to have set up go parameters by now
    finish_go_params(&go_params);
    go_model = new go_model_info();
    go_model->create_contact_map(top->nres,top->resinfo,top->natom,initcoords,&go_params,aaregion_res);

    //Create nonbond list object.  Set listcutoff=cutoff or less to disable the nb list.
    /*use_nb_list=(listcutoff>cutoff);
    if (use_nb_list) frag_nblist=new fragment_nblist(top->nfrag,listcutoff);*/

    cutoff2=cutoff*cutoff;
    if (pbc) {
        printf("PBC is on. Box size =      %.2f A\n",boxsize);
    } else {
        printf("PBC is off.\n");
    }
    if (rdie) {
        printf("Distance dependent dielectric will be used.\n");
    }
    /*if (use_nb_list) {
        printf("Nonbond list will be used.\n");
        printf("Spherical/list cutoff         %.2f %.2f\n",cutoff,listcutoff);
    } else {*/
        printf("Nonbond list will not be used.\n");
        printf("Spherical cutoff              %.2f\n",cutoff);
    //}
    printf("Dielectric constant:          %.2f\n",eps);
    print_go_params(go_params);
    if (seed==0) {
       	   	seed=time(NULL);
#ifdef UNIX
            seed^=getpid(); //To ensure uniqueness even if many are started at same time.
#endif
           	printf("Initialized seed based on time.  Seed used is %ld\n",seed);
    } else {
           	printf("Seed specified in file.  Seed used is %ld\n",seed);
    }
    init_genrand(seed);
    top->generate_backbone_moves(&backbone_moves);
    top->generate_sidechain_moves(&sidechain_moves);
    top->generate_backrub_moves(&backrub_moves);
    if ((prob[MOVE_SIDECHAIN]>0) && (sidechain_moves.size()<=0)) {
        //this usually
        printf("You cannot have sidechain moves with an empty AA region.\n");
        die();
    }
    printf("Number of monte carlo steps: %ld\n",nmcstep);
    printf("Save and print frequencies:  %ld %ld\n",nsave,nprint);


    initialized=true;
}


void simulation::prepare_docking(double trans_size, double rot_size, double bond_rot_size, double * coords)
{
    double com_aa_region[3],com_ligand[3],disp[3],mass_aa,mass_ligand,x[3],q[4],rmsd;
    int iatom,k,imove;
    double angle;
    double * private_coords;
    double * weight;
    private_coords = (double *) checkalloc(3*top->natom,sizeof(double));
    weight = (double *) checkalloc(top->natom,sizeof(double));
    printf("Preparing for docking -- random displacement %.2f A, random rotation %.2f deg, random bond rotation %.2f deg\n",
        trans_size, rot_size*RAD_TO_DEG, bond_rot_size*RAD_TO_DEG);
    for (iatom=0; iatom<top->natom; iatom++) {
        //weight is to be used for RMSD analysis of heavy-atom RMSD in ligand.
        if (top->ligand[iatom] && top->atoms[iatom].atomicNum>1) weight[iatom]=1.0; else weight[iatom]=0.0;
        for (k=0; k<3; k++) private_coords[3*iatom+k]=coords[3*iatom+k];
    }
    //Move the ligand COM to superimpose over the COM of the all-atom region.
    for (k=0; k<3; k++) com_aa_region[k]=0.0;
    for (k=0; k<3; k++) com_ligand[k]=0.0;
    for (iatom=0; iatom<top->natom; iatom++) {
        if ((top->atoms[iatom].is_in_aa_region) && (!top->ligand[iatom])) {
            mass_aa+=top->atoms[iatom].mass;
            for (k=0; k<3; k++) com_aa_region[k]+=top->atoms[iatom].mass*private_coords[3*iatom+k];
        }
        if (top->ligand[iatom]) {
            mass_ligand+=top->atoms[iatom].mass;
            for (k=0; k<3; k++) com_ligand[k]+=top->atoms[iatom].mass*private_coords[3*iatom+k];
        }
    }
    for (k=0; k<3; k++) {
        com_aa_region[k]/=mass_aa;
        com_ligand[k]/=mass_ligand;
        disp[k]=com_aa_region[k]-com_ligand[k];
    }
    for (iatom=0; iatom<top->natom; iatom++) if (top->ligand[iatom]) {
            for (k=0; k<3; k++) private_coords[3*iatom+k]+=disp[k];
    }
    //give a random displacement and orientation
    do_ligand_trans(trans_size,private_coords);
    do_ligand_rot(rot_size,private_coords);
    //rotate by random extents around rotatable bonds
    for (imove=0; imove<sidechain_moves.size(); imove++)
        //in theory there should be no bonds connecting the ligand to anything, but for safety we check both
        if (top->ligand[sidechain_moves[imove].iaxis] && top->ligand[sidechain_moves[imove].jaxis]) {
            angle=(2.0*genrand_real3()-1.0)*bond_rot_size;
            rotate_atoms_by_axis(&sidechain_moves[imove],angle,private_coords);
        }
    rmsd_fit(top->natom,weight,private_coords,coords,&x[0],&q[0],&rmsd);
    printf("Initial ligand RMSD: %.3f\n",rmsd);
    for (iatom=0; iatom<top->natom; iatom++) {
        for (k=0; k<3; k++) coords[3*iatom+k]=private_coords[3*iatom+k];
    }
    free(private_coords);
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
                strlower(fname);
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

/*void simulation::do_energies(char * type, char * fname,  char * enfname)
{
    long int frame,istep;
    int ifrag;
    FILE * input;
    FILE * enoutput;
    double * rmsds;
    double energies[EN_TERMS],etot;
    top->chaincodes[0]='P'; //hack to make compatible with pdb files derived from catdcd
    rmsds=(double *) checkalloc(top->nfrag,sizeof(double));
    initcoords=(double *) checkalloc(3*top->natom,sizeof(double));
    en_by_table=(double *) checkalloc(top->nfragtypes*top->nfragtypes,sizeof(double));
    if ((enfname==NULL) || (strlen(enfname)==0)) enoutput=stdout; else enoutput=fopen(enfname,"w");
    if (enoutput==NULL) {
        printf("Could not open energies file %s\n",enfname);
        die();
    }
    if (strcasecmp(type,"rest")==0) {  //Calculate energies from restart file. Needed for replica exchange.
        read_restart(fname);
        for (ifrag=0; ifrag<top->nfrag; ifrag++) top->update_coords(ifrag,oldcenter,oldorient,oldcoords);
        total_energy(oldcenter,oldorient,oldcoords,energies,&etot);{
    FILE * f;
    unsigned long seed;
    double listcutoff;
    double ptot,mctemp,p,size;
    double total_table_size;
    char buffer[255],word[255],tablefmt[255],covtablefmt[255],fragfmt[255],structf3name[255],psffname[255];
    char * seq;
    char * pch;
    bool start,restart;
    int itype,ifrag,i,jty
        print_energies(stdout,false,"Energy: ",0,energies,etot);
        return;
    }
    input=fopen(fname,"r");
    if (input==NULL) {
        printf("Could not open trajectory file %s\n",fname);
        die();
    }
    frame=1;
    while (!feof(input)) {
        printf("Reading frame %ld\n",frame);
        if (strcasecmp(type,"pdb")==0) {
            //No need for diagnostic messages.  Later we can implement statistics on rmsds.
            top->read_pdb_stream(input,initcoords);
            top->fit_all_fragments(initcoords,oldcenter,oldorient,oldcoords,rmsds);
        } else if (strcasecmp(type,"dat")==0) {
            read_frame_quat(input,&istep,oldcenter,oldorient);
            printf("Step number %ld\n",istep);
            for (ifrag=0; ifrag<top->nfrag; ifrag++) top->update_coords(ifrag,oldcenter,oldorient,oldcoords);
        }
        if (feof(input)) break;
        total_energy(oldcenter,oldorient,oldcoords,energies,&etot);
        print_energies(enoutput,(frame==1),"Energy: ",frame,energies,etot);
        if (use_std_tables) print_energies_by_table();
        frame++;
    }
    fclose(input);
    fclose(enoutput);
    free(rmsds);
}*/

