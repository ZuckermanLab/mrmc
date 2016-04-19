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
simulation::simulation(const char * command, const char * fname, int _mynod, int _numnod)
#else
simulation::simulation(const char * command, const char * fname)
#endif
{
    FILE * f;
    unsigned long seed;
    double listcutoff;
    double ptot,mctemp,p,size;
    double total_table_size;
    char buffer[255],word[255],tablefmt[255],covtablefmt[255],fragfmt[255],structf3name[255],psffname[255];
    char * seq;
    char * pch;
    bool start,restart,do_mc,do_energy;
    int itype,ifrag,i,jtype,move,nres;
    double energies[EN_TERMS],etot;
#if defined(PARALLEL) || defined(EXCHANGE)
    mynod=_mynod;
    numnod=_numnod;
#endif
    do_mc=(strcasecmp(command,"run")==0);
    do_energy=(strcasecmp(command,"energy")==0);
    /*tables=NULL;
    en_by_table=NULL;
    frag_nblist=NULL;*/
    printf("Reading control file: %s\n",fname);
    f=fopen(fname,"r");
    if (f==NULL) {
        printf("FATAL ERROR: file %s is not found\n",fname);
        die();
    }
    //The force field
    fgets(buffer,sizeof(buffer),f);
    sscanf(buffer,"%s %s\n",deffname,forcefieldfname);
    trim_string(forcefieldfname);
    printf("Loading parameter file %s.\n",forcefieldfname);
    ffield = new forcefield(forcefieldfname);
    trim_string(deffname);
    printf("Loading definitions file %s.\n",deffname);
    top = new topology(deffname,ffield);
/*#ifndef EXCHANGE //otherwise, this will be done in exchange_init
    if (tables_lambda<0.0) tables_lambda=0.0;
    if (tables_lambda>1.0) tables_lambda=1.0;
    //0 = fully exact, 1 = fully tabulated.
    use_std_tables=(tables_lambda>0.0);
    if ((tables_lambda>0.0) && (tables_lambda<1.0)) {
        printf("Will mix exact/tabulated energies in this simulation.  Fraction %.2f exact and %.2f tabulated.\n",(1.0-tables_lambda),tables_lambda);
    } else if (tables_lambda==0.0) {
        printf("Will calculate all nonbonded interactions exactly.\n");
    } else if (tables_lambda==1.0) {
        printf("Will use noncovalent tables in this simulation.\n");
    }
    if (use_cov_tables) printf("Will use covalent tables for peptide groups in this simulation.\n"); else printf("Will calculate all peptide covalent interactions exactly.\n");
#endif*/
    //Read the sequence and assemble all the topology info.
    //For now, read only one chain.  There is support for multiple chains, but need to work on chain "names", etc.
    //provide support for multiple lines
    seq=read_multiline(f);
    printf("Sequence: %s\n",seq);
    //we need to precount the number of residues so we can set up the all-atom region for residues
    nres=count_words(seq);
    aaregion_res.init(nres);
    fgets(buffer,sizeof(buffer),f);
    aaregion_res.parse_int_list(buffer);
    top->add_segment(' ',seq,aaregion_res); //temporary
    free(seq);
    //top->link_fragments();
    top->create_angle_dihedral_lists(false);
    top->create_non_tab_list(false,&non_tab_list);
    ffield->find_parameters(top->natom,top->atoms);
    top->create_improper_dihedral_lists(false,ffield);
    aaregion_res.print("All atom region residues ");
#ifdef DEBUG_NON_TABULATED
    top->print_detailed_info(aaregion_res);
#else
    top->print_summary_info();
#endif

    //Allocate all the coordinate arrays.
    //oldcenter=(double *) checkalloc(3*top->nfrag,sizeof(double));
    //oldorient=(double *) checkalloc(4*top->nfrag,sizeof(double));
    oldcoords=(double *) checkalloc(3*top->natom,sizeof(double));
    //newcenter=(double *) checkalloc(3*top->nfrag,sizeof(double));
    //neworient=(double *) checkalloc(4*top->nfrag,sizeof(double));
    newcoords=(double *) checkalloc(3*top->natom,sizeof(double));

    //if (do_mc) {
        //Read start/restart and initial structure.
    fgets(buffer,sizeof(buffer),f);
    sscanf(buffer,"%s %s %s\n",word,structfname,struct2fname);
    initcoords=NULL;
    start=(strncasecmp("START",word,5)==0);
    restart=(strncasecmp("RESTART",word,7)==0);
    if (start) {
        //We're starting from a PDB file.
        printf("Will read PDB file %s.\n",structfname);
        initcoords=(double *) checkalloc(3*top->natom,sizeof(double));
        top->read_pdb_file(structfname,initcoords);
        //In place of rebuilding the structure from fragments, need to copy the coordinates from initcoords to oldcoords.
        for (i=0; i<3*top->natom; i++) oldcoords[i]=initcoords[i];
        for (i=0; i<3*top->natom; i++) newcoords[i]=oldcoords[i];
        /*top->assemble_fragments(initcoords,oldcenter,oldorient,oldcoords);
        if (strlen(struct2fname)>0) {
            printf("Writing fitted structure to file %s.\n",struct2fname);
            top->write_pdb_file(struct2fname,oldcoords);
        }*/
        nprevstep=0;
        /*} else if (restart) {
            initcoords=NULL;
#ifdef EXCHANGE
            snprintf(structf3name,sizeof(structf3name),structfname,mynod+1);//1-based replicas
            printf("Will read restart file %s.\n",structf3name);
            read_restart(structf3name);
#else
            printf("Will read restart file %s.\n",structfname);
            read_restart(structfname);
#endif
            printf("Number of previous steps: %ld\n",nprevstep);
            for (ifrag=0; ifrag<top->nfrag; ifrag++) top->update_coords(ifrag,oldcenter,oldorient,oldcoords);
            if (strlen(struct2fname)>0) {
                printf("Reading initial coordinates from PDB file %s.\n",struct2fname);
                initcoords=(double *) checkalloc(3*top->natom,sizeof(double));
                top->read_pdb_file(struct2fname,initcoords);
            }*/
    } else {
        printf("Must specify start or restart in input file.\n");
        die();
    }        //for (ifrag=0; ifrag<top->nfrag; ifrag++) top->copy_frag(ifrag,oldcenter,oldorient,oldcoords,newcenter,neworient,newcoords);
        //read the rest of MC related parameters
    if (do_mc) {

    } //end of if (do_mc)
    fgets(buffer,sizeof(buffer),f);
    sscanf(buffer,"%d %lg   %lg %lg %lg %d\n",&pbc,&boxsize,&cutoff,&listcutoff,&eps,&rdie);
    //Create the go model! (Must do this after coordinates are read.)
    fgets(buffer,sizeof(buffer),f);
    read_go_params(buffer,&go_params);
    print_go_params(go_params);
    go_model = new go_model_info();
    go_model->create_contact_map(top->nres,top->resinfo,top->natom,initcoords,&go_params,aaregion_res);

    //Create nonbond list object.  Set listcutoff=cutoff or less to disable the nb list.
    /*use_nb_list=(listcutoff>cutoff);
    if (use_nb_list) frag_nblist=new fragment_nblist(top->nfrag,listcutoff);*/
    halfboxsize=boxsize*0.5;
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
    //keyword, (relative) probability, and size
    if (do_mc) {
     	fgets(buffer,sizeof(buffer),f);
#ifdef EXCHANGE
        sscanf(buffer,"%ld %ld %ld %ld %ld %d\n",&nmcstep,&nsave,&nprint,&ncheck,&seed, &enwrite);
#else
    	sscanf(buffer,"%ld %ld %ld %ld %lg\n",&nmcstep,&nsave,&nprint,&seed,&mctemp);
#endif
    	if (!restart) { //If restarting, we do not reinitialize the random number generator; read_restart put back the state.
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
        }
        for (i=1; i<=NUM_MOVES; i++) prob[i]=0.0;
    	for(;;) {
            fgets(buffer,sizeof(buffer),f);
            sscanf(buffer,"%s %lg %lg\n",word,&p,&size);
            move=-1;
            for (i=1; i<=NUM_MOVES; i++)
                if (strncasecmp(mc_move_names[i],word,strlen(mc_move_names[i]))==0) move=i;
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
#ifdef EXCHANGE
        fflush(stdout);
        exchange_init(f,fragfmt); //.will read the rest of the input file
    } //if (do_mc)
#else
       printf("Temperature:                 %.2f K\n",mctemp);
       beta=1/(KBOLTZ*mctemp);

       fgets(buffer,sizeof(buffer),f);
       sscanf(buffer,"%s %s %s\n",psffname,xyzfname,restartfname);
       trim_string(psffname);
       top->write_psf_file(psffname,ffield);
       trim_string(xyzfname);
       //trim_string(quatfname);
       trim_string(restartfname);
    } //also if (do_mc)
    //And now... load all the tables.
    /*if (use_std_tables || enwrite) {
        fgets(tablefmt,sizeof(tablefmt),f); //A format string for file names.
        tables = (table * *) checkalloc(top->nfragtypes*top->nfragtypes,sizeof(table *));
        for (i=0; i<top->nfragtypes*top->nfragtypes; i++) tables[i]=NULL;
        top->load_tables(tablefmt,fragfmt,tables);
    } //else tablefmt[0]='\0';
    if (use_cov_tables || enwrite) {
        fgets(covtablefmt,sizeof(covtablefmt),f); //A format string for file names.
        covalent_tables = (covalent_table * *) checkalloc(top->nfragtypes*top->nfragtypes,sizeof(covalent_table *));
        for (i=0; i<top->nfragtypes*top->nfragtypes; i++) covalent_tables[i]=NULL;
        top->load_covalent_tables(covtablefmt,covalent_tables);
    } //else covtablefmt[0]='\0';
    total_table_size=0.0;
    if (use_std_tables || enwrite) for (i=0; i<top->nfragtypes*top->nfragtypes; i++) if (tables[i]!=NULL) total_table_size+=tables[i]->getsize();
    if (use_cov_tables || enwrite) for (i=0; i<top->nfragtypes*top->nfragtypes; i++) if (covalent_tables[i]!=NULL) total_table_size+=covalent_tables[i]->getsize();
    if (use_std_tables || use_cov_tables || enwrite) printf("Total table size:     %.2f MB\n",total_table_size);
    if (enwrite) {
        printf("Will calculate both exact and table-based energies and write them to energy.dat.\n");
        energy_output=fopen("energy.dat","w");
        pairs_output=fopen("pairs.pdb","w");
    }*/
/*#ifdef DEBUG_NON_TABULATED
    if (!do_energy) {
    	total_energy(initcoords,energies,&etot);
    	print_energies(stdout,true,"init:",0,energies,etot);
    }
#endif*/
    //don't need this anymore, we will do the energies in do_energies
    if (do_energy) {
       total_energy(oldcoords,energies,&etot);
       print_energies(stdout,TRUE,"Energy: ",0,energies,etot);
    }
#endif //EXCHANGE
    fclose(f); //Close control file. We're done with it.
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
        total_energy(oldcenter,oldorient,oldcoords,energies,&etot);
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

