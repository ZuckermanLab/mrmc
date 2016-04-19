#include <cstdio>
#include <cctype>
//#include "fragments.h"
#include "topology.h"
#include "ffield.h"
#include "rotations.h"
#include "util.h"

/*
int topology::frag_type_by_name(const char * name)
{
    int itype;
    for (itype=0; itype<nfragtypes; itype++)
        if (strncasecmp(name,fragtypes[itype]->fragname,sizeof(fragtypes[itype]->fragname))==0) return itype;
    return -1; //not found
}

//Tests to see if a file name contains the fragment name.  Used by table's "multiload" constructor to identify the fragment types.
int topology::frag_type_by_file(const char * fname)
{
    int itype;
    char tempname[MAX_FRAGMENT_NAME];
    char tempfile[255];
    strncpy(tempfile,fname,sizeof(tempfile));
    strlower(tempfile);
    for (itype=0; itype<nfragtypes; itype++) {
        strncpy(tempname,fragtypes[itype]->fragname,sizeof(tempname));
        strlower(tempname);
        if (strstr(tempfile,tempname)!=NULL) return itype;
    }
    return -1;
}
*/

int topology::resdefbyname(const char * name)
{
    int iresdef;
    for (iresdef=0; iresdef<nresdef; iresdef++)
        if (strncasecmp(name,resdef[iresdef].name,sizeof(resdef[iresdef].name))==0) return iresdef;
    return -1; //not found
}

//Reads commands from definitions file.
//FRAG fragname filename -- read fragment file
//RESI resname peptide-bond sidechain1,.... -- declares a residue and identifies the fragments of which it is composed
//#... comment line, skip
//Eventually:
//ATOM resname atomname fragmentname atom-in-fragment -- Adds an atom to a residue definition, and associates it with a corresponding atom of a fragment.
//BOND resname atomname atomname ROTATE/NOROTATE -- Identifies two atoms which are bonded and specifies if this bond is one about which rotations should be allowed.
topology::topology(const char * commandfile, forcefield * ffield)
{
    FILE * f;
    char command[255],command2[255],fname[255],name[MAX_FRAGMENT_NAME];
    char * token;
    char * begin;
    const char * delim = " \t\n";
    int i,iatom;
    f=fopen(commandfile,"r");
    /*fragtypes=NULL;
    frags=NULL;*/
    resdef=NULL;
    resinfo=NULL;
    atoms=NULL;
    //nfragtypes=0;
    nresdef=0;
    //nfrag=0;
    natom=0;
    nseg=0;
    nres=0;
    nscrot=0;
    segstart=NULL;
    segend=NULL;
    first_main_chain_frag=NULL;
    chaincodes=NULL;
    iscrot=NULL;
    jscrot=NULL;
    qsystem=0.0;
    /*whichseg=NULL;
    sequence=NULL;*/
    /*itype=0;
    iresdef=0;*/
    if (f==NULL) {
        printf("Could not open definitions file %s\n",commandfile);
        die();
    }
    while (!feof(f)) {
        fgets(command,sizeof(command),f);
        if (feof(f)) break;
        for (i=0; i<strlen(command); i++) {
            command[i]=toupper(command[i]);
            if (command[i]=='#') command[i]='\0';
        }
   //for (i=0; i<strlen(command); i++) command[i]=toupper(command[i]);
        token=strtok(command,delim);
        if (token==NULL) continue; //blank line
        //else if (*token=='#') continue; //comment
/*
        else if (strncmp("FRAG",token,4)==0) { //load a fragment
            //sscanf(token,"%*s %8s %255s\n",name,fname);
            token=strtok(NULL,delim);
            strncpy(name,token,sizeof(name));
            snprintf(fname,sizeof(fname),fragfmt,name);
            strlower(fname); //the fragment file will be all lowercase
            trim_string(fname);
            fragtypes=(fragmenttype **) checkrealloc(fragtypes,nfragtypes+1,sizeof(fragmenttype *));
            fragtypes[nfragtypes]=new fragmenttype(name,fname,ffield);
            nfragtypes++;
        }
*/
         else if (strncmp("RESI",token,4)==0) {
            resdef=(residuedef *)checkrealloc(resdef,nresdef+1,sizeof(residuedef));
            token=strtok(NULL,delim);
            strncpy(resdef[nresdef].name,token,sizeof(resdef[nresdef].name));
            read_residue_definition(f,&resdef[nresdef]);
#ifdef DEBUG_NON_TABULATED
            printf("Processed residue %.4s in definitions file\n",resdef[nresdef].name);
#endif
            /*for (;;) {
                fgets(command,sizeof(command),f);
                if (strncasecmp("atom",command,4)==0) {
                    sscanf("%*s %s %s*/
            nresdef++;
        } else {
            printf("Unrecognized command in definitions file: %s\n",token);
            die();
        }
    }
    fclose(f);
    /*if ((iresdef!=nresdef) || (itype!=nfragtypes)) {
        printf("Did not define all fragment or residue types.\n");
        die();
    }*/
}

//Commands within residue definiton.
//ATOM name [BRANCH]
void topology::read_residue_definition(FILE * f, residuedef * def)
{
    char command[255],command2[255],atomres[6],atomres2[6],atomfrag[6];
    char * token;
    const char * delim = " \t\n";
    int i,is_branch,itype,iatom,iatomfrag;
    def->natom=0;
    def->nbond=0;
    //def->nfrag=0;
    //def->branchname[0]='\0';
    def->branchatom=-1;
    while (!feof(f)) {
        fgets(command,sizeof(command),f);
        for (i=0; i<strlen(command); i++) {
            command[i]=toupper(command[i]);
            if (strcmp(def->name,"TRP")==0) {
                printf(""); //debug lander
            }
            if (command[i]=='#') command[i]='\0';
        }
        strncpy(command2,command,sizeof(command2));
        token=strtok(command2,delim);
        if (token==NULL) continue; //blank line
        //else if (*token=='#') continue; //comment
        else if (strncmp("END",token,3)==0) return; //end of the residue definition
        else if (strncmp("ATOM",token,4)==0) {
            //sscanf(command,"%*s %s",atomnames[natom]);
            token=strtok(NULL,delim);
            strncpy(def->atomnames[def->natom],token,sizeof(def->atomnames[def->natom]));
            token=strtok(NULL,delim); //read atom type
            sscanf(token,"%d",&def->atomtypes[def->natom]);
            if (strstr(command,"BRANCH")!=NULL) {
                if (def->branchatom>0) {
                    printf("You specified more than one branch atom name for residue %s\n",def->name);
                    die();
                }
                def->branchatom=def->natom;
            }
            /*if (strlen(def->branchname)>0) {

            } else strncpy(def->branchname,def->atomnames[def->natom],sizeof(def->branchname));*/
            def->natom++;
        } else if (strncmp("BOND",token,4)==0) {
            //sscanf(command,"%*s %s %s",iname[nbond],namebuf);
            token=strtok(NULL,delim);
            strncpy(def->iname[def->nbond],token,sizeof(def->iname[def->nbond]));
            token=strtok(NULL,delim);
            if (token[0]=='-') {
                def->joffset[def->nbond]=-1;
                strncpy(def->jname[def->nbond],token+1,sizeof(def->jname[def->nbond]));
            } else {
                def->joffset[def->nbond]=0;
                strncpy(def->jname[def->nbond],token,sizeof(def->jname[def->nbond]));
            }
            if (strstr(command,"ROTATABLE")!=NULL) {
                if (strstr(command,"BACKBONE")!=NULL) {
                    def->rottype[def->nbond]=RT_BACKBONE;
                } else if (strstr(command,"SIDECHAIN")!=NULL) {
                    def->rottype[def->nbond]=RT_SIDECHAIN;
                } else {
                    printf("You must specify whether the rotatable bond in residue %s involving atoms %s and %s is part of the SIDECHAIN or BACKBONE.\n",
                        def->name,def->iname[def->nbond],def->jname[def->nbond]);
                    die();
                }
            } else def->rottype[def->nbond]=RT_NOT_ROTATABLE;
            def->nbond++;
/*
           //Read fragment to residue defintiion.
            token=strtok(NULL,delim);
            itype=frag_type_by_name(token);
            if (itype<0) {
                printf("Unrecognized fragment type %s in residue %s\n",token,def->name);
                die();
            }
            def->frags[def->nfrag].fragtype=itype;
            if (strstr(command,"SIDE")!=NULL) def->frags[def->nfrag].is_side_chain=TRUE;
            else if (strstr(command,"BACK")!=NULL) def->frags[def->nfrag].is_side_chain=FALSE;
            else {
                printf("You must specify whether a fragment in residue %s is a sidechain or backbone fragment\n",def->name);
                die();
            }
            for (iatom=0; iatom<fragtypes[itype]->natom; iatom++) def->frags[def->nfrag].atomnames[iatom][0]='\0';
            //Begin reading the atom-name fragment-name pairs.
            for (;;) {
                if (feof(f)) {
                    printf("unexpected end of file in definitions file\n");
                    die();
                }
                fgets(command,sizeof(command),f);
                if (strstr(command,"END")!=NULL) break;
                sscanf(command,"%s %s\n",atomres,atomfrag);
                //find atom within fragment
                iatomfrag=fragtypes[itype]->findatombyname(atomfrag);
                if (iatomfrag<0) {
                    printf("unrecognized fragment atom in fragment definition in residue %s: %s\n",def->name,atomfrag);
                    die();
                }
                //check for '-' syntax to designate previous residue
                if (atomres[0]=='-') {
                    def->frags[def->nfrag].offset[iatomfrag]=-1;
                    strncpy(def->frags[def->nfrag].atomnames[iatomfrag],atomres+1,sizeof(def->atomnames[iatomfrag]));
                } else {
                    def->frags[def->nfrag].offset[iatomfrag]=0;
                    strncpy(def->frags[def->nfrag].atomnames[iatomfrag],atomres,sizeof(def->atomnames[iatomfrag]));
                }

            }
            for (iatom=0; iatom<fragtypes[itype]->natom; iatom++) if(strlen(def->frags[def->nfrag].atomnames[iatom])==0) {
                printf("You failed to specify a corresponding residue atom for fragment atom %s of fragment type %s in residue %s\n",
                    fragtypes[itype]->names[iatom],fragtypes[itype]->fragname,def->name);
                die();
            }
            def->nfrag++;
*/
        } else if (strncmp("FRAG",token,4)==0) { //so that we can ignore fragment definitions
            for (;;) {
                if (feof(f)) {
                    printf("unexpected end of file in definitions file\n");
                    die();
                }
                fgets(command,sizeof(command),f);
                if (strstr(command,"END")!=NULL) break;
            }

        } else {
            printf("Unrecognized command in definitions file: %s\n",token);
            die();
        }
    }
    //if we got here, it must be because of an unexpected end of file
    if (feof(f)) {
        printf("unexpected end of file in definitions file\n");
        die();
    }
}


void topology::print_detailed_info(subset aaregion_res)
{
    int iatom,iactualatom,ibond,jatom,ifrag,jfrag,ires,itype;
    //Residue info:
    printf("residue number, segment, res. type (res. name), branch atom number, start, end, in AA region?\n");
    for (ires=0; ires<nres; ires++) {
        printf("%d, %d, %d (%.4s), %d, %d, %d, %c\n",ires,resinfo[ires].whichseg,resinfo[ires].restype,
        resdef[resinfo[ires].restype].name,resinfo[ires].branchatom, resinfo[ires].atomstart,resinfo[ires].atomend, yesno(aaregion_res[ires]));
        for (ibond=0; ibond<resinfo[ires].nbbrot; ibond++) printf("Backbone bond %d-%d is rotatable.\n",resinfo[ires].ibbrot[ibond],resinfo[ires].jbbrot[ibond]);

        printf("\n");
    }
    printf("\n");
    for (ibond=0; ibond<nscrot; ibond++) printf("Side chain bond %d-%d is rotatable.\n",iscrot[ibond],jscrot[ibond]);
    printf("\n");
    //Fragment info:
/*
    for (ifrag=0; ifrag<nfrag; ifrag++) {
        printf("Fragment %d: \n",ifrag);
        itype=frags[ifrag].type;
        printf("Type %d \"%s\"\n",itype,fragtypes[itype]->fragname);
        printf("Side chain?  %c\n",yesno(frags[ifrag].is_side_chain));
        if (!frags[ifrag].is_side_chain) printf("N, CA, and C atoms: %d %d %d\n",frags[ifrag].atn,frags[ifrag].atca,frags[ifrag].atc);
        printf("Main chain previous/next: %d %d\n",frags[ifrag].main_chain_prev,frags[ifrag].main_chain_next);
        printf("Side chain previous/next: %d %d\n",frags[ifrag].side_chain_prev,frags[ifrag].side_chain_next);
        printf("Atoms -- index in fragment, name in fragment, index in system, res. number, res.name,  name in system\n");
        for (iatom=0; iatom<fragtypes[itype]->natom; iatom++) {
            iactualatom=frags[ifrag].atoms[iatom];
            printf("%d %.4s %d %d %s %s\n",iatom,fragtypes[itype]->names[iatom],iactualatom,atoms[iactualatom].resNum,atoms[iactualatom].resName,atoms[iactualatom].name);
        }
        printf("\n");
    }
*/
    //Atom info
    printf("index, res. number, atom name, type, class, list of bonded atoms\n");
    for (iatom=0; iatom<natom; iatom++) {
        printf("%d, %d, %s, %d, %d, ",iatom,atoms[iatom].resNum,atoms[iatom].name,atoms[iatom].type,atoms[iatom].classx);
        for (ibond=0; ibond<atoms[iatom].numOfBondedAtoms; ibond++) printf("%d ",atoms[iatom].bondedAtomList[ibond]);
        printf("\n");
    }

}

void topology::print_summary_info(void)
{
    //int ires, naares, iatom, naa_atom;
    //for (ires=0; ires<nres; ires++) if (aaregion_res[ires]) naares++;
    printf("Number of segments:                     %d\n",nseg);
    printf("Number of residues:                     %d\n",nres);
    printf("Number of atoms:                        %d\n",natom);

    printf("Number of side chain rotatable bonds: %d\n",nscrot);
    printf("Total system charge:                  %.2f\n",qsystem);
}

void topology::insert_residue(const char * res, subset aaregion_res)
{
    int nnewfrag,nnewatom,nnewscrot,ifrag,restype,iatom,jatom,ibond,itype,ires,iactualatom,nfragbonded;
    //int fragbonded[6]; //max bonds per atom, should be a constant
    //fragmenttype * ftype;
    restype=resdefbyname(res);
    if (restype<0) {
        printf("Unknown residue type %s\n",res);
        die();
    }
    //Current segment = nseg-1.\
    //New residue nres.
    nnewatom=resdef[restype].natom;
    //nnewfrag=resdef[restype].nsidechain+resdef[restype].nmainchain;
    //Extend "whichseg" and "sequence" arrays by 1, and fill new element in.
    resinfo=(reslookup *) checkrealloc(resinfo,nres+1,sizeof(reslookup));
    resinfo[nres].whichseg=nseg;
    resinfo[nres].restype=restype;
    resinfo[nres].atomstart=natom;
    resinfo[nres].atomend=natom+nnewatom-1;
    resinfo[nres].nbbrot=0;
    //resinfo[nres].nscrot=0;
    //Insert all atoms into array.
    atoms=(ATOMS *) checkrealloc(atoms,natom+nnewatom,sizeof(ATOMS));
    for (iatom=0; iatom<nnewatom; iatom++) {
        strncpy(atoms[natom].name,resdef[restype].atomnames[iatom],sizeof(atoms[natom].name));
        strncpy(atoms[natom].resName,res,sizeof(atoms[natom].resName));
        atoms[natom].resNum=nres;
        atoms[natom].numOfBondedAtoms=0;
        atoms[natom].fragment=-1;
        atoms[natom].fragatom=-1;
        atoms[natom].type=resdef[restype].atomtypes[iatom];//The atom class will be identified by forcefield::find_parameters.
        atoms[natom].is_in_aa_region=(aaregion_res[nres]);//this may be redundant but is convenient
        atoms[natom].is_backbone=((strncmp(atoms[natom].name,"N",sizeof(atoms[natom].name))==0) ||
                                    (strncmp(atoms[natom].name,"CA",sizeof(atoms[natom].name))==0) ||
                                    (strncmp(atoms[natom].name,"C",sizeof(atoms[natom].name))==0));
        natom++;
    }
/*
    //Insert all fragments.
    //How many new fragments will there be?  Insert only fragments that are unique to this residue or that overlap with the previous.
    //We'll get to fragments that overlap with next on the next go around.
    //We have to do something special at the beginning and end of the sequence but I'm not sure what yet.
    //nnewfrag=0;
    //for (ifrag=0; ifrag<resdef[restype].nfrag; ifrag++) if (!(resdef[restype].frags[ifrag].overlap)) nnewfrag++;
    nnewfrag=resdef[restype].nfrag;
    //Extend fragments array by nnewfrag.
    frags=(fragment *)checkrealloc(frags,nfrag+nnewfrag,sizeof(fragment));
    for (ifrag=0; ifrag<nnewfrag; ifrag++) {
        itype=resdef[restype].frags[ifrag].fragtype;
        frags[nfrag].type=itype;
        fragtypes[itype]->n_used++;
        qsystem+=fragtypes[itype]->qtot;
        frags[nfrag].is_side_chain=resdef[restype].frags[ifrag].is_side_chain;
        //keep track of the first main chain fragment, so we can walk along peptide backbone
        if ((!frags[nfrag].is_side_chain) && (first_main_chain_frag[nseg]<0)) first_main_chain_frag[nseg]=nfrag;
        //For now, set all the linking indices to -1.
        frags[nfrag].main_chain_next=-1;
        frags[nfrag].main_chain_prev=-1;
        frags[nfrag].side_chain_next=-1;
        frags[nfrag].side_chain_prev=-1;
        frags[nfrag].atn=-1;
        frags[nfrag].atca=-1;
        frags[nfrag].atc=-1;     //ftype=fragtypes[itype];
        //For each atom within this fragment (that intersects the residue, figure out the actual residue number of the atom,
        //and the actual number of the atom within the whole system.  Use this to construct an atom assignment for each fragment,
        //and also to set the fragment field in each atom structure.
        //This will only work if the "offset" is less than zero, otherwise we haven't added the residue yet.
        for (iatom=0; iatom<fragtypes[itype]->natom; iatom++) {
            ires=nres+resdef[restype].frags[ifrag].offset[iatom]; //This is the actaul residue number in which this atom can be found.
            if (ires<0) {
                printf("You cannot specify a previous residue on the first residue of the chain.\n");
                die();
            }
            //Find the atom with this name in residue ires (either this or the previous residue)
            iactualatom=find_atom(ires,resdef[restype].frags[ifrag].atomnames[iatom]);
            if (iactualatom<0) {
                printf("Failed to find residue atom %s in residue %s %d as a partner for fragment atom %s in fragment of type %s\n",
                    resdef[restype].frags[ifrag].atomnames[iatom],resdef[restype].name,nres+1,fragtypes[itype]->names[iatom],fragtypes[itype]->fragname);
                die();
            }
            frags[nfrag].atoms[iatom]=iactualatom;
            atoms[iactualatom].fragment=nfrag;
            atoms[iactualatom].fragatom=iatom;
            atoms[iactualatom].type=fragtypes[itype]->types[iatom];
            atoms[iactualatom].mass=fragtypes[itype]->mass[iatom];
            atoms[iactualatom].is_side_chain=frags[nfrag].is_side_chain;
            if (!frags[nfrag].is_side_chain) {
                if (strcmp(fragtypes[itype]->names[iatom],"N")==0) frags[nfrag].atn=iactualatom;
                if (strcmp(fragtypes[itype]->names[iatom],"CA")==0) frags[nfrag].atca=iactualatom;
                if (strcmp(fragtypes[itype]->names[iatom],"C")==0) frags[nfrag].atc=iactualatom;
            }
        }
        nfrag++; //finished with this fragment
    }
*/
    //Insert bonding info from the residue definition.
    for (ibond=0; ibond<resdef[restype].nbond; ibond++) {
        //Find the two atoms which are bonded. "joffset" must be less than 0, otherwise we havent' added
        iatom=find_atom(nres,resdef[restype].iname[ibond]);
        if (iatom<0) {
            printf("Failed to find bonded residue atom %s in residue %s %d\n",
                    resdef[restype].iname[ibond],resdef[restype].name,nres);
            die();
        }
        ires=nres+resdef[restype].joffset[ibond]; //the residue to which the bonded atom belongs
        if (ires<0) {
            printf("You cannot specify a previous residue on the first residue of the chain.\n");
            die();
        }
        jatom=find_atom(ires,resdef[restype].jname[ibond]);
        if (jatom<0) {
            printf("Failed to find bonding partner %s for atom %s in residue %s %d while adding residue %s %d\n",
                    resdef[restype].jname[ibond],resdef[restype].iname[ibond],resdef[resinfo[ires].restype].name,ires,resdef[restype].name,nres);
            die();
        }
        //Only add backbone bonds if not in the all-atom region.
        if (aaregion_res[nres] || (atoms[iatom].is_backbone && atoms[jatom].is_backbone)) {
            atoms[iatom].bondedAtomList[atoms[iatom].numOfBondedAtoms]=jatom;
            atoms[iatom].numOfBondedAtoms++;
            atoms[jatom].bondedAtomList[atoms[jatom].numOfBondedAtoms]=iatom;
            atoms[jatom].numOfBondedAtoms++;
        }
        //If the bond is rotatable, include it in appropriate list.  This information is used to generate Monte Carlo moves.
        if (resdef[restype].rottype[ibond]==RT_BACKBONE) {
            resinfo[nres].ibbrot[resinfo[nres].nbbrot]=iatom;
            resinfo[nres].jbbrot[resinfo[nres].nbbrot]=jatom;
            resinfo[nres].nbbrot++;
        } else if (aaregion_res[nres] && resdef[restype].rottype[ibond]==RT_SIDECHAIN) { //only add side chains if in the all-atom region
            /*resinfo[nres].iscrot[resinfo[nres].nscrot]=iatom;
            resinfo[nres].jscrot[resinfo[nres].nscrot]=jatom;*/
            iscrot=(int *)checkrealloc(iscrot,nscrot+1,sizeof(int));
            jscrot=(int *)checkrealloc(jscrot,nscrot+1,sizeof(int));
            iscrot[nscrot]=iatom;
            jscrot[nscrot]=jatom;
            nscrot++;
        }

        /*if (resdef[restype].rotatable[ibond]) {
            resinfo[nres].irotatable[resinfo[nres].nrotatable]=iatom;
            resinfo[nres].jrotatable[resinfo[nres].nrotatable]=jatom;
            //The bond is part of the side chain if either of the atoms in it is.
            resinfo[nres].isscrotatable[resinfo[nres].nrotatable]=(atoms[iatom].is_side_chain || atoms[jatom].is_side_chain);
            resinfo[nres].nrotatable++;
            if (resinfo[nres].isscrotatable[resinfo[nres].nrotatable]) resinfo[nres].nscrotatable++; else resinfo[nres].nbbrotatable++;
        }*/
    }
    //this was a bug in the tablemc-proteins code as well -- this condition was inside the "ibond" loop,
    //so that the CB-CD bond was added once for every bond in the pro residue, not just once as it should have been.
    if (aaregion_res[nres] && strcasecmp(resdef[restype].name,"PRO")==0) {
        iatom=find_atom(nres,"CB");
        jatom=find_atom(nres,"CD");
        iscrot=(int *)checkrealloc(iscrot,nscrot+1,sizeof(int));
        jscrot=(int *)checkrealloc(jscrot,nscrot+1,sizeof(int));
        iscrot[nscrot]=iatom;
        jscrot[nscrot]=jatom;
        nscrot++;
    }
    resinfo[nres].branchatom=resinfo[nres].atomstart+resdef[restype].branchatom;
    segend[nseg]=nres;
    //Need to construct atoms based on the definition.
    nres++;
}


int topology::find_atom(int actual_res, const char * aname)
{
    int iatom;
    for (iatom=resinfo[actual_res].atomstart; iatom<=resinfo[actual_res].atomend; iatom++)
        if (strncasecmp(atoms[iatom].name,aname,sizeof(atoms[iatom].name))==0) return iatom;
    return -1; //not found
}

//res is 1-based here.
int topology::find_atom(char chain, int res, const char * aname)
{
    int iseg,foundseg;
    foundseg=-1;
    for (iseg=0; iseg<nseg; iseg++) if (chaincodes[iseg]==chain) foundseg=iseg;
    if (foundseg<0) return -1; else return find_atom(segstart[foundseg]+res-1,aname);

}

/*
void topology::link_fragments(void)
{
    int iatom,jatom,ifrag,jfrag,ibond,imainchain,isidechain,branchatom;
    int fragbonded[6];
    for (iatom=0; iatom<natom; iatom++) {
        ifrag=atoms[iatom].fragment;
        if (ifrag<0) {
            printf("Atom %s in residue %d not in any fragment\n",atoms[iatom].name,atoms[iatom].resNum);
            die();
        }
        for (ibond=0; ibond<atoms[iatom].numOfBondedAtoms; ibond++) {
            jatom=atoms[iatom].bondedAtomList[ibond];
            jfrag=atoms[jatom].fragment;
            //This assumes the fragments are added from beginning to end of the chain and down each side chain.
            if (jfrag==ifrag) continue; //This bond is within the fragment.
            if (frags[ifrag].is_side_chain && frags[jfrag].is_side_chain) { //Both part of the side chain.
                if (jfrag>ifrag) {
                    frags[ifrag].side_chain_next=jfrag;
                    frags[jfrag].side_chain_prev=ifrag;
                } else {
                    frags[ifrag].side_chain_prev=jfrag;
                    frags[jfrag].side_chain_next=ifrag;
                }
            } else if ((!frags[ifrag].is_side_chain) && (!frags[jfrag].is_side_chain)) { //Both part of the main chain.
                if (jfrag>ifrag) {
                    frags[ifrag].main_chain_next=jfrag;
                    frags[jfrag].main_chain_prev=ifrag;
                } else {
                    frags[ifrag].main_chain_prev=jfrag;
                    frags[jfrag].main_chain_next=ifrag;
                }
            } else { //Bond from main-chain to side-chain
                if (frags[ifrag].is_side_chain) {
                    imainchain=jfrag;
                    isidechain=ifrag;
                } else {
                    imainchain=ifrag;
                    isidechain=jfrag;
                }
                frags[imainchain].side_chain_next=isidechain;
                frags[isidechain].side_chain_prev=imainchain;
            }
        }
    }
}

int topology::is_bonded(int ifrag, int jfrag)
{
    int result;
    result=(frags[ifrag].main_chain_next==jfrag) || (frags[ifrag].main_chain_prev==jfrag) || (frags[ifrag].side_chain_next==jfrag) || (frags[ifrag].side_chain_prev==jfrag);
    result=result || (frags[jfrag].main_chain_next==ifrag) || (frags[jfrag].main_chain_prev==ifrag) || (frags[jfrag].side_chain_next==ifrag) || (frags[jfrag].side_chain_prev==ifrag);
    return result;

}
*/

void topology::add_segment(char chain, const char * sequence, subset aaregion_res)
{
    char * token;
    char * buf;
    const char * delim = " \t\n";
    //Allocate enough space for all those new residues.
    segstart=(int *)checkrealloc(segstart,(nseg+1),sizeof(int));
    segend=(int *)checkrealloc(segend,(nseg+1),sizeof(int));
    first_main_chain_frag=(int *)checkrealloc(first_main_chain_frag,(nseg+1),sizeof(int));
    first_main_chain_frag[nseg]=-1;
    chaincodes=(char *) checkrealloc(chaincodes,(nseg+1),sizeof(char));
    segstart[nseg]=nres;
    chaincodes[nseg]=chain;
    buf=(char*)checkalloc(strlen(sequence)+1,sizeof(char));
    strcpy(buf,sequence);
    buf[strlen(sequence)]='\0';
    token=strtok(buf,delim);
    while (token!=NULL) {
        insert_residue(token,aaregion_res);
        token=strtok(NULL,delim);
    }
    segend[nseg]=nres-1;
    nseg++;
    free(buf);
}

/*
//calls fit_all_fragments, prints diagnostic messages
void topology::assemble_fragments(double * orig_coords, double * center, double * orient, double * new_coords)
{
    double * mass;
    double * rmsds;
    int iatom,ifrag,k;
    double rmsd,rmsd2,disp[3],q[4];
    printf("Assembling structure from fragments by RMSD fitting.\n");
    rmsds=(double *) checkalloc(nfrag,sizeof(double));
    fit_all_fragments(orig_coords,center,orient,new_coords,rmsds);
    for (ifrag=0; ifrag<nfrag; ifrag++) printf("Fitted fragment %d (type %s) to reference geometry. Mass-weighted RMSD = %.3f A\n",ifrag,fragtypes[frags[ifrag].type]->fragname,rmsds[ifrag]);
    mass=(double *) checkalloc(natom,sizeof(double));
    for (iatom=0; iatom<natom; iatom++) mass[iatom]=atoms[iatom].mass;
    //Find mass-weighted RMSD over whole system.
    rmsd_fit(natom,mass,orig_coords,new_coords,disp,q,&rmsd);
    //find heavy atom RMSD
    for (iatom=0; iatom<natom; iatom++) if (mass[iatom]<2.0) mass[iatom]=0.0; else mass[iatom]=1.0;
    rmsd_fit(natom,mass,orig_coords,new_coords,disp,q,&rmsd2);
    printf("Assembly complete.\n");
    printf("Overall mass-weighted RMSD = %.3f A\n",rmsd);
    printf("Heavy atom RMSD            = %.3f A\n",rmsd2);
    free(mass);
    free(rmsds);
}

//Fits all the fragments. to the coordinates specified by orig_coords.
//Reconstructs the coordinates from the fragments, placing into new coordinates
//This subroutine assumes all the coordinates are known.  Probably desirable to fill in hydrogen atoms/minimize with an MD program
//prior to setting up in here.
//A quiet routine.
void topology::fit_all_fragments(double * orig_coords, double * center, double * orient, double * new_coords, double * rmsds)
{
    int ifrag,iatom,iactualatom,nknownatoms,k;
    fragmenttype * fragtype;
    double tempcoords[3*MAX_ATOMS_PER_FRAGMENT];
    for (ifrag=0; ifrag<nfrag; ifrag++) {
        fragtype=fragtypes[frags[ifrag].type];
        //Rearrange atoms from coords into the fragment order inside tempcoords.
        //Also figure out how many atoms have known coordinates.
        for (iatom=0; iatom<fragtype->natom; iatom++) {
            iactualatom=frags[ifrag].atoms[iatom];
            for (k=0; k<3; k++) tempcoords[3*iatom+k]=orig_coords[3*iactualatom+k];
            //mass[iactualatom]=atoms[iactualatom].mass;
        }
        fragtype->fit_fragment(tempcoords,&center[3*ifrag],&orient[4*ifrag],&rmsds[ifrag]);
        update_coords(ifrag,center,orient,new_coords);
    }
}

void topology::update_coords(int ifrag, double * center, double * orient, double * coords)
{
    int iatom,iactualatom,k;
    fragmenttype * fragtype;
    double tempcoords[3*MAX_ATOMS_PER_FRAGMENT];
    fragtype=fragtypes[frags[ifrag].type];
    fragtype->get_coords(&center[3*ifrag],&orient[4*ifrag],tempcoords);
    for (iatom=0; iatom<fragtype->natom; iatom++) {
        iactualatom=frags[ifrag].atoms[iatom];
        for (k=0; k<3; k++) coords[3*iactualatom+k]=tempcoords[3*iatom+k];
    }
}
//Copy fragment i's center, orientation, and coordinates from "1" to "2".
void topology::copy_frag(int ifrag, double * center1, double * orient1, double * coords1, double * center2, double * orient2, double * coords2)
{
    int iatom,k,iactualatom;
    fragmenttype * fragtype;
    for (k=0; k<3; k++) center2[3*ifrag+k]=center1[3*ifrag+k];
    for (k=0; k<4; k++) orient2[4*ifrag+k]=orient1[4*ifrag+k];
    fragtype=fragtypes[frags[ifrag].type];
    for (iatom=0; iatom<fragtype->natom; iatom++) {
        iactualatom=frags[ifrag].atoms[iatom];
        for (k=0; k<3; k++) coords2[3*iactualatom+k]=coords1[3*iactualatom+k];
    }
}

//There is a covalent table for each pair of peptide fragment type.
bool topology::use_covalent_table(int itype, int jtype)
{
    //return ((strstr(fragtypes[itype]->fragname,"PEPTIDE")!=NULL) && (strstr(fragtypes[jtype]->fragname,"PEPTIDE")!=NULL));
    return (fragtypes[itype]->has_covalent_tables) && (fragtypes[jtype]->has_covalent_tables);
}

//A term is represented in the covalent tables if it represents exactly two fragments, and there is a covalent table for the two fragments

bool topology::term_in_covalent_tables(int iatom, int jatom)
{
    return ((atoms[iatom].fragment!=atoms[jatom].fragment) && use_covalent_table(frags[atoms[iatom].fragment].type,frags[atoms[jatom].fragment].type));
}

bool topology::term_in_covalent_tables(int iatom, int jatom, int katom)
{
    bool flag;
    flag=((atoms[iatom].fragment==atoms[jatom].fragment) && (atoms[katom].fragment!=atoms[jatom].fragment)
        && use_covalent_table(frags[atoms[jatom].fragment].type,frags[atoms[katom].fragment].type));
    flag=flag || ((atoms[katom].fragment==atoms[jatom].fragment) && (atoms[iatom].fragment!=atoms[jatom].fragment)
        && use_covalent_table(frags[atoms[jatom].fragment].type,frags[atoms[katom].fragment].type));
    return flag;
}

bool topology::term_in_covalent_tables(int iatom, int jatom, int katom)
{
    int count, i,frag[3],types[3];
    count=1;
    frag[0]=atoms[iatom].fragment;
    frag[1]=atoms[jatom].fragment;
    frag[2]=atoms[katom].fragment;
    if (frag[1]!=frag[0]) count++;
    if ((frag[2]!=frag[0]) && (frag[2]!=frag[1])) count++;
    if (count!=2) return false;
    for (i=0; i<3; i++) types[i]=frags[frag[i]].type;
    if (frag[1]!=frag[0]) return use_covalent_table(types[0],types[1]);
    if (frag[2]!=frag[1]) return use_covalent_table(types[1],types[2]);
    printf("Error in covalent term test.\n"); //shouldn't get here.
    die();
    return false; //definitely won't get here.
}

bool topology::term_in_covalent_tables(int iatom, int jatom, int katom, int latom)
{
    int count, i,frag[4],types[4];
    count=1;
    frag[0]=atoms[iatom].fragment;
    frag[1]=atoms[jatom].fragment;
    frag[2]=atoms[katom].fragment;
    frag[3]=atoms[latom].fragment;
    if (frag[1]!=frag[0]) count++;
    if ((frag[2]!=frag[0]) && (frag[2]!=frag[1])) count++;
    if ((frag[3]!=frag[0]) && (frag[3]!=frag[1]) && (frag[3]!=frag[2])) count++;
    if (count!=2) return false;
    for (i=0; i<4; i++) types[i]=frags[frag[i]].type;
    if (frag[1]!=frag[0]) return use_covalent_table(types[0],types[1]);
    if (frag[2]!=frag[1]) return use_covalent_table(types[1],types[2]);
    if (frag[3]!=frag[2]) return use_covalent_table(types[2],types[3]);
    printf("Error in covalent term test.\n"); //shouldn't get here.
    die();
    return false; //can't happen
}
*/

//Fills the angleAtomList and bonded14AtomList fields in the atoms.  Also looks for fragments that are 1-3 or 1-4 and sets the closefragments flags
void topology::create_angle_dihedral_lists(bool using_cov_tables)
{
    int iatom, jatom, katom, matom, oatom, i, j, k, m,num,ifrag,jfrag,kfrag,mfrag;
    int iclass,jclass,kclass,mclass;
    int flag2,ii;
    //closefragments = (bool *) checkalloc(nfrag*nfrag,sizeof(bool));
    //for (i=0; i<nfrag*nfrag; i++) closefragments[i]=false;
    for(iatom=0;iatom<natom;iatom++){
        atoms[iatom].numOfAngles=0;
        //atoms[iatom].numOfBonded12Atoms=0;
        //atoms[iatom].numOfAngles=0;
        atoms[iatom].numOfImprops=0;
        atoms[iatom].numOfBonded14Atoms=0;
        //atoms[iatom].numOfNonBondedAtoms=0;
    }
    for(iatom=0;iatom<natom;iatom++){
        for(j=0;j<atoms[iatom].numOfBondedAtoms;j++){
            jatom=atoms[iatom].bondedAtomList[j];
/*
      ifrag=atoms[iatom].fragment;
      jfrag=atoms[jatom].fragment;
      closefragments[ifrag*nfrag+jfrag]=true;
      closefragments[jfrag*nfrag+ifrag]=true;
*/
            for(k=0;k<atoms[jatom].numOfBondedAtoms;k++){
                katom=atoms[jatom].bondedAtomList[k];
                if(katom!=iatom){
/*
	    if (using_cov_tables && term_in_covalent_tables(iatom,jatom,katom)){
#ifdef DEBUG_NON_TABULATED
            printf("Skipping angle term involving atoms %d %d %d\n",iatom,jatom,katom);
#endif
            continue;
	    }
*/
        //regular angles
                    num = atoms[iatom].numOfAngles;
                    atoms[iatom].angleAtomList[3*num]   = iatom;
                    atoms[iatom].angleAtomList[3*num+1] = jatom;
                    atoms[iatom].angleAtomList[3*num+2] = katom;
                    atoms[iatom].numOfAngles++;
                    kfrag=atoms[katom].fragment;
	  /*closefragments[ifrag*nfrag+kfrag]=true;
	  closefragments[kfrag*nfrag+ifrag]=true;*/
                }
            }
        }
    }
  /*
  //print out bonded 1-3 atoms
  for(iatom=0;iatom<numOfAtoms;iatom++){
    printf("atom: %d, has bonded 1-3 atoms: ",iatom);
    for(j=0;j<atoms[iatom].numOfAngles;j++){
      printf("%d %d %d   ",atoms[iatom].angleAtomList[3*j],atoms[iatom].angleAtomList[3*j+1],atoms[iatom].angleAtomList[3*j+2]);
    }
    printf("\n");
  }
  */

    //now take care of 1-4 bonded atoms
  for(iatom=0;iatom<natom;iatom++){
    for(j=0;j<atoms[iatom].numOfBondedAtoms;j++){
      jatom=atoms[iatom].bondedAtomList[j];
      for(k=0;k<atoms[jatom].numOfBondedAtoms;k++){
        katom=atoms[jatom].bondedAtomList[k];
        if(katom!=iatom){
            for(m=0;m<atoms[katom].numOfBondedAtoms;m++){
                matom=atoms[katom].bondedAtomList[m];
                if(matom!=jatom){
                            //check for Pro ring if bonded1-4 is actually 1-3
                            //in Pro ring bonded1-4 are used for dihedrals but not for Coulomb/LJ
                            //therefore set bonded1-4Flag to 0 only for Pro ring atoms, so that it won't calculate Coulomb/LJ for them
                            /*flag2=1;
                            for(ii=0;ii<atoms[iatom].numOfAngles;ii++)
                            {
                                oatom=atoms[iatom].angleAtomList[3*ii+2];
                                if(matom==oatom)
                                {
                                    flag2=0;
                                    break;
                                }
                            }
                            // for 6-ring systems check if the last atom in the quartet is already present in the 1-4 list
                            // this is needed for correct vdw and coulomb to avoid double counting interactions
                            if(flag2==1)
                            {
                                for(ii=0;ii<atoms[iatom].numOfBonded14Atoms;ii++)
                                {
                                    oatom=atoms[iatom].bonded14AtomList[4*ii+3];
                                    if(matom==oatom)
                                    {
                                        flag2=0;
                                        break;
                                    }
                                }
                            }
	    if (flag2==0) continue;*/
/*
            if (using_cov_tables && term_in_covalent_tables(iatom,jatom,katom,matom)){
#ifdef DEBUG_NON_TABULATED
                printf("Skipping dihedral term involving atoms %d %d %d %d\n",iatom,jatom,katom,matom);
#endif
                continue;
            }
*/
            num = atoms[iatom].numOfBonded14Atoms;
            atoms[iatom].bonded14AtomList[4*num]   = iatom;
            atoms[iatom].bonded14AtomList[4*num+1] = jatom;
            atoms[iatom].bonded14AtomList[4*num+2] = katom;
            atoms[iatom].bonded14AtomList[4*num+3] = matom;
            atoms[iatom].numOfBonded14Atoms++;
/*
            ifrag=atoms[iatom].fragment;
            mfrag=atoms[matom].fragment;
            closefragments[ifrag*nfrag+mfrag]=true;
            closefragments[mfrag*nfrag+ifrag]=true;
*/
	    }
	  }
	}
      }
    }
  }


}

void topology::create_improper_dihedral_lists(bool using_cov_tables, forcefield * ffield)
{
    int i,j,k,m,num,iatom, jatom, katom, matom,iclass,jclass,kclass,mclass;
     // identify improper dihedrals
    for(iatom=0;iatom<natom;iatom++){
    if(/*(!(atoms[iatom].is_side_chain)) &&*/ atoms[iatom].numOfBondedAtoms==3){// identify backbone triangular centers
      iclass = atoms[iatom].classx;
      for(j=0;j<=2;j++){// generate all possible permutations of three terminal atoms
	for(k=0;k<=2;k++){
	  if(k!=j){
	    for(m=0;m<=2;m++){
	      if(m!=j && m!=k){
		jatom = atoms[iatom].bondedAtomList[j];
		katom = atoms[iatom].bondedAtomList[k];
		matom = atoms[iatom].bondedAtomList[m];
		//printf("j: %d, k: %d, m: %d, iatom: %d, jatom: %d, katom: %d, matom: %d\n",j,k,m,iatom,jatom,katom,matom);

        /*if (using_cov_tables && term_in_covalent_tables(iatom,jatom,katom,matom)) {
#ifdef DEBUG_NON_TABULATED
            //printf("Skipping improper dihedral: %d %d %d %d\n",iatom,jatom,katom,matom);
#endif
            continue;
        }*/
		jclass = atoms[jatom].classx;
		kclass = atoms[katom].classx;
		mclass = atoms[matom].classx;

		for(i=0;i<ffield->numOfImpropParams;i++){
        //In the CHARMM 19 parameter set, the FIRST atom is the central atom.
#ifdef CHARMM19
          if((iclass==ffield->impropParams[i].atomClass[0] && jclass==ffield->impropParams[i].atomClass[1]
		&& mclass==ffield->impropParams[i].atomClass[2] && kclass==ffield->impropParams[i].atomClass[3])
              && ((jclass!=mclass) || (jatom>matom))){ //prevent duplication, ensure proper stereochemistry on Leu and Val
#elif AMBER
        //in the AMBER set it is the THIRD atom.
          if((iclass==ffield->impropParams[i].atomClass[2] && jclass==ffield->impropParams[i].atomClass[0]
		&& mclass==ffield->impropParams[i].atomClass[1] && kclass==ffield->impropParams[i].atomClass[3])){
                //&& ((kclass!=mclass) || (katom>matom))){ //prevent duplication, ensure proper stereochemistry on Leu and Val
#endif
	     //(iclass==ffield->impropParams[i].atomClass[0] && jclass==ffield->impropParams[i].atomClass[2]
		//&& mclass==ffield->impropParams[i].atomClass[1] && kclass==ffield->impropParams[i].atomClass[3])){
		  //if(jclass==ffield->impropParams[i].atomClass[0] && kclass==ffield->impropParams[i].atomClass[1] && iclass==ffield->impropParams[i].atomClass[2] && mclass==ffield->impropParams[i].atomClass[3]){

		    num = atoms[iatom].numOfImprops;
		    atoms[iatom].impropAtomList[4*num]   = jatom;
		    atoms[iatom].impropAtomList[4*num+1] = katom;
		    atoms[iatom].impropAtomList[4*num+2] = iatom;
		    atoms[iatom].impropAtomList[4*num+3] = matom;
		    atoms[iatom].impropParamType[num]    = i;
		    atoms[iatom].numOfImprops++;

		    num = atoms[jatom].numOfImprops;
		    atoms[jatom].impropAtomList[4*num]   = jatom;
		    atoms[jatom].impropAtomList[4*num+1] = katom;
		    atoms[jatom].impropAtomList[4*num+2] = iatom;
		    atoms[jatom].impropAtomList[4*num+3] = matom;
		    atoms[jatom].impropParamType[num]    = i;
		    atoms[jatom].numOfImprops++;

		    num = atoms[katom].numOfImprops;
		    atoms[katom].impropAtomList[4*num]   = jatom;
		    atoms[katom].impropAtomList[4*num+1] = katom;
		    atoms[katom].impropAtomList[4*num+2] = iatom;
		    atoms[katom].impropAtomList[4*num+3] = matom;
		    atoms[katom].impropParamType[num]    = i;
		    atoms[katom].numOfImprops++;

		    num = atoms[matom].numOfImprops;
		    atoms[matom].impropAtomList[4*num]   = jatom;
		    atoms[matom].impropAtomList[4*num+1] = katom;
		    atoms[matom].impropAtomList[4*num+2] = iatom;
		    atoms[matom].impropAtomList[4*num+3] = matom;
		    atoms[matom].impropParamType[num]    = i;
		    atoms[matom].numOfImprops++;

#ifdef DEBUG_NON_TABULATED
		    printf("improper: iatom: %d, jatom: %d, katom: %d, matom: %d, i: %d\n",iatom,jatom,katom,matom,i);
#endif
		    break;
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
}

/*void topology::total_internal_interaction_energy(forcefield * ffield, double eps, int rdie, double * evdw_internal, double * eelec_internal)
{
    int ifrag;
    double evdw_int1, eelec_int1;
    *evdw_internal = 0.0;
    *eelec_internal = 0.0;
    for (ifrag=0; ifrag<nfrag; ifrag++) {
        fragtypes[frags[ifrag].type]->internal_interaction_energy(ffield,eps,rdie,&evdw_int1,&eelec_int1);
        *evdw_internal += evdw_int1;
        *eelec_internal += eelec_int1;
    }
}*/


//Creates a list of the pairs of atoms that need to be evaluated exactly, without using the nonbond list.
//We need to include (and mark) the 1,2 and 1,3 pairs so they will be counted in the GB sum.
void topology::create_non_tab_list(bool using_cov_tables, std::vector<atom_nb_entry> * atom_nb_list)
{
    int ifrag,jfrag,j,ii,jj,k,l,m,iatom,jatom,temp;
    bool is12, is13, is14;
    atom_nb_entry newentry;
    atom_nb_list->clear();
    //For every pair of fragments, including interactions between atoms in the same fragment.
    /*for (ifrag=0; ifrag<nfrag; ifrag++)
        for (jfrag=ifrag; jfrag<nfrag; jfrag++) {
            //if they are not close and we are not doing an exact simulation, they will be calculated through the tables.
            if (!(closefragments[ifrag*nfrag+jfrag])) continue;
            //for every pair of atoms belonging to the two fragments...
            for (ii=0; ii<fragtypes[frags[ifrag].type]->natom; ii++) {
                for (jj=0; jj<fragtypes[frags[jfrag].type]->natom; jj++) {
		    iatom=frags[ifrag].atoms[ii];
                    jatom=frags[jfrag].atoms[jj];
                    if ((ifrag==jfrag) && (iatom>=jatom)) continue;
                    if (iatom>jatom) {
                        temp=iatom;
                        iatom=jatom;
                        jatom=temp;
                    }*/

    for (iatom=0; iatom<natom; iatom++)
        for (jatom=iatom+1; jatom<natom; jatom++) {
                    //include in list only if both atoms don't belong to the CG region
                    if (!atoms[iatom].is_in_aa_region && !atoms[jatom].is_in_aa_region) continue;
                    //iatom<jatom and check to see if 1-2, 1-3, or 1-4.
                    is12=false;
                    for (k=0; k<atoms[iatom].numOfBondedAtoms; k++)
                        if (atoms[iatom].bondedAtomList[k]==jatom) {
                            is12=true;
                            break;
                        }
                    //if (is12) continue;
                    is13=false;
                    for (k=0; k<atoms[iatom].numOfAngles; k++)
                        if (atoms[iatom].angleAtomList[3*k+2]==jatom) {
                            is13=true;
                            break;
                        }
                    //if (is13) continue;
                    is14=false;
                    for (k=0; k<atoms[iatom].numOfBonded14Atoms; k++)
                        if (atoms[iatom].bonded14AtomList[4*k+3]==jatom) {
                            is14=true;
                            break;
                        }
                    //no 1-2 or 1-3 interactions
                    if (is12 || is13) continue;
                    /*if (using_cov_tables && term_in_covalent_tables(iatom,jatom)) {
#ifdef DEBUG_NON_TABULATED
                        printf("Skipping atom pair: %d %d\n",iatom,jatom);
#endif
                        continue;
                    }*/
#ifdef DEBUG_NON_TABULATED
                    printf("Add atom pair %d %d %c\n",iatom,jatom,yesno(is14));
#endif
                    newentry.iatom=iatom;
                    newentry.jatom=jatom;
                    //newentry.is12_or_13=is12 || is13;
                    newentry.is14=is14;
                    atom_nb_list->push_back(newentry);
                }
        /*    }
        }*/
    printf("Number of atom pairs to be evaluated exactly: %ld\n",atom_nb_list->size());
}


topology::~topology()
{
    int itype;
/*
    for (itype=0; itype<nfragtypes; itype++) delete fragtypes[itype];
    free(fragtypes);
    free(frags);
*/
    free(atoms);
    free(resdef);
    free(resinfo);
    free(segstart);
    free(segend);
    free(first_main_chain_frag);
    free(chaincodes);
    free(iscrot);
    free(jscrot);
    //free(closefragments);
}
