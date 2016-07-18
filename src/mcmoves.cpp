#include <cstdio>
#include <vector>
#include <cctype>
//#include "fragments.h"
#include "topology.h"
#include "rotations.h"
#include "mc.h"
#include "ffield.h"
#include "mt.h"
#include "util.h"
//Rotate a selection of fragments about an axis through two atoms.  Does not support PBC.
void simulation::rotate_atoms_by_axis(mc_move * move, const double angle, double * coords)
{
  int iatom,k;
  double axis[3],point[3],axism,quat[4],newquat[4];

  //Compute the axis
  axism=0.0;
  for (k=0; k<3; k++) {
      axis[k]=coords[3*move->jaxis+k]-coords[3*move->iaxis+k];
      axism+=axis[k]*axis[k];
      point[k]=coords[3*move->iaxis+k];
  }
  axism=sqrt(axism);
  for (k=0; k<3; k++) axis[k]/=axism;
  axisangle_to_quat(angle,axis,quat);
  rotate_atoms_by_point(move->movedatoms,&quat[0],&point[0],coords);
}

void simulation::rotate_atoms_by_point(subset atoms, const double * quat, const double * point, double * coords)
{
    double rotmatrix[3][3];
    double disp[3],newdisp[3];
    int iatom,k;
    quat_to_matrix(quat,&rotmatrix[0][0]);
    for (iatom=0; iatom<top->natom; iatom++) if (atoms[iatom]) {
      //Rotate center about axis.
      for (k=0; k<3; k++) disp[k]=coords[3*iatom+k]-point[k];
      matmul(&rotmatrix[0][0],disp,newdisp);
      for (k=0; k<3; k++) coords[3*iatom+k]=point[k]+newdisp[k];
  }
}
//return the subset of all atoms conected to atom iatom through bonds, except those in the exclude subset
subset topology::follow_bonds(int iatom, const subset exclude)
{
    int i, bondedatom, jatom;
    subset atomset, additional_atoms,newexclude;
    atomset.init(natom);
    for (i=0; i<atoms[iatom].numOfBondedAtoms; i++) {
        bondedatom=atoms[iatom].bondedAtomList[i];
        atomset+=bondedatom;
    }
    atomset/=exclude;
    newexclude = exclude;
    newexclude += iatom; // newexclude = exclude U {iatom}
    additional_atoms.init(natom);
    //For each bonded atom not in esclude
    for (jatom=0; jatom<natom; jatom++) if (atomset[jatom]) {
#ifdef DEBUG
        if (atoms[jatom].is_backbone) {
            printf("Warning:  backbone atom %d %s to set\n",atoms[jatom].resNum+1,atoms[jatom].name);
        }
#endif
        additional_atoms|=follow_bonds(jatom,newexclude);
    }
    atomset |= additional_atoms;
    atomset += iatom;
    return atomset;
}

void topology::print_atom_subset(subset atomset) {
    int iatom, ires;
    ires=-1;
    for (iatom=0; iatom<natom; iatom++) if (atomset[iatom]) {
        if (atoms[iatom].resNum!=ires) {
            ires=atoms[iatom].resNum;
            if (iatom!=0) printf("\n");
            printf("%s %d ",atoms[iatom].resName,atoms[iatom].resNum+1);
        }
        printf("%s ",atoms[iatom].name);
    }
    printf("\n");
}

//We cannot use "follow_bonds" to generate backbone or backrub moves
//because they will not include side chains in the coarse grained region
void topology::generate_backbone_moves(vector<mc_move> * backbone_moves)
{
    int ires,iatom;
    mc_move phi_move,psi_move;
    bool is_ace, is_nme;
    //subset phiatoms, psiatoms;
    backbone_moves->clear();
    for (ires=0; ires<nres; ires++) if (ires!=ligand_res) {
        is_ace=(strcmp(resdef[resinfo[ires].restype].name,"ACE")==0);
        is_nme=(strcmp(resdef[resinfo[ires].restype].name,"NME")==0);
        //need to watch the pointers and copying of objects
        //phiatoms.init(natom);
        phi_move.movedatoms.init(natom);
        phi_move.iaxis=-1;
        phi_move.jaxis=-1;
        //psiatoms.init(natom);
        psi_move.movedatoms.init(natom);
        psi_move.iaxis=-1;
        psi_move.jaxis=-1;
        for (iatom=0; iatom<natom; iatom++) if (resinfo[atoms[iatom].resNum].whichseg==resinfo[ires].whichseg) {
            //atoms to be moved for phi move: all atoms in residues N-terminal of the one to be moved, plus N and H from the current residue
            //for psi move: C and O from the current residue, plus all atoms in residues C-terminal of it
            if (atoms[iatom].resNum<ires) {
                phi_move.movedatoms+=iatom;
            } else if (atoms[iatom].resNum>ires) {
                psi_move.movedatoms+=iatom;
            } else if (strcmp(atoms[iatom].name,"N")==0) {
                phi_move.iaxis=iatom;
                phi_move.movedatoms+=iatom;
            } else if ((atoms[iatom].name[0]=='H') && (atoms[iatom].name[1]=='\0' || (!is_nme && isdigit(atoms[iatom].name[1])))) {
                //to get H, H1, H2, or H3 (also H4, etc. but these are not normal)  -- don't get H1/H2/H3 for NME
                phi_move.movedatoms+=iatom;
            } else if ((strcmp(atoms[iatom].name,"CA")==0) || (strcmp(atoms[iatom].name,"CH3")==0)) {
                //does not need to be moved, but record index
                phi_move.jaxis=iatom;
                psi_move.iaxis=iatom;
            } else if (strcmp(atoms[iatom].name,"C")==0) {
                psi_move.movedatoms+=iatom;
                psi_move.jaxis=iatom;
            } else if ((strcmp(atoms[iatom].name,"O")==0) || (strcmp(atoms[iatom].name,"OXT")==0)) {
                psi_move.movedatoms+=iatom;
            } //all other atoms don't need to be processed
        }
	//ACE has no phi move possible, and NME has no psi move possible.
        if (!is_ace) backbone_moves->push_back(phi_move);
        if (!is_nme) backbone_moves->push_back(psi_move);
#ifdef DEBUG
        printf("Backbone phi move for residue %d %s:\n",ires+1,resdef[resinfo[ires].restype].name);
        print_atom_subset(phi_move.movedatoms);
        printf("Backbone psi move for residue %d %s:\n",ires+1,resdef[resinfo[ires].restype].name);
        print_atom_subset(psi_move.movedatoms);
#endif
    }
}

void topology::generate_backrub_moves(vector<mc_move> * backrub_moves)
{
    int ires, jres, iatom;
    mc_move move;
    backrub_moves->clear();
    for (ires=0; ires<nres; ires++)
        for (jres=ires+1; jres<nres; jres++)
            //both start and end residues must belong to same chain, and we can't end on a proline
            if ((resinfo[ires].whichseg==resinfo[jres].whichseg) && (strcmp(resdef[resinfo[jres].restype].name,"PRO")!=0)) {
                move.movedatoms.init(natom);
                move.iaxis=resinfo[ires].branchatom;
                move.jaxis=resinfo[jres].branchatom;
                for (iatom=0; iatom<natom; iatom++) {
                    //want C,O from ires, N, H, from jres, all atoms in between
                    if (((atoms[iatom].resNum>ires) && (atoms[iatom].resNum<jres)) ||
                        ((atoms[iatom].resNum==ires) && ((strcmp(atoms[iatom].name,"C")==0) || (strcmp(atoms[iatom].name,"O")==0))) ||
                        ((atoms[iatom].resNum==jres) && ((strcmp(atoms[iatom].name,"N")==0) || (strcmp(atoms[iatom].name,"H")==0)))) {
                            move.movedatoms+=iatom;
                        }
                }
                backrub_moves->push_back(move);
#ifdef DEBUG
                    printf("Adding backrub move %d %d:\n",ires,jres);
                    print_atom_subset(move.movedatoms);
#endif
            }
}


void topology::generate_sidechain_moves(vector<mc_move> * sidechain_moves)
{
    int ibond,iatom,nligand,n;
    mc_move newmove;
    subset exclude;
    subset temp;
    sidechain_moves->clear();
    exclude.init(natom);
    temp.init(natom);
    nligand=0;
    for (iatom=0; iatom<natom; iatom++) if (ligand[iatom]) nligand++;
    //We use the rotatable bond information collected while inserting the residues (topology::insert_residue)
    for (ibond=0; ibond<nscrot; ibond++) {
        newmove.iaxis=iscrot[ibond];
        newmove.jaxis=jscrot[ibond];
        exclude.init(natom);
        exclude+=newmove.iaxis;
        //need to prevent walking up the backbone for CB-CD "bond" in proline
        if (strcmp(atoms[newmove.iaxis].resName,"PRO")==0) {
            exclude+=find_atom(atoms[newmove.iaxis].resNum,"N");
        }
        newmove.movedatoms=follow_bonds(newmove.jaxis,exclude);
        //If this move involves the ligand, then check to see if it involves more than half the atoms in the ligand.
        //If so, invert it with respect to the ligand, so that the majority of atoms will not be moved
        if (ligand[newmove.iaxis] && ligand[newmove.jaxis]) {
            n=0;
            for (iatom=0; iatom<natom; iatom++) if (newmove.movedatoms[iatom]) n++;
            if (n>(nligand/2)) {
                temp=ligand;
                temp/=newmove.movedatoms;
                newmove.movedatoms=temp;
            }
        }
        sidechain_moves->push_back(newmove);
#ifdef DEBUG
        printf("Adding sidechain move:\n");
        print_atom_subset(newmove.movedatoms);
#endif
    }
}

void simulation::do_ligand_trans(double movesize_small, double large_dist_frac, double movesize_large, double * coords)
{
    double disp[3],m,r;
    int k,iatom;
    //construct a random vector within the unit sphere, and multiply by movesize
    //to give a random displacement whose magnitude is no more than movesize.
    r=genrand_real3();
    if (r<large_dist_frac) {
        rand_trans_vector(movesize_large,&disp[0]);
    } else {
        rand_trans_vector(movesize_small,&disp[0]);
    }
    //translate all ligand atoms
    for (iatom=0; iatom<top->natom; iatom++) if (top->ligand[iatom]) {
        for (k=0; k<3; k++) coords[3*iatom+k]+=disp[k];
    }
}

void simulation::do_ligand_rot(double movesize_small, double large_dist_frac, double movesize_large, double * coords)
{
    double quat[4],com[3],mass,totmass,r;
    int k, iatom;
    //find the center of mass of the ligand.
    totmass=0.0;
    for (k=0; k<3; k++) com[k]=0.0;
    for (iatom=0; iatom<top->natom; iatom++) if (top->ligand[iatom]) {
        mass=top->atoms[iatom].mass;
        for (k=0; k<3; k++) com[k]+=mass*coords[3*iatom+k];
        totmass+=mass;
    }
    for (k=0; k<3; k++) com[k]/=totmass;
    r=genrand_real3();
    if (r<large_dist_frac) {
        rand_small_quat(movesize_large,&quat[0]);
    } else {
        rand_small_quat(movesize_small,&quat[0]);
    }
    rotate_atoms_by_point(top->ligand,&quat[0],&com[0],coords);
}

//randomly move a heavy atom, togehter with any hydrogen atoms that may be bonded to it
void simulation::heavy_atom_trans(subset * movedatoms, double movesize, double * coords)
{
    int iatom,iheavy,j,jatom,k;
    double disp[3];
    //choose a heavy atom at random
    for (;;) {
        iheavy=int(genrand_real3()*top->natom);
        if (top->atoms[iheavy].is_in_aa_region && top->atoms[iheavy].atomicNum>1) break;
    }
    *movedatoms+=iheavy;
    for (j=0; j<top->atoms[iheavy].numOfBondedAtoms; j++) {
        jatom=top->atoms[iheavy].bondedAtomList[j];
        if (top->atoms[jatom].atomicNum==1) *movedatoms+=jatom;
    }
    rand_trans_vector(movesize,&disp[0]);
    for (iatom=0; iatom<top->natom; iatom++) if ((*movedatoms)[iatom]) {
        for (k=0; k<3; k++) coords[3*iatom+k]+=disp[k];
    }
}

void simulation::heavy_atom_rot(subset * movedatoms, double movesize, double * coords) 
{
    int iatom,iheavy,j,jatom,k;
    double q[4];
    //choose a heavy atom at random
    for (;;) {
        iheavy=int(genrand_real3()*top->natom);
        if (top->atoms[iheavy].is_in_aa_region && top->atoms[iheavy].atomicNum>1) break;
    }
    *movedatoms+=iheavy;
    for (j=0; j<top->atoms[iheavy].numOfBondedAtoms; j++) {
        jatom=top->atoms[iheavy].bondedAtomList[j];
        if (top->atoms[jatom].atomicNum==1) *movedatoms+=jatom;
    }
    rand_small_quat(movesize,&q[0]);
    rotate_atoms_by_point(*movedatoms,&q[0],&coords[3*iheavy],coords);
}

void simulation::mcmove(int * movetype, subset * movedatoms, double * coords)
{
    mc_move * movelist;
    int moveindex,move, movecount;
    double angle,r;
    movedatoms->init(top->natom);

    r=genrand_real3();
    for (move=1; move<=NUM_MOVES; move++) {
        if (r<=cumprob[move]) break;
    }
    *movetype=move;
    //choose one of the lists from which to select the move; getting a pointer to the vector object itself
    //results in memory errors
    switch (move) {
        case MOVE_BACKBONE:
            movelist=&backbone_moves[0];
            movecount=backbone_moves.size();
            break;
        case MOVE_SIDECHAIN:
            movelist=&sidechain_moves[0];
            movecount=sidechain_moves.size();
            break;
        case MOVE_BACKRUB:
            movelist=&backrub_moves[0];
            movecount=backrub_moves.size();
            break;
        case MOVE_LIGAND_TRANS:
            do_ligand_trans(movesize[MOVE_LIGAND_TRANS],large_dist_frac[MOVE_LIGAND_TRANS],movesize_large[MOVE_LIGAND_TRANS],coords);
            *movedatoms=top->ligand;
            break;
        case MOVE_LIGAND_ROT:
            do_ligand_rot(movesize[MOVE_LIGAND_ROT],large_dist_frac[MOVE_LIGAND_ROT],movesize_large[MOVE_LIGAND_ROT],coords);
            *movedatoms=top->ligand;
            break;
        case MOVE_HEAVY_TRANS:
            heavy_atom_trans(movedatoms,movesize[MOVE_HEAVY_TRANS],coords);
            break;
        case MOVE_HEAVY_ROT:
            heavy_atom_rot(movedatoms,movesize[MOVE_HEAVY_ROT],coords);
            break;
        default:
            printf("Error in switch statement.\n");
            die();
    }
    if ((move==MOVE_BACKBONE) || (move==MOVE_SIDECHAIN) || (move==MOVE_BACKRUB)) {
        //pick a random move from the list
        moveindex=(int) (genrand_real3()*movecount);
        *movedatoms=movelist[moveindex].movedatoms;
        //top->print_atom_subset(*movedatoms);
        angle=(2.0*genrand_real3()-1.0)*movesize[move];
        rotate_atoms_by_axis(&movelist[moveindex],angle,coords);
    }
}
