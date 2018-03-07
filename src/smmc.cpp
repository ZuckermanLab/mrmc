#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
//#include "tables.h"
#include "mc.h"
#include "mt.h"
#include "rotations.h"

//Reads structures from PDB file, generate list of transforamtions,
void simulation::generate_smmc_trans(char * pdbfname, char * outfname)
{
    FILE * input;
    FILE * output;
    subset valid_coords;
    int ipose, jpose,index, iatom, k, idih,ndih;
    double backbone_rmsd,ligand_rmsd,dih1,dih2;
    smmc_nposes = 0;
    smmc_structures=NULL;
    input = fopen(pdbfname,"r");
    printf("Reading poses from file %s.\n",pdbfname);
    for (;;) {
        smmc_structures=(double *) checkrealloc(smmc_structures,((smmc_nposes+1)*3*top->natom),sizeof(double));
        //not all atoms may be read -- we may use only the backbone atoms
        for (iatom=0; iatom<top->natom; iatom++) for (k=0; k<3; k++) smmc_structures[3*top->natom*smmc_nposes+3*iatom+k]=0.0;
        top->read_pdb_stream(input,&smmc_structures[3*top->natom*smmc_nposes],valid_coords);
    	if (feof(input)) break;
        smmc_nposes++;
    }
    fclose(input);
    if (smmc_nposes<2) {
        printf("Error: you must have at least 2 poses in order to use SMMC.  You only have %d.\n",smmc_nposes);
        die();
    } else {
        printf("Total of %d SMMC poses read.\n",smmc_nposes);
    }
    smmc_trans = (smmc_trans_info *) checkalloc(smmc_nposes*smmc_nposes,sizeof(smmc_trans_info));
    output = fopen(outfname,"w");
    fprintf(output,"%d\n",smmc_nposes);
    for (index=0; index<smmc_nposes*smmc_nposes; index++) smmc_trans[index].dih_diff=NULL;
    ndih=sidechain_moves.size();
    for (ipose=0; ipose<smmc_nposes; ipose++)
        for (jpose=0; jpose<smmc_nposes; jpose++) if (ipose!=jpose) {
            index=ipose*smmc_nposes+jpose;
            printf("Aligning poses %d %d\n",ipose,jpose);
            get_aligned_ligand_rmsd(&smmc_structures[3*top->natom*ipose],&smmc_structures[3*top->natom*jpose],&backbone_rmsd,
                    &smmc_trans[index].disp[0],&smmc_trans[index].rot[0],&ligand_rmsd);
            fprintf(output,"%d %d %.10f %.10f %.10f %.10f %.10f %.10f %.10f\n",ipose+1,jpose+1,
                    smmc_trans[index].disp[0],smmc_trans[index].disp[1],smmc_trans[index].disp[2],
                    smmc_trans[index].rot[0],smmc_trans[index].rot[1],smmc_trans[index].rot[2],smmc_trans[index].rot[3]);
            smmc_trans[index].dih_diff=(double *) checkalloc(ndih,sizeof(double));
            for (idih=0; idih<ndih; idih++) {
                dih1=dihedral_angle(&smmc_structures[3*top->natom*ipose],sidechain_moves[idih].iaxis_prev,sidechain_moves[idih].iaxis,
                                    sidechain_moves[idih].jaxis,sidechain_moves[idih].jaxis_succ);
                dih2=dihedral_angle(&smmc_structures[3*top->natom*jpose],sidechain_moves[idih].iaxis_prev,sidechain_moves[idih].iaxis,
                                    sidechain_moves[idih].jaxis,sidechain_moves[idih].jaxis_succ);
                smmc_trans[index].dih_diff[idih]=dih2-dih1;
            }
    }
    fclose(output);
}



void simulation::perform_smmc_move(subset * movedatoms, double * coords)
{
    int ipose, ipose_closest, jpose, index,iatom,k,idih;
    double * wt_backbone;
    double disp1[3],q1[4],matrix[3][3], ligand_com[3],origin[3];
    double backbone_disp[3],backbone_q[4],backbone_q_conj[4],backbone_rmsd,ligand_rmsd,closest_ligand_rmsd;
#ifdef DEBUG
    FILE * output;
    output=fopen("smmc-test.pdb","w");
    fprintf(output,"MODEL %6d\n",1);
    top->write_pdb_stream(output,coords,NULL);
    fprintf(output,"ENDMDL\n");
#endif
    //use voronoi procedure to choose starting pose
    closest_ligand_rmsd=10000;
    for (ipose=0; ipose<smmc_nposes; ipose++) {
        get_aligned_ligand_rmsd(&smmc_structures[3*top->natom*ipose],coords,&backbone_rmsd,
            &disp1[0],&q1[0],&ligand_rmsd);
        if (ligand_rmsd<closest_ligand_rmsd) {
            closest_ligand_rmsd=ligand_rmsd;
            ipose_closest=ipose;
        }
    }
    //choose a random transformation from smmc_trans
    do {
        jpose=int(smmc_nposes*genrand_real3());
    } while (ipose_closest==jpose);
#ifdef DEBUG
    printf("Shift move pose %d %d\n",ipose_closest,jpose);
#endif
    index=ipose_closest*smmc_nposes+jpose;
    //the transformation is in the frame of reference of structure "ipose".
    //We need to transform it to the frame of reference of hte current structure.
    //This requires we perform a backbone alignment on the current set of coordinates and
    //double backbone_disp[3],backbone_q[4], backbone_q2[4],backbone_center[3], nbackbone,m[3][3],temp[3]; //_total;
    wt_backbone=(double *) checkalloc(top->natom,sizeof(double));
    for (iatom=0; iatom<top->natom; iatom++) {
        if (strstr("N CA C",top->atoms[iatom].name)!=NULL) wt_backbone[iatom]=1.0; else wt_backbone[iatom]=0.0;
    }
    rmsd_fit(top->natom,wt_backbone,&smmc_structures[3*top->natom*ipose_closest],coords,&backbone_disp[0],&backbone_q[0],&backbone_rmsd);
    //backbone_disp and backbone_q are the transformation needed to take structure ipose to the current coordinates.
    //smmc_trans[index] contains the transformation needed in the frame of reference of structure ipose.
    conjugate_quat(&backbone_q[0],&backbone_q_conj[0]);
    //we need to construct a quaternion representing a rotation of the same size as smmc_trans.rot,
    //but has been rotated through backbone_q_conj
    //the quaternion algebra version of this is (backbone_q_conj) * (smmc_trans.rot) * (backbone_q)
    //but we can also do this by calling rotate_vector_by_quat on the imaginary components of smmc_trans.rot
    q1[0] = smmc_trans[index].rot[0];
    rotate_vector_by_quat(&backbone_q_conj[0],&smmc_trans[index].rot[1],&q1[1]);
    //we do the same for the translational part of smmc_trans
    rotate_vector_by_quat(&backbone_q_conj[0],&smmc_trans[index].disp[0],&disp1[0]);
    //disp1 and q1 are the transformation we need, in our current frame of reference.
    //rotate the ligand through q1 about center of mass
    get_ligand_com(coords,&ligand_com[0]);
    rotate_atoms_by_point(top->ligand,&smmc_trans[index].rot[0],&ligand_com[0],coords);
    //translate by disp1
    for (iatom=0; iatom<top->natom; iatom++) if (top->ligand[iatom]) {
        for (k=0; k<3; k++) coords[3*iatom+k]+=disp1[k];
    }
    *movedatoms=top->ligand;
    //do sidechain rotations
    if (smmc_sidechain_rotations) {
        for (idih=0; idih<sidechain_moves.size(); idih++) {
            //the minus sign here corrects for a difference in the sense of rotate_axis_by_angle
            //compared to dihedral_angle
            rotate_atoms_by_axis(sidechain_moves[idih],-smmc_trans[index].dih_diff[idih],coords);
            (*movedatoms)|=sidechain_moves[idih].movedatoms;
        }
    }
    //we're done
    free(wt_backbone);
#ifdef DEBUG
    fprintf(output,"MODEL %6d\n",2);
    top->write_pdb_stream(output,coords,NULL);
    fprintf(output,"ENDMDL\n");
    fclose(output);
#endif

}

