#include <cstdio>
//#include "fragments.h"
#include "topology.h"
#include "rotations.h"
#include "mc.h"
#include "ffield.h"
#include "mt.h"
#include "util.h"
//Rotate a selection of fragments about an axis through two atoms.  Does not support PBC.
void simulation::rotate_fragments_by_axis(const bool * moved, const int atom1, const int atom2, const double angle, double * coords)
{
  int iatom,k;
  double axis[3],point[3],axism,quat[4],newquat[4];
  double rotmatrix[3][3];
  double disp[3],newdisp[3];
  //Compute the axis
  axism=0.0;
  for (k=0; k<3; k++) {
      axis[k]=oldcoords[3*atom2+k]-oldcoords[3*atom1+k];
      axism+=axis[k]*axis[k];
      point[k]=oldcoords[3*atom1+k];
  }
  axism=sqrt(axism);
  for (k=0; k<3; k++) axis[k]/=axism;
  axisangle_to_quat(angle,axis,quat);
  quat_to_matrix(quat,&rotmatrix[0][0]);
  for (iatom=0; iatom<top->natom; iatom++) if (moved[iatom]) {
      //Rotate center about axis.
      for (k=0; k<3; k++) disp[k]=coords[3*iatom+k]-point[k];
      matmul(&rotmatrix[0][0],disp,newdisp);
      for (k=0; k<3; k++) coords[3*iatom+k]=point[k]+newdisp[k];
  }
}



