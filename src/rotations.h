#ifndef ROTATIONS_H_INCLUDED
#define ROTATIONS_H_INCLUDED

#include <math.h>

#define FALSE           0
#define TRUE            1
#define TRIG_TABLE_SIZE 8192
#define TINY            1.0/TRIG_TABLE_SIZE
#define TWO_PI          2.0*M_PI
#define ONESQRT2        1.0/sqrt(2.0)
#define MAX_MATRIX      10
#define RAD_TO_DEG                  57.2957795//converts radians to degrees
#define DEG_TO_RAD                  1.74532925e-2//converts degrees to radians

#ifndef NO_TRIG_TABLES
//Fast versions of the trig functions, for quat_to_euler and cart_to_sph
static double atantable[TRIG_TABLE_SIZE+1];
static double acostable[TRIG_TABLE_SIZE+1];


inline double table_atan2(const double y, const double x)
{
        double z,result;
        int mask,flip,index;
        long long int ly, lx;
        //quadrants I, II, IV, III in that order -- designed to work with bit masks
        const double offset[4] = {0.0,M_PI,2.0*M_PI,M_PI};
        const double sign[4] = {1.0,-1.0,-1.0,1.0};
        //const double offset[8] = {0.0,M_PI_2,M_PI,M_PI_2,TWO_PI,1.5*M_PI,M_PI,1.5*M_PI};
        //const double sign[8] = {1.0,-1.0,-1.0,1.0,-1.0,1.0,1.0,-1.0};
        //masks 1, 2, 4, and 7 require negation (an odd number of bits)
        mask=0;
        flip=FALSE;
        if (x<0) mask |= 1;
        if (y<0) mask |= 2;
        z=fabs(y/x); // do not flip
        if (z>1.0) {flip=TRUE; z=1.0/z;}
        index=(int) (z*TRIG_TABLE_SIZE);
        //if (index<TRIG_TABLE_SIZE) {
                result=atantable[index];
        //} else {*/
        //      result=M_PI_4;
        //}
        if (flip) result=M_PI_2-result;
        result = offset[mask] + sign[mask] * result;
        return result;
}

inline double table_acos(const double x)
{
        double result;
        int pos = (x>0);
        int index =(int) (fabs(x)*TRIG_TABLE_SIZE);
        if (index<0) {
                result=M_PI_2;
        } else if (fabs(x)>0.99) {
                result=sqrt(2*(1-fabs(x))); /*First term in series expansion of acos(x) about x=1*/
        } else result=acostable[index];
        if (!pos) result=M_PI-result;
        return result;
}

#else
//replace the table_XXX trig functions by standard implementations
#define table_atan2 atan2
#define table_acos acos
#endif
//Assumes all quaternions have four entries  Differs from wikipedia in that q[0] represents the real part,
//whereas in wikipedia q_4 is the real part.  Steve's code also stores the real part first.
//Convention used is "x-convention 3-1-3" as used in the Wikipedia article "Rotation formalisms in three dimensions."
//The convention is as follows:
//1. Rotate by phi radians about the z axis
//2. Rotate by theta radians about the x'-axis (the x-axis after the first rotation)
//3. Rotate by psi radians about the z''-axis.
//Each rotation is counterclockwise as seen from the positive end of the corresponding axis.

//Convert quaternion to Euler angles.
//we adopt the convention that when theta < TINY, phi - psi = 0, whereas when (pi - theta) < TINY, phi + psi = 0 (or 2*pi).
//Different convention: when we have a gimbal lock, psi is always 0 -- this avoids issues with rounding.
//this requires a normalized quaternion
// 4 multiplciations, 4 additions, 3 inverse trig functions.
//attempting to inline quat_to_euler seems to change the results
/*inline void quat_to_euler(const double * q, double * phi, double * theta, double * psi)
{
    double _phi, _theta, _psi;
    double alpha, beta; // alpha = (phi - psi ) / 2; beta = (phi + psi) / 2;
    //*theta = acos(q[3]*q[3]+q[0]*q[0]-q[1]*q[1]-q[2]*q[2]);
    double costheta, tmp;
    costheta=1.0-2.0*(q[1]*q[1]+q[2]*q[2]);
    //adopt consistent convention for gimbal lock cases
    if ( costheta > (1.0 - TINY) ) {
        _theta = 0.0;
        _phi = 2.0*table_atan2(q[3],q[0]);
        _psi = 0.0;
    } else if ( costheta < (-1.0 + TINY)) {
        _theta = M_PI;
        _phi = 2.0*table_atan2(q[2],q[1]);
        _psi = 0.0;
    } else {
        _theta = table_acos(costheta);
        alpha = table_atan2(q[2],q[1]);
        beta = table_atan2(q[3],q[0]);
        _phi = beta + alpha;
        _psi = beta - alpha;
    };
    //Using the [0,2*pi] range for phi and psi will make indexing the array easier.
    if ( _phi > TWO_PI) _phi = _phi - TWO_PI;
    if ( _phi < 0.0) _phi = _phi + TWO_PI;
    if ( _psi > TWO_PI) _psi = _psi - TWO_PI;
    if ( _psi < 0.0) _psi = _psi + TWO_PI;
    *phi = _phi;
    *theta = _theta;
    *psi = _psi;
}*/

//write q = (w, r)
//v2 = v + 2r x (r x v + wv
//v2 = v + 2(v x r + wv) x r to be consistent with other subroutines
inline void rotate_vector_by_quat(const double * q, const double * v1, double * v2)
{
    double a[3];
    //r=q[1],q[2],q[3] v1=v1[0],v1[1],v1[2]
    a[0]=v1[1]*q[3]-v1[2]*q[2]+q[0]*v1[0];
    a[1]=v1[2]*q[1]-v1[0]*q[3]+q[0]*v1[1];
    a[2]=v1[0]*q[2]-v1[1]*q[1]+q[0]*v1[2];
    v2[0]=v1[0]+2.0*(a[1]*q[3]-q[2]*a[2]);
    v2[1]=v1[1]+2.0*(a[2]*q[1]-q[3]*a[0]);
    v2[2]=v1[2]+2.0*(a[0]*q[2]-q[1]*a[1]);
}


inline void multiply_quat(const double * qa, const double * qb, double * qc)
{
     //take Grossmann product to combine two rotations into one qd=qa*qb
  //rotate through the first vector. qd is our product
  double result[4];
  result[0] = qa[0]*qb[0] - qa[1]*qb[1] - qa[2]*qb[2] - qa[3]*qb[3];
  result[1] = qa[0]*qb[1] + qa[1]*qb[0] + qa[2]*qb[3] - qa[3]*qb[2];
  result[2] = qa[0]*qb[2] - qa[1]*qb[3] + qa[2]*qb[0] + qa[3]*qb[1];
  result[3] = qa[0]*qb[3] + qa[1]*qb[2] - qa[2]*qb[1] + qa[3]*qb[0];
  qc[0]=result[0];
  qc[1]=result[1];
  qc[2]=result[2];
  qc[3]=result[3];
}

//qb=conj(qa)
inline void conjugate_quat(const double * qa, double * qb)
{
    qb[0]=qa[0];
    qb[1]=-qa[1];
    qb[2]=-qa[2];
    qb[3]=-qa[3];
}


//Requires the distance r to be pre-computed.
inline void cart_to_sph(const double * x, const double r, double * sphtheta, double * sphphi)
{
    double costheta,_sphtheta,_sphphi;
    costheta = x[2]/r;
    if (costheta > (1.0-TINY)) {
        _sphtheta = 0.0;
        _sphphi = 0.0;
    } else if (costheta < (-1.0+TINY)) {
        _sphtheta = M_PI;
        _sphphi = 0.0;
    } else {
        _sphtheta = table_acos(costheta);
        _sphphi = table_atan2(x[1],x[0]);
    }
    //not necessary, table_atan2 designed for [0,2pi] range
    //if (_sphphi<0) _sphphi +=TWO_PI;
    //if (_sphphi>=TWO_PI) _sphphi -=TWO_PI;
    *sphtheta=_sphtheta;
    *sphphi=_sphphi;
}

void fill_trig_tables(void);
//void prefetch_trig_tables(void);
void quat_to_euler(const double * q, double * phi, double * theta, double * psi);
void euler_to_quat(double phi, double theta, double psi, double * q);
void quat_to_matrix(double * q, double * r);
//void rotate_vector_by_quat(const double * q, const double * v1, double * v2);
void euler_to_matrix(double phi, double theta, double psi, double * r);
void axisangle_to_quat(double alpha, double * axis, double * q);
void quat_to_axisangle(double * q, double * alpha, double * axis);
void normalize_quat(double * q);
//void multiply_quat(const double * qa, const double * qb, double * qc);
void conjugate_quat(const double * qa, double * qb);
//void divide_quat(double * qa, double * qb, double * qc);
void multiply_conj_quat(double * qa, double * qb, double * qc);
void matmul(double * a, double * b, double * c);
void matmul2(double a[3][3], double b[3][3], double c[3][3]);
void rand_unif_quat(double *q);
void rand_small_quat(double delta, double *q);
//void cart_to_sph(const double * x, const double r, double * sphtheta, double * sphphi);
void sph_to_cart(double r, double sphtheta, double sphphi, double * x);
void rmsd_fit(int natom, double * weight, double * coords1, double * coords2, double * center, double * orient, double * rmsd);
void create_spherical_diffusion_kernel(const int tablesize, const int lmax, const double alpha, double values[], double * work1);
void create_so3_diffusion_kernel(const int tablesize, const int jmax, const double alpha, double values[], double * work1);
void build_atom(double * coords, int i, int j, int k, int l, double bond, double angle, double dihedral);

#endif // ROTATIONS_H_INCLUDED
