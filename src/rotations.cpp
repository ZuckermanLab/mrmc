#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "rotations.h"
#include "mt.h"
#include "util.h"

#define FALSE           0
#define TRUE            1
#define TRIG_TABLE_SIZE 8192
#define TINY            1.0/TRIG_TABLE_SIZE
#define TWO_PI          2.0*M_PI
#define ONESQRT2        1.0/sqrt(2.0)
#define MAX_MATRIX      10

//Fast versions of the trig functions, for quat_to_euler and cart_to_sph
double atantable[TRIG_TABLE_SIZE+1];
double acostable[TRIG_TABLE_SIZE+1];


double table_atan2(const double y, const double x)
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
        if ((x==0) && (y>=0)) return M_PI_2;
        if ((x==0) && (y<0)) return 1.5*M_PI;
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

double table_acos(const double x)
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

void fill_trig_tables(void)
{
        int i;
        double x;
        for (i=0; i<TRIG_TABLE_SIZE; i++) {
                x=((double) i + 0.5)/TRIG_TABLE_SIZE;
                atantable[i]=atan(x);
                acostable[i]=acos(x);
        }
	atantable[TRIG_TABLE_SIZE]=M_PI_4; //this way, if we overrun, we get pi/4
	acostable[TRIG_TABLE_SIZE]=0;
}


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
void quat_to_euler(const double * q, double * phi, double * theta, double * psi)
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
}


//convert Euler angles to quaternion; makes sure positive real part
//6 trig functions, 2 additions/subtractions, 3 multiplications by half, 4 arbitrary multiplications, 4 negations in half the cases
//the quaternion produced should be normalized based on the mathematical properties of trig functions, but normalization is not explicitly done
void euler_to_quat(double phi, double theta, double psi, double * q)
{
    double alpha = 0.5 * (phi - psi);
    double beta = 0.5 * (phi + psi);
    double chtheta = cos(0.5 * theta);
    double shtheta = sin(0.5 * theta);
    q[1] = cos(alpha) * shtheta;
    q[2] = sin(alpha) * shtheta;
    q[3] = sin(beta) * chtheta;
    q[0] = cos(beta) * chtheta;
    if (q[0]<0) {
        q[1] = -q[1];
        q[2] = -q[2];
        q[3] = -q[3];
        q[0] = -q[0];
    }
}

//r is indexed 11, 12, 13, 21, 22, 23, 31, 32, 33
//assumes that r is in fact a 3x3 array, we treat it as single dimensional
//assumes normalized quaternion

//from wikipedia article, first formula (with q[0] as real part, not q_4)
//It is the same as Steve's LBMC code.
void quat_to_matrix(double * q, double * r)
{
    r[0]=1.0-2.0*(q[2]*q[2]+q[3]*q[3]);
    r[1]=2.0*(q[1]*q[2]-q[3]*q[0]);
    r[2]=2.0*(q[1]*q[3]+q[2]*q[0]);
    r[3]=2.0*(q[1]*q[2]+q[3]*q[0]);
    r[4]=1.0-2.0*(q[1]*q[1]+q[3]*q[3]);
    r[5]=2.0*(q[2]*q[3]-q[1]*q[0]);
    r[6]=2.0*(q[1]*q[3]-q[2]*q[0]);
    r[7]=2.0*(q[2]*q[3]+q[1]*q[0]);
    r[8]=1.0-2.0*(q[1]*q[1]+q[2]*q[2]);
}

//Needed for processing crystallographic transformations.
//Used to be bottom end of rmsd_fit (eq. 10 from paper by Coutsias et al.)
void matrix_to_quat(const double r[3][3], double * q, double * maxev)
{
    double f[4][4], qq[4][4], q2[4],ev;
    int iev,k;
    f[0][0]=r[0][0]+r[1][1]+r[2][2];
    f[0][1]=r[1][2]-r[2][1];
    f[0][2]=r[2][0]-r[0][2];
    f[0][3]=r[0][1]-r[1][0];
    f[1][1]=r[0][0]-r[1][1]-r[2][2];
    f[1][2]=r[0][1]+r[1][0];
    f[1][3]=r[0][2]+r[2][0];
    f[2][2]=-r[0][0]+r[1][1]-r[2][2];
    f[2][3]=r[1][2]+r[2][1];
    f[3][3]=-r[0][0]-r[1][1]+r[2][2];
    //The matrix is symmetric.
    f[1][0]=f[0][1];
    f[2][0]=f[0][2];
    f[3][0]=f[0][3];
    f[2][1]=f[1][2];
    f[3][1]=f[1][3];
    f[3][2]=f[2][3];
    //Calculate the eigenvalues with Jacobi's algorithm.
    jacobi(4,&f[0][0],&qq[0][0]);
    //Find the maximum eigenvalue and which one it is.
    ev=-1.0e20;
    for (k=0; k<4; k++) if (f[k][k]>ev) {
        ev=f[k][k];
        iev=k;
    }
    //The eigenvectors are stored as columns of qq. (Check this.)
    //The quaternion is the eigenvector corresponding to the largest eigenvalue.
    for (k=0; k<4; k++) q2[k]=qq[k][iev];
    normalize_quat(q2); //Just to make sure.  Among other things ensures positive real part of quaternion.
    conjugate_quat(q2,q);//For some reason we need this.  It probably compensates for having swapped something earlier.
    *maxev=ev;
}

//overload not to have to provide for the maximum eigenvalue
void matrix_to_quat(const double r[3][3], double * q)
{
    double maxev;
    matrix_to_quat(r,q,&maxev);
}
//write q = (w, r)
//v2 = v + 2r x (r x v + wv
//v2 = v + 2(v x r + wv) x r to be consistent with other subroutines
void rotate_vector_by_quat(const double * q, const double * v1, double * v2)
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
//formulas from MathWorld article "Euler Angles", first part. They define their rotations
//in the opposite sense to the Wikipedia article, so the matrix has been transposed to compensate.
//possibly less efficient; for clarity, directness, testing
//more efficient to go via quaternions
void euler_to_matrix(double phi, double theta, double psi, double * r)
{
    double cphi,sphi,ctheta,stheta,cpsi,spsi;
    cphi=cos(phi);
    sphi=sin(phi);
    ctheta=cos(theta);
    stheta=sin(theta);
    cpsi=cos(psi);
    spsi=sin(psi);
    r[0]=cpsi*cphi-ctheta*sphi*spsi;
    r[1]=-spsi*cphi-ctheta*sphi*cpsi;
    r[2]=stheta*sphi;
    r[3]=cpsi*sphi+ctheta*cphi*spsi;
    r[4]=-spsi*sphi+ctheta*cphi*cpsi;
    r[5]=-stheta*cphi;
    r[6]=spsi*stheta;
    r[7]=cpsi*stheta;
    r[8]=ctheta;
}

void axisangle_to_quat(const double alpha, const double * axis, double * q)
{
    double chalpha,shalpha;
    chalpha=cos(0.5*alpha);
    shalpha=sin(0.5*alpha);
    q[0]=chalpha;
    q[1]=shalpha*axis[0];
    q[2]=shalpha*axis[1];
    q[3]=shalpha*axis[2];
}

void quat_to_axisangle(double * q, double * alpha, double * axis)
{
    double s;
    *alpha = 2.0*acos(q[0]); // assumes positive q[0]!
    s = 1/sqrt(q[1]*q[1]+q[2]*q[2]+q[3]*q[3]); // = sin(alpha/2)
    axis[0]=q[1]*s;
    axis[1]=q[2]*s;
    axis[2]=q[3]*s;
}
//qc=qa*qb
void normalize_quat(double * q)
{
    double qnorm;
    qnorm=1/sqrt(q[0]*q[0]+q[1]*q[1]+q[2]*q[2]+q[3]*q[3]);
    if (!is_quat_ok(q)) qnorm=-qnorm; //ensure positive real part
    q[0]=q[0]*qnorm;
    q[1]=q[1]*qnorm;
    q[2]=q[2]*qnorm;
    q[3]=q[3]*qnorm;
}

void multiply_quat(const double * qa, const double * qb, double * qc)
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

//The distance between two quaternions, cos(alpha/2) = Re(qa conj(qb))

//need "fabs" here in case quaternions turn out to be negatives of each other
double dist(const double * qa, const double * qb)
{
    double result = fabs(qa[0]*qb[0] + qa[1]*qb[1] + qa[2]*qb[2] + qa[3]*qb[3]);
    return result;
}

bool is_quat_ok(const double * q)
{
    int k;
    k=0;
    //possibly change this to skip tiny values?
    while ((k<4) && (fabs(q[k])<1e-6)) k++;
    if (k==4) return false; //zero quat.
    return (q[k]>0);
}

//qc=qa*qx^{-1} = qa * conj(qx) / |qx|^2
//qc=qa*conj(qb) = qa * qb^-1 if qb is a unit quaternion
/*void multiply_conj_quat(double * qa, double * qb, double * qc)
{
    //double qnorm;
    //qnorm=qb[0]*qb[0]+qb[1]*qb[1]+qb[2]*qb[2]+qb[3]*qb[3];
  //rotate through the first vector. qd is our product
    qc[0] = ( qa[0]*qb[0] + qa[1]*qb[1] + qa[2]*qb[2] + qa[3]*qb[3]);
    qc[1] = (-qa[0]*qb[1] + qa[1]*qb[0] - qa[2]*qb[3] + qa[3]*qb[2]);
    qc[2] = (-qa[0]*qb[2] + qa[1]*qb[3] + qa[2]*qb[0] - qa[3]*qb[1]);
    qc[3] = (-qa[0]*qb[3] - qa[1]*qb[2] + qa[2]*qb[1] + qa[3]*qb[0]);
  /*qc[0] = (qx[0]*qa[0] + qb[1]*qa[1] + qb[2]*qa[2] + qb[3]*qa[3])/qnorm;
  qc[1] = (qb[0]*qa[1] - qb[1]*qa[0] - qb[2]*qa[3] - qb[3]*qa[2])/qnorm;
  qc[2] = (qb[0]*qa[2] + qb[1]*qa[3] - qb[2]*qa[0] - qb[3]*qa[1])/qnorm;
  qc[3] = (qb[0]*qa[3] - qb[1]*qa[2] + qb[2]*qa[1] - qb[3]*qa[0])/qnorm;*/
//}

//qb=conj(qa)
void conjugate_quat(const double * qa, double * qb)
{
    qb[0]=qa[0];
    qb[1]=-qa[1];
    qb[2]=-qa[2];
    qb[3]=-qa[3];
}
//c=A.b, A a 3x3 matrix and b a vector
//Transposed A to be consistent with quat_to_matrix.
void matmul(double * a, double * b, double * c)
{
    c[0]=a[0]*b[0]+a[3]*b[1]+a[6]*b[2];
    c[1]=a[1]*b[0]+a[4]*b[1]+a[7]*b[2];
    c[2]=a[2]*b[0]+a[5]*b[1]+a[8]*b[2];
}


//the following two routines are adapted from Steve's code

void rand_unif_quat(double *q){
  //generates a random quaternion by first generating random angles for rotation
  double s;
  double sigma1,sigma2;
  double theta1,theta2;

  s = genrand_real3();
  sigma1 = sqrt(1-s);
  sigma2 = sqrt(s);
  theta1 = 2*M_PI*genrand_real3();
  theta2 = 2*M_PI*genrand_real3();

  q[0] = cos(theta2)*sigma2;
  q[1] = sin(theta1)*sigma1;
  q[2] = cos(theta1)*sigma1;
  q[3] = sin(theta2)*sigma2;
  if (q[0]<0) {
      q[0]=-q[0];
      q[1]=-q[1];
      q[2]=-q[2];
      q[3]=-q[3];
  }
}

//the maximum rotation produced by this is delta * pi radians, or delta * 180 degrees
void rand_small_quat(double delta, double *q){
  double dot;
  double alpha,gamma;
  double r,s;
  double qa[4],qb[4];

  //set up quaternion corresponding to identity matrix
  qa[0] = 1.0;
  qa[1] = 0;
  qa[2] = 0;
  qa[3] = 0;

  //generate random uniform quaternion
  rand_unif_quat(&qb[0]);


  //the next few sections are merely to ensure that quaternion rotations always have positive dot products - we don't want to flip; we want to rotate.
  //take dot product between two quaternions
  //dot = qa[0]*qb[0] + qa[1]*qb[1] + qa[2]*qb[2] + qa[3]*qb[3];

  //wrap around qb for angles larger than 90 deg because q=-q
  //corrects if a reflection is occuring
  /*f(qb[0]<0){
    qb[0] = -qb[0];
    qb[1] = -qb[1];
    qb[2] = -qb[2];
    qb[3] = -qb[3];
    //dot = -dot;
  }*/


  //This uses something called "SLERP" to control how large an angle of rotation can occur
  //the 'delta' acts as a cap much the same way that a hard angle cap would.
  //Delta should be though of a unitarion i.e. its the same as a radian except measured out of 1 instead of out of 2*pi since 2*pi doesn't really make sense in 3-D
  //use "slerp" algorithm (spherical interpolation) between qa and qb
  alpha = acos(qb[0]);
  gamma = 1/sin(alpha);
  r = sin((1-delta)*alpha)*gamma;
  s = sin(delta*alpha)*gamma;

  //reuse qa as spherical interpolation between qa and qb
  //this finishes the slerp interpolation
  qa[0] = r*qa[0] + s*qb[0];
  qa[1] = r*qa[1] + s*qb[1];
  qa[2] = r*qa[2] + s*qb[2];
  qa[3] = r*qa[3] + s*qb[3];
  normalize_quat(qa);
  //normalize qa and write to q
  //dot = qa[0]*qa[0] + qa[1]*qa[1] + qa[2]*qa[2] + qa[3]*qa[3];
  q[0] = qa[0];
  q[1] = qa[1];
  q[2] = qa[2];
  q[3] = qa[3];
  }

void matmul2(double a[3][3], double b[3][3], double c[3][3])
{
    int i,j,k;
    for (i=0; i<3; i++)
        for (j=0; j<3; j++) {
                c[i][j]=0.0;
                for (k=0; k<3; k++) c[i][j]+=a[i][k]*b[k][j];
        }
}

/*Requires the distance r to be pre-computed.*/
void cart_to_sph(const double * x, const double r, double * sphtheta, double * sphphi)
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
        if (fabs(x[0])>TINY) {
        	_sphphi = table_atan2(x[1],x[0]);
	} else {
		if (x[1]<0) _sphphi=3*M_PI/2;
		if (x[1]>0) _sphphi=M_PI/2;
	}
    }
    //not necessary, table_atan2 designed for [0,2pi] range
    //if (_sphphi<0) _sphphi +=TWO_PI;
    //if (_sphphi>=TWO_PI) _sphphi -=TWO_PI;
    *sphtheta=_sphtheta;
    *sphphi=_sphphi;
}

void sph_to_cart(double r, double sphtheta, double sphphi, double * x)
{
    x[0] = r * sin(sphtheta) * cos(sphphi);
    x[1] = r * sin(sphtheta) * sin(sphphi);
    x[2] = r * cos(sphtheta);
}

//This calculates the spherical diffusion kernel, sum exp(-l(l+1)*alpha)*Y_l0(theta,phi) where x=cos theta is represented in the table.
void create_spherical_diffusion_kernel(const int tablesize, const int lmax, const double alpha, double values[], double * work1)
{
    double x,dx,coeff,pl,plm1,plm2,fact,oldfact,l2,xx;
    int l,i;
    for (i=0; i<tablesize; i++) values[i]=0.0;
    dx=2.0/tablesize;
    //work1=(double *) checkalloc((lmax+1)*sizeof(double));
    //work2=(double *) checkalloc(YSIZE[lmax]*sizeof(double));
    for (l=0; l<=lmax; l++) {
	x=l*(l+1)*alpha;
	xx=sqrt(((double)(2*l+1))/(4*M_PI));
	if (x<700.0) work1[l]=xx*exp(-x); else work1[l]=0.0;
    }
    for (i=0; i<tablesize; i++) {
        x=-1.0+dx*(i+0.5);
        //l=0
        plm2=1/sqrt(4*M_PI);
        values[i]+=plm2*work1[0];
        //l=1
        plm1=sqrt(3/(4*M_PI))*x;
        values[i]+=plm1*work1[1];
        oldfact=sqrt(3.0);
        for (l=2; l<=lmax; l++) {
            //eq. 6.7.9, p. 294, Numerical Recipes 3rd ed. with m-0 substituted
	    l2=(double) l;
	    l2=l2*l2;
            fact=sqrt(4.0-1.0/l2);
            pl=fact*(x*plm1-plm2/oldfact);
            values[i]+=pl*work1[l];
            plm2=plm1;
            plm1=pl;
            oldfact=fact;
        }
	if (values[i]<=-1e-6) {
	    printf("error: %d %.4f\n",i,values[i]);
	    die();
	}
	if (values[i]<0) values[i]=0; //for small values only
    }
    x=values[tablesize-1];
    for (i=0; i<tablesize; i++) values[i]/=x;
}


//The diffusion kernel for SO(3) is sum_j exp(-j(j+1)*alpha) (2j+1) sin((j+1/2)theta)/sin(theta/2)
//where thete is the overall angle for the rotation.
void create_so3_diffusion_kernel(const int tablesize, const int jmax, const double alpha, double values[], double * work1)
{
    double x,dx,sintheta2,costheta2,sintheta,costheta,sinj12theta,cosj12theta,sinj32theta,cosj32theta;
    int j,i,k;
    for (i=0; i<tablesize; i++) values[i]=0.0;
    dx=2.0/tablesize;
    //work1=(double *) checkalloc((lmax+1)*sizeof(double));
    //work2=(double *) checkalloc(YSIZE[lmax]*sizeof(double));
    for (j=0; j<=jmax; j++) {
        x=j*(j+1)*alpha;
        if (x<700.0) work1[j]=(2*j+1)*exp(-x); else work1[j]=0.0;
    }
    for (i=0; i<tablesize; i++) {
        costheta=-1.0+dx*(i+0.5);
        sintheta2=sqrt((1-costheta)/2); //sin(theta/2) = (1-cos(theta))/2
	costheta2=sqrt((1+costheta)/2); //cos(theta/2) = (1+cos(theta))/2
	sintheta=sqrt(1-costheta*costheta);
	sinj12theta=sintheta2;
	cosj12theta=costheta2;
        for (j=0; j<=jmax; j++) {
	    values[i]+=work1[j]*(sinj12theta/sintheta2);
	    cosj32theta=cosj12theta*costheta-sinj12theta*sintheta;
	    sinj32theta=cosj12theta*sintheta+sinj12theta*costheta;
	    cosj12theta=cosj32theta;
	    sinj12theta=sinj32theta;
        }
        if (values[i]<=-0.01) {
            printf("error: %d %.4f\n",i,values[i]);
            die();
	}
	if (values[i]<0) values[i]=0; //for small values only
    }
    x=values[tablesize-1];
    for (i=0; i<tablesize; i++) values[i]/=x;
}

//THESE SUBROUTINES FOR THE RMSD OVERLAY
//Perform a givens rotation to zero out a[n*i+j].  Multiply b by the same givens rotaton.
//Assume i<j.
void givens(int n, int i, int j, double * a, double * b)
{
    double aij,aii,ajj,aik,ajk,bki,bkj,tan2theta,cos2theta,sin2theta,costheta2,sintheta2,costheta,sintheta;
    int k;
    aij=a[n*i+j];
    aii=a[n*i+i];
    ajj=a[n*j+j];
    if (fabs(ajj-aii)<TINY) { //theta = pi/4
        cos2theta=0.0;
        sin2theta=1.0;
        costheta=ONESQRT2;
        sintheta=ONESQRT2;
	//it may not be numerically accurate, but sign should be correct
	//this seems not to help
	/*if ((aij*(ajj-aii))<0.0) { //theta = -pi/4!
		sin2theta=-1.0;
		sintheta=-ONESQRT2;
	}*/
        costheta2=0.5;
        sintheta2=0.5;
    } else {
        tan2theta=2.0*aij/(ajj-aii);
        cos2theta=1/sqrt(1+tan2theta*tan2theta);
        sin2theta=tan2theta*cos2theta;
        costheta2=(1.0+cos2theta)/2.0;
        costheta=sqrt(costheta2);
        sintheta2=(1.0-cos2theta)/2.0;
        sintheta=sqrt(sintheta2);
        if (tan2theta<0.0) sintheta=-sintheta;
    }
    a[n*i+i]=costheta2*aii-sin2theta*aij+sintheta2*ajj;
    a[n*j+j]=sintheta2*aii+sin2theta*aij+costheta2*ajj;
    a[n*i+j]=0.0;
    a[n*j+i]=0.0;
    for (k=0; k<n; k++) if ((k!=i) && (k!=j)) {
        aik=a[n*i+k];
        ajk=a[n*j+k];
        a[n*i+k]=costheta*aik-sintheta*ajk;
        a[n*j+k]=sintheta*aik+costheta*ajk;
        a[n*k+i]=a[n*i+k];
        a[n*k+j]=a[n*j+k];
    }
    for (k=0; k<n; k++) {
        bki=b[n*k+i];
        bkj=b[n*k+j];
        b[n*k+i]=costheta*bki-sintheta*bkj;
        b[n*k+j]=sintheta*bki+costheta*bkj;
    }
}

//Find maximum off-diagonal element in row i.
void maxind(int n, int i, double * a, int * jmax, double * elmax)
{
    int j;
    *elmax=0.0;
    for (j=i+1; j<n; j++) {
        if (fabs(a[n*i+j])>*elmax) {
            *jmax=j;
            *elmax=fabs(a[n*i+j]);
        }
    }
}
//Diagonalize a, provide eigenvectors in q
void jacobi(int n, double * a, double * q)
{
    int i,j,k,imax,jjmax;
    int jmax[MAX_MATRIX];
    double maxoffdiag[MAX_MATRIX],overallmaxoffdiag;
    //Initialize q to the identity matrix.
    for (i=0; i<n*n; i++) q[i]=0.0;
    for (i=0; i<n; i++) q[(n+1)*i]=1.0;
    //Find the maximum off diagonal element in each column.
    for (i=0; i<n; i++) maxind(n,i,a,&jmax[i],&maxoffdiag[i]);
    for (;;) {
        //Find the maximum off diagonal element.
        overallmaxoffdiag=0.0;
        for (i=0; i<n; i++) if (maxoffdiag[i]>overallmaxoffdiag) {
            imax=i;
            overallmaxoffdiag=maxoffdiag[i];
        }
        if (overallmaxoffdiag<1e-12) break; //converged!
        //it's imax, jmax[imax].  Do a givens rotation to make it zero.
        jjmax=jmax[imax];
        givens(n,imax,jjmax,a,q);
        //now update maximum indices for rows imax, jjmax
        maxind(n,imax,a,&jmax[imax],&maxoffdiag[imax]);
        maxind(n,jjmax,a,&jmax[jjmax],&maxoffdiag[jjmax]);
    }
}

//Best RMSD fit of coords2 relative to coords1.  r is the displacement vector, q the quaternion.
//Locates position of center of coords1 in the group of coordinates coords2.
void rmsd_fit(int natom, double * weight, double * coords1, double * coords2, double * center, double * orient, double * rmsd)
{
    int iatom,k,l,iev;
    double r[3][3],f[4][4],qq[4][4],ev,aux,totalweight,c1center[3],c2center[3],c1center2[3],q[4];
    totalweight=0.0;
    for (iatom=0; iatom<natom; iatom++) totalweight+=weight[iatom];
    //First find the center of each group of atoms, and the difference between them, which is net displacement.
    for (k=0; k<3; k++) {
        c1center[k]=0.0;
        c2center[k]=0.0;
    }
    for (iatom=0; iatom<natom; iatom++)
        for (k=0; k<3; k++) {
            c1center[k]+=weight[iatom]*coords1[3*iatom+k];
            c2center[k]+=weight[iatom]*coords2[3*iatom+k];
        }
    for (k=0; k<3; k++) {
        c1center[k]/=totalweight;
        c2center[k]/=totalweight;
    }
    //Now temporarily subtract off the center.
    for (iatom=0; iatom<natom; iatom++)
        for (k=0; k<3; k++) {
            coords1[3*iatom+k]-=c1center[k];
            coords2[3*iatom+k]-=c2center[k];
        }
    //Calculate covariance matrix, eq. 5 in Coutsias et al. JCC 25, 1849 (2004), with a weight weighting.
    //In the case of equal weights the matrix is divided by N compared to what's in the paper.
    //For some reason we need to transpose R.]
    for (k=0; k<3; k++) for (l=0; l<3; l++) r[k][l]=0.0;
    for (iatom=0; iatom<natom; iatom++) for (k=0; k<3; k++) for (l=0; l<3; l++)
        r[k][l]+=weight[iatom]*coords1[3*iatom+k]*coords2[3*iatom+l];
    for (k=0; k<3; k++) for (l=0; l<3; l++) r[k][l]/=totalweight;
    //The first term for best fit RMSD (equation after eq. 11)
    aux=0.0;
    for (iatom=0; iatom<natom; iatom++) for (k=0; k<3; k++) {
        aux+=weight[iatom]*coords1[3*iatom+k]*coords1[3*iatom+k];
        aux+=weight[iatom]*coords2[3*iatom+k]*coords2[3*iatom+k];
    }
    aux/=totalweight;
    matrix_to_quat(r,&orient[0],&ev);
    normalize_quat(orient); //Just to make sure.  Among other things ensures positive real part of quaternion.
    //normalize_quat(q);
    //conjugate_quat(q,orient);//For some reason we need this.  It probably compensates for having swapped something earlier.
    *rmsd=sqrt(fabs(aux-2*ev));
    //Restore original coordinates.
    for (iatom=0; iatom<natom; iatom++)
        for (k=0; k<3; k++) {
            coords1[3*iatom+k]+=c1center[k];
            coords2[3*iatom+k]+=c2center[k];
        }
    //The center of coords2 is given by (r2-r1) - R(orient) * (r1), where r1 and r2 are the centers of coords1 and coords2.
    rotate_vector_by_quat(orient,c1center,c1center2);
    for (k=0; k<3; k++) center[k]=c2center[k]-c1center2[k];
}

//Constructs coordinates for atom "l" using atoms i, j, k, and the bond length (k-l), angle (j-k-l), and dihedral (i-j-k-l)
void build_atom(double * coords, int i, int j, int k, int l, double bond, double angle, double dihedral)
{
    double rij[3],rjk[3],mat[3][3],rkl[3],rrjk,rrij,x[3],dot;
    int m;
    rrjk=0.0;
    rrij=0.0;
    for (m=0; m<3; m++) {
        rjk[m]=coords[3*k+m]-coords[3*j+m];
        rrjk+=rjk[m]*rjk[m];
        rij[m]=coords[3*j+m]-coords[3*i+m];
        rrij+=rij[m]*rij[m];
    }
    rrjk=sqrt(rrjk);
    rrij=sqrt(rrij);
    //Set up a coordinate system so that the x axis is a unit vector along j-k,
    //the x axis is perpendicular to j-k but in the i-j-k plane,
    //and the y axis is perpendicular to both.
    //Set up vectors as columns since matmul transposes the matrix.
    dot=0.0;
    for (m=0; m<3; m++) {
        rjk[m]=rjk[m]/rrjk;
        rij[m]=rij[m]/rrij;
        dot+=rij[m]*rjk[m];
    }
    for (m=0; m<3; m++) {
        mat[2][m]=rjk[m];
        mat[0][m]=(rij[m]-dot*rjk[m])/sqrt(1-dot*dot);
    }
    mat[1][0]=mat[2][2]*mat[0][1]-mat[2][1]*mat[0][2];
    mat[1][1]=mat[2][0]*mat[0][2]-mat[2][2]*mat[0][0];
    mat[1][2]=mat[2][1]*mat[0][0]-mat[2][0]*mat[0][1];
    /*mat[0][1]=mat[1][2]*mat[2][0]-mat[2][2]*mat[1][0];
    mat[1][1]=mat[2][2]*mat[0][0]-mat[0][2]*mat[2][0];
    mat[2][1]=mat[0][2]*mat[1][0]-mat[1][2]*mat[0][0];*/
    //Convert (bond,pi-angle,dihedral) from spherical to cartesian coordinates:
    //pi-angle b/c an angle of 180 degrees corresponds to the "north pole"
    //this is displacement vector in the new system.
    //Then multiply by mat^T (matmul does the transposition) to get to
    sph_to_cart(bond,M_PI-angle,M_PI-dihedral,&x[0]);
    matmul(&mat[0][0],x,rkl);
    for (m=0; m<3; m++) coords[3*l+m]=coords[3*k+m]+rkl[m];
}






