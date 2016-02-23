#ifndef ROTATIONS_H_INCLUDED
#define ROTATIONS_H_INCLUDED

#define RAD_TO_DEG                  57.2957795//converts radians to degrees
#define DEG_TO_RAD                  1.74532925e-2//converts degrees to radians

void fill_trig_tables(void);
//void prefetch_trig_tables(void);
void quat_to_euler(const double * q, double * phi, double * theta, double * psi);
void euler_to_quat(double phi, double theta, double psi, double * q);
void quat_to_matrix(double * q, double * r);
void matrix_to_quat(const double r[3][3],double * q);
void rotate_vector_by_quat(const double * q, const double * v1, double * v2);
void euler_to_matrix(double phi, double theta, double psi, double * r);
void axisangle_to_quat(const double alpha, const double * axis, double * q);
void quat_to_axisangle(double * q, double * alpha, double * axis);
void normalize_quat(double * q);
void multiply_quat(const double * qa, const double * qb, double * qc);
double dist(const double * qa, const double * qb);
void conjugate_quat(const double * qa, double * qb);
bool is_quat_ok(const double * q);
//void divide_quat(double * qa, double * qb, double * qc);
void multiply_conj_quat(double * qa, double * qb, double * qc);
void matmul(double * a, double * b, double * c);
void matmul2(double a[3][3], double b[3][3], double c[3][3]);
void rand_unif_quat(double *q);
void rand_small_quat(double delta, double *q);
void cart_to_sph(const double * x, const double r, double * sphtheta, double * sphphi);
void sph_to_cart(double r, double sphtheta, double sphphi, double * x);
void jacobi(int n, double * a, double * q);
void rmsd_fit(int natom, double * weight, double * coords1, double * coords2, double * center, double * orient, double * rmsd);
void create_spherical_diffusion_kernel(const int tablesize, const int lmax, const double alpha, double values[], double * work1);
void create_so3_diffusion_kernel(const int tablesize, const int jmax, const double alpha, double values[], double * work1);
void build_atom(double * coords, int i, int j, int k, int l, double bond, double angle, double dihedral);


/*inline double pbc_distance2(int pbc, double halfboxsize, double boxsize, double * icoords, double * jcoords, double * dx)
{
    double dx[3];
    int k;
            for (k=0; k<3; k++) {
                dx[k]=icoords[3*ifragatom+k]-jcoords[3*jfragatom+k];
                if (pbc) {
                    if (dx[k]<-halfboxsize) dx[k]+=boxsize;
                    if (dx[k]>halfboxsize) dx[k]-=boxsize;
                }
            }
            r2=dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2];
    return r2;
}*/

inline double pbc_distance2(int pbc, double halfboxsize, double boxsize, double * icoords, double * jcoords)
{
    double dx[3],r2;
    int k;
            for (k=0; k<3; k++) {
                dx[k]=icoords[k]-jcoords[k];
                if (pbc) {
                    if (dx[k]<-halfboxsize) dx[k]+=boxsize;
                    if (dx[k]>halfboxsize) dx[k]-=boxsize;
                }
            }
            r2=dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2];
    return r2;
}

#endif // ROTATIONS_H_INCLUDED
