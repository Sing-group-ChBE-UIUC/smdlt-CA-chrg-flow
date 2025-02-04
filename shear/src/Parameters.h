#ifndef Parameters_h
#define Parameters_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
// #include <complex.h>
// #include <stdbool.h>
// #include <stdint.h>
// #include <assert.h>
#include <string.h>
// #include <limits.h>
// #include <unistd.h>
// #include <ctype.h>
// #include <sys/time.h>
// #include <sys/types.h>
#include <omp.h>
#include <mkl.h>

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
#define NDIV (1+IMM1/NTAB)
#define max(a,b) ({ __typeof__ (a) _a = (a); __typeof__ (b) _b = (b); _a > _b ? _a : _b; })
#define min(a,b) ({ __typeof__ (a) _a = (a); __typeof__ (b) _b = (b); _a < _b ? _a : _b; })
#define EP 0.9624
#define THETA0 0

/////////////////////////////
// Global Variables
/////////////////////////////
///Charge Variables///
void coulombForce(); //applies coulomb interaction
void chargeHop() ; //charges hop with coulomb interactions
void chargeHop_nonint(); //charges hop without coloumb interactions
void calcMSD(); //calculates MSD
void calcDisFD(); //calculates displacement of FD charges
void calcDisHI(); //calculates displacement of charges with HI
void calcDis(); //calculates displacement of charges
void calcCH();
void initPos();//gets the initial position of charges
void calcCOM();//calculates center of mass position of each unwrapped chain
void calcCOM_r();//calculates center of mass position of each chain
void calcMSDcom();//calculates MSD of C.O.M
void calcRg();//calculates radius of gyration tensor
int *Charge;  //0 if uncharged //1 if charged
int *Charge_track ;  //tracks charge
int *Charge_old ; //tracks charge
int *Charge_indices ; //tracks charge
int *check_ch ; //checks if charge hopped thru PBC in dx,dy, and dz | 0 - regular displacement 1 - PBC displacement
int Ncharges ; //Number of charged beads
int i,cnumber,counter1,counter2,counter3, counter4, counter5, counter6;
double lambda_d ; //debye screening length
double lambda_b ; //bjerrum length
double kappa_d ; //inverse debye screening length(1/lambda_d)
int chargehop_type ; //0 - no charge hop, 1 - interacting charge hop, 2 - non-interacting charge hop
int MStep ; // do MC every MSteps
int *i_values ;
double barrier ; //barrier for adjacent hopping
double barrier2 ; //barrier for non-adjacent hopping
char* str ;
char* str2 ;
char* map;
FILE *outputfile1, *chfile;
double Lbar[9]; // cell vectors {L1,L2,L3}
char *outp1;
FILE *datafile,*outputfile2;
int printprops;
double MSD_ch;
double MSD_com;
double conc; // concentration normalized by c_star
long double *dxt, *dyt, *dzt;
double *rbi; //initial position of charge at timestep
int printprops;
int tcount;
unsigned long tmax_temp;
int chargehop_type_temp;
double E_initial; //initial coulomb energy before charge hop
FILE *mapfile; //file that stores mapped positions
double *sx,*sy,*sz; //KRBC positions mapped to a unit cube
double matrix[3][3], inv[3][3] ; //hmap and hmap^-1
double *rbic; //initial position of charge
double *comx_r,*comy_r,*comz_r ;  //center of mass of real bead positions
int flowtype ; 
int MSDstart ; 
char spring[100];
double tauxy;
char *rgt;
FILE *rgtfile;
//////////////////////
double c_norm; // concentration normalized by c_star
int Nb; // number of beads per chain
int Nc; // number of chains
int N; // total number of beads
int Ncm; // Number of centers of mass per chain
unsigned long tmax,t,ttemp,tstart,equil_start,t_relax,t_cess,t_m_tcess,tau_relax; // Max number of time steps
double tau_FD;
double epsilon; // LJ interaction strength
double kappas; // Spring constant
double kappab; // Bending potential constant
double Rg; // Root mean squared radius of gyration
double dt,dt_eq; // Time step
int restart,fixSeed; // Restart condition. 0 -> new run
double c_star; // Overlap concentration
double box_volume; // Volume of the box given N and c_norm
double L; // Length of the box from the volume
double box_side; // Location of the side of the box
double sigma; // LJ hard-core overlap distance = 2a
double rc; // LJ cutoff radius
double rv,rv2; // Verlet list cutoff radius
double r2cut;
double rnew; // Distance which determines when to update neighbor list
int *nlist; // Number of neighbors in neighbor list for each bead
int *list; // Neighbor list in format (bead,neighbor indices)
double *xo,*yo,*zo;
double *dxv,*dyv,*dzv;
double *dxu,*dyu,*dzu;
double *drb;
// double *dr; // Displace from bead position the last time neighbor tables were updated
double *rx,*ry,*rz; // Position arrays, eg rx[0] = x-pos of bead 0, ry[4] = y-pos of bead 4
double *rb;
double *px,*py,*pz;
double *rgx,*rgy,*rgz,*rg;
double *comx,*comy,*comz;
double *reex,*reey,*reez;
double *rext;
double tauxx,tauyy;
double *fx,*fy,*fz; // Force arrays
double *f;
double *R,*Z,*Y; // Gaussian dist random velocities, R[0:2] = x,y,z random velocities of bead 0
double p,pc,pc_eq; // Width of the Gaussian dist random velocity to satisfy fluctuation-dissipation
int srpy;
long *idum; // Random number generator seed
char *xyz; // xyz filename
char *outp; // output filename
char *clustero;
char *mmatrix;
char *decomp;
char *ext;
char *grid;
char *visc;
char *tim;
char *beta;
FILE *outputfile, *xyzfile,*clusterout,*MMfile,*DCfile,*extfile,*gridfile,*viscfile,*timfile,*betafile;
int num_threads;
int trace;
int printperiod;
int spp;
unsigned long *tprint;
double ljtimetest,ljtime;
double flowrate;
double point[3][2], point0[3][2];
double L1[2], L2[2]; //vector of the box sides
int tp;
typedef struct {
    double x,y,z;
} Vector3D_t;
double L1xp,L1yp,L2xp,L2yp;
double detL; // Determinant of the lattice basis matrix
double thetabox,L2pmag;
int ebins; // Number of extension bins for GW average, interval = 1.0/ebin
double ext_avg,ext_interval; // Average extension for a given time step - determines which GW parameters to use
int sbins,samp_per_inter,strain_sample,rbins,r_interval,r_sample; // stain bins, samples per interval
double strain_ss,strain_interval; // steady state strain and interval
double strain,total_strain;
double **C,*Ct; // Average TEA decomposition C_i's from previous iteration
double **D; // Grid-space average HI from pre-calculation gridfile
double ***An;
double *Ds; // Intrachain HI average from the previous iteration
double *Dsrun; // Running intrachain HI average for the current iteration
double *Drun; // Running grid-space HI average for the current iteration
double **Crun; // Running average TEA decomposition C_i's for the current iteration
double *rowsum;
double *M; // Full diffusion tensor for use in Langevin update
int *count;
int selfcount;
int *gwcount;
double eps,eps2,*bp,*brun,bt;
double bin_size_x,bin_size_y,bin_size_z,bin_size;
int num_bins_x,num_bins_y,num_bins_z,num_bins,nhx,nhy,nhz;
int itermax,iterstart,itercount;
int xyzstart;
int *start_ij;
double mult_time,bin_time;
int nt_avg;
double total_time;
double H,qmax,qmax2;
int red_ts;

int overstretch; // 0 -> time step accepted, 1 -> reject time step, bond overstretched

#endif // Parameters.h
