#ifndef grdHI_h
#define grdHI_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>

#define EP 0.9624
#define THETA0 31.7

typedef struct {
    double x,y,z;
} Vector3D_t;

void ParseInput(int argc, char * argv[]); // Parse command line inputs for variables
void readInput(); // Read constant inputs from file, eg Input.txt
void initBox(); // Define the box dimensions given number of beads and concentration
void initLattice();
void allocate();
void updateLattice(unsigned long t);
Vector3D_t invbin();

unsigned long t,t1,t2;
double rce;
int nt_avg;
double M_HI;
int nmax;
double alpha,alpha2,alpha3,alpha4,alpha5,alpha7;
double selfconst,rtpi;
int kxmax,kymax,kzmax;
double kcoeff;
double rkcut,rkkcut;
double L,L1xp,L1yp,L2xp,L2yp,detL,box_volume,box_side;
double point[3][2], point0[3][2];
double L1[2], L2[2]; //vector of the box sides
double c_norm,flowrate,Rg,c_star,dt;
int N,Nb,Nc,num_threads,trace,restart,tp;
unsigned long tmax;
FILE *inputfile,*mmfile;
char *mm;
double **D;
double xmax,ymax,zmax,theta,bin_size,bin_size_x,bin_size_y,bin_size_z;
int nhx,nhy,nhz,num_bins,origin,num_bins_x,num_bins_y,num_bins_z;
int treal;

#endif