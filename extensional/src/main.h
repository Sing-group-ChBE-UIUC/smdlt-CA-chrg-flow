#ifndef main_h
#define main_h

#include "Parameters.h"
#include "Initialization.c"

void initChains();

struct bin {int x; int y; int z;};

void resetForce(); // Resets forces on beads after each time step

float gasdev(long *idum); // Box-Mueller Transform output from ran1 (uniform) to a Gaussian Dist

void bondforce(); // Bonded forces. Normally stretching and bending given by kappas, kappab.

void LJforce(); // LJ EV forces given by epsilon

void verletlist(); // Function for creating neighbor lists. See Frenkel & Smit p. 545 for details

void ewaldBin();

void CATEA();

void printTEA();

void getNoise(); // Get Brownian random velocity from the Gaussian dist function

void updateEquil();

void updateFD();

void updateHI();

void checkVerlet(); // Check if neighbor list needs to be updated. Again see See Frenkel & Smit p. 545 for details

void printTrajectory(); // Save trajectory coordinates to the xyz file.

void printRestart();

void calcExt();

void printExt();

void updateLattice(unsigned long t);

void applyPBC();

Vector3D_t getNID(int i,int j);

Vector3D_t getNIDbin(int i,int j);

void transformCoords();

void inverseCoords();

int binij(int i,int j);

void calcHI();

void resetAverage();

void printVisc();

void printTiming();

#endif // main_h
