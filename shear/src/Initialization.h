#ifndef Initialization_h
#define Initialization_h

#include "Parameters.h"
#include <getopt.h>

void Initialization(int argc, char * argv[]); // Perform initialization operations

void ParseInput(int argc, char * argv[]); // Parse command line inputs for variables

void readInput(); // Read constant inputs from file, eg Input.txt

void initBox(); // Define the box dimensions given number of beads and concentration

void initVerlet(); // Define the verlet neighbor list cutoffs

void initCell(); // Initial linked list cell variables

void initEwald();

void allocate(); // Allocate memory for the position, force, mobility matrix etc arrays.
// Also other miscellaneous constant definitions here.

void printOutput(); // Prints simulation parameters to the output file

void initLattice();

pid_t getpid(void); // Just here for the initRan function

long initRan(); // Initialize the RNG.

float ran1(long *idum); // Uniform RNG [0,1) from numerical recipes textbook

#endif // Initialization_h
