#include <stdlib.h>
#ifndef Constants_h
#define Constants_h

static const int SIZEPOP = 480;
static const int GENERATIONS = 4000;
static const int AVERAGE_OVER = 10;

static const int nTarget = 3; //number of targets to be optimized (1. Computational Error, 2. Number of connections, 3. Regeneration rate)

//Networks' topology
static const int MAXSIZE = 15; //maximum number of units in each hidden layer
static const int INPUT = 3; //nº input units
static const int OUTPUT = 1; //nº output units
static const int LAYERS = 2; //Hidden neurons layers

#endif
