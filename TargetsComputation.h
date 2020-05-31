#ifndef TargetsComputation_h
#define TargetsComputation_h

#include <stdio.h>

double ExpectedObservedOutputAccuracy(int twotoinput, int time, int *obs_out, int *exp_out);
double CountConnections (int rows, int cols, int *feas_neurons_rows, int *feas_neurons_cols, int *weights_matrix);

#endif 
