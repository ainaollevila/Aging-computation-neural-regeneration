#include "TargetsComputation.h"
#include <stdlib.h>

double ExpectedObservedOutputAccuracy(int twotoinput, int time, int *obs_out, int *exp_out){
    double target=0,fit_aux=0,factor1;
    int i,j;
    for (i=0;i<time;i++)
        for (j=0;j<twotoinput;j++)
            fit_aux += abs(exp_out[j] - obs_out[i*twotoinput+j]);
    factor1 = (double)1/(twotoinput*time);
    target = fit_aux*factor1;
    return target;
};

double CountConnections (int rows, int cols, int *feas_neurons_rows, int *feas_neurons_cols, int *weights_matrix){
    int i,j;
    double num_connections = 0;
    for(i=0;i<rows;i++)
        for(j=0;j<cols;j++)
            if (feas_neurons_rows[i] == 1 && feas_neurons_cols[j]==1 && weights_matrix[i*cols+j] !=0)
                num_connections++;
    return num_connections;
};

