#ifndef Feasibility_h
#define Feasibility_h

#include <stdio.h>
#include "Structs.h"

void InputsFeasibility(int hid1, int *weights_inhid1, int *feas_hid1, int totalfeas1);
void ValidNeurons(int layer, int n_pre, int n_actual, int n_post, int *weights_preactual, int *weights_actualpost, int *non_valid_neurons);
void CheckFeasibility(int twotoinput, network *network1);
void CheckFeasibilityPTP(int twotoinput, network **network1);
#endif 
