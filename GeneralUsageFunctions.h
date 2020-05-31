#ifndef __RegenerationNN__GeneralUsageFunctions__
#define __RegenerationNN__GeneralUsageFunctions__

#include <stdio.h>

void SortSmalltoBig(int *vec_position, int *vec_to_order, int *vec_ordered, int size);
void StepFunction(int isize, int jsize, double *u, int *v, int *v_inputs, int *omega, double *theta);

#endif 
