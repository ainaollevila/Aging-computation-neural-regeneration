#ifndef Aging_h
#define Aging_h

#include <stdio.h>

void Damage(int rows, int cols, int *conn_mat, int *conn_mat_damaged, int *feas_neurons_origin, int *feas_neurons_dest);
void Regeneration(int rows, int cols, int *conn_mat, int *conn_mat_damaged, double reg_rate);
#endif 
