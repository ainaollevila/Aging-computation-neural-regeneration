#include "Aging.h"
#include "GlobalVariables.h"
#include <stdlib.h>
void Damage(int rows, int cols, int *conn_mat, int *conn_mat_damaged, int *feas_neurons_origin, int *feas_neurons_dest){
    int i,j;
    double rand1;
    for (i=0;i<rows;i++)
        for (j=0;j<cols;j++)
        {
            if (conn_mat[i*cols+j] != 0 && feas_neurons_origin[i]==1 && feas_neurons_dest[j]==1)
            {
                rand1 = (double)rand()/RAND_MAX;
                if (rand1 < damage_rate)
                {
                    conn_mat_damaged[i*cols+j] = conn_mat[i*cols+j];
                    conn_mat[i*cols+j] = 0;
                }
            }
        }
};

void Regeneration(int rows, int cols, int *conn_mat, int *conn_mat_damaged, double reg_rate){
    int i,j;
    double rand1;
    for (i=0;i<rows;i++)
        for (j=0;j<cols;j++)
        {
            if (conn_mat_damaged[i*cols+j]!=0)
            {
                rand1 = (double)rand()/RAND_MAX;
                if (rand1 < reg_rate)
                {
                    conn_mat[i*cols+j] = conn_mat_damaged[i*cols+j];
                    conn_mat_damaged[i*cols+j] = 0;
                }
                
            }
        }
};
