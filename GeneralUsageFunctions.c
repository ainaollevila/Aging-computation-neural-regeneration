#include "GeneralUsageFunctions.h"
#include "Constants.h"

void SortSmalltoBig(int *vec_position, int *vec_to_order, int *vec_ordered, int size){
    int i,j;
    int aux1,aux2;
    
    for (i=0;i<SIZEPOP;i++)
        vec_ordered[i] = vec_to_order[i];
    
    for (i=0;i<SIZEPOP;i++)
        vec_position[i] = i;
    
    for (i=0;i<SIZEPOP;i++)
        for(j=i+1;j<SIZEPOP;j++)
        {
            if (vec_ordered[i]>vec_ordered[j])
            {
                aux1 = vec_ordered[j];
                vec_ordered[j] = vec_ordered[i];
                vec_ordered[i] = aux1;
                
                aux2 = vec_position[j];
                vec_position[j] = vec_position[i];
                vec_position[i] = aux2;
            }
        }
};

void StepFunction(int rows, int cols, double *u, int *v, int *v_inputs, int *omega, double *theta){
    int i,j;
    for (j=0;j<cols;j++)
    {
        for (i=0;i<rows;i++)
            u[j] += v_inputs[i]*omega[i*cols+j];
        u[j]-= theta[j];
        if (u[j]>=0) v[j] = 1; else v[j]=0;
    }
};
