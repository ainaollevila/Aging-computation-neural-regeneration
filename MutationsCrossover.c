#include "MutationsCrossover.h"
#include <stdlib.h>
#include "Constants.h"
#include "GlobalVariables.h"
#include <stdbool.h>

void Crossover (int *matrix1, int i1, int j1, int *matrix2, int i2, int j2){
    int i,j,k,l;
    int *auxmatrix;
    auxmatrix = (int*)malloc(sizeof(int)*MAXSIZE*MAXSIZE);
    for (i=0;i<MAXSIZE*MAXSIZE;i++)
        auxmatrix[i] = 0;
    int common_hid1, common_hid2;
    int dif_row, dif_col;
    int rand_row, rand_col;
    
    
    if (i1==i2 && j1 == j2)
        for (i=0;i<i1*j1;i++)
        {
            auxmatrix[i]= matrix1[i];
            matrix1[i] = matrix2[i];
            matrix2[i] = auxmatrix[i];
        }
    
    if (i1==i2 && j1 > j2)
    {
        dif_col = j1-j2;
        rand_col = rand()%(dif_col+1);
        common_hid1 = i1;
        
        for (i=0;i<common_hid1;i++)
        {
            k=0;
            for(j=rand_col;j<(rand_col+j2);j++)
            {
                auxmatrix[i*j2+k] = matrix1[i*j1+j];
                matrix1[i*j1+j] = matrix2[i*j2+k];
                matrix2[i*j2+k] = auxmatrix[i*j2+k];
                k++;
            }
        }
    }
    
    if (i1==i2 && j1 < j2)
    {
        dif_col = j2-j1;
        rand_col = rand()%(dif_col+1);
        common_hid1 = i1;
        
        for (i=0;i<common_hid1;i++)
        {
            k=0;
            for(j=rand_col;j<(rand_col+j1);j++)
            {
                auxmatrix[i*j1+k] = matrix2[i*j2+j];
                matrix2[i*j2+j] = matrix1[i*j1+k];
                matrix1[i*j1+k] = auxmatrix[i*j1+k];
                k++;
            }
        }
    }
    
    if (i1>i2 && j1 == j2)
    {
        dif_row = i1-i2;
        rand_row = rand()%(dif_row+1);
        common_hid2 = j1;
        l=0;
        for (i=rand_row;i<(rand_row+i2);i++)
        {
            for(j=0;j<common_hid2;j++)
            {
                auxmatrix[l*common_hid2+j] = matrix1[i*common_hid2+j];
                matrix1[i*common_hid2+j] = matrix2[l*common_hid2+j];
                matrix2[l*common_hid2+j] = auxmatrix[l*common_hid2+j];
            }
            l++;
        }
    }
    
    if (i1<i2 && j1 ==j2)
    {
        dif_row = i2-i1;
        rand_row = rand()%(dif_row+1);
        common_hid2 = j1;
        l=0;
        for (i=rand_row;i<(rand_row+i1);i++)
        {
            for(j=0;j<common_hid2;j++)
            {
                auxmatrix[l*common_hid2+j] = matrix2[i*common_hid2+j];
                matrix2[i*common_hid2+j] = matrix1[l*common_hid2+j];
                matrix1[l*common_hid2+j] = auxmatrix[l*common_hid2+j];
                
            }
            
            l++;
        }
    }
    
    if (i1 < i2 && j1 < j2)
    {
        dif_col = j2-j1;
        dif_row = i2-i1;
        
        rand_col = rand()%(dif_col+1);
        rand_row = rand()%(dif_row+1);
        
        l=0;
        for(i=rand_row;i<(rand_row+i1);i++)
        {
            k=0;
            for(j=rand_col;j<(rand_col+j1);j++)
            {
                auxmatrix[l*j1+k]=matrix2[i*j2+j];
                matrix2[i*j2+j]=matrix1[l*j1+k];
                matrix1[l*j1+k] = auxmatrix[l*j1+k];
                k++;
            }
            l++;
        }
    }
    
    if (i1 > i2 && j1 > j2)
    {
        dif_col = j1-j2;
        dif_row = i1-i2;
        
        rand_col = rand()%(dif_col+1);
        rand_row = rand()%(dif_row+1);
        
        l=0;
        for(i=rand_row;i<(rand_row+i2);i++)
        {
            k=0;
            for(j=rand_col;j<(rand_col+j2);j++)
            {
                auxmatrix[l*j2+k]=matrix1[i*j1+j];
                matrix1[i*j1+j]=matrix2[l*j2+k];
                matrix2[l*j2+k] = auxmatrix[l*j2+k];
                k++;
            }
            l++;
        }
    }
    
    if (i1 < i2 && j1 > j2)
    {
        dif_col = j1 - j2;
        dif_row = i2 - i1;
        
        rand_col = rand() % (dif_col+1);
        rand_row = rand() % (dif_row+1);
        
        l=0;
        for(i=rand_row;i<(rand_row+i1);i++)
        {
            k=0;
            for (j=rand_col;j<(rand_col+j2);j++)
            {
                auxmatrix[l*j2+k] = matrix1[l*j1+j];
                matrix1[l*j1+j] = matrix2[i*j2+k];
                matrix2[i*j2+k]= auxmatrix[l*j2+k];
                k++;
            }
            l++;
        }
    }
    
    if (i1 > i2 && j1 < j2)
    {
        dif_col = j2 - j1;
        dif_row = i1 - i2;
        
        rand_col = rand() % (dif_col+1);
        rand_row = rand() % (dif_row+1);
        
        l=0;
        for(i=rand_row;i<(rand_row+i2);i++)
        {
            k=0;
            for (j=rand_col;j<(rand_col+j1);j++)
            {
                auxmatrix[l*j1+k] = matrix2[l*j2+j];
                matrix2[l*j2+j] = matrix1[i*j1+k];
                matrix1[i*j1+k]= auxmatrix[l*j1+k];
                k++;
            }
            l++;
        }
    }
};

void NeuronsMutation(network **network1, double In, double Dn){
    int i;
    bool validmutations;
    int *neurons_layer, *new_neurons_layer;
    neurons_layer = (int*)malloc(sizeof(int)*(LAYERS+2));
    new_neurons_layer = (int*)malloc(sizeof(int)*(LAYERS+2));
    for (i=0;i<LAYERS+2;i++)
    {
        neurons_layer[i]=0;
        new_neurons_layer[i]=0;
    }
    int *duplicated_neurons, *eliminated_neurons;
    duplicated_neurons = (int*)malloc(sizeof(int)*MAXSIZE*(LAYERS+2));
    eliminated_neurons = (int*)malloc(sizeof(int)*MAXSIZE*(LAYERS+2));
    for(i=0;i<MAXSIZE*(LAYERS+2);i++)
    {
        duplicated_neurons[i] = 0;
        eliminated_neurons[i] = 0;
    }
    int *protected_neurons;
    protected_neurons = (int*)malloc(sizeof(int)*(LAYERS+2));
    for(i=0;i<LAYERS+2;i++) protected_neurons[i] = 0;
    
    protected_neurons[0]=(*network1)->vec_feas_paths[(*network1)->protected_path].input;
    protected_neurons[1]=(*network1)->vec_feas_paths[(*network1)->protected_path].hid1;
    protected_neurons[2]=(*network1)->vec_feas_paths[(*network1)->protected_path].hid2;
    protected_neurons[3]=(*network1)->vec_feas_paths[(*network1)->protected_path].output;
    
    neurons_layer[0] = INPUT;
    neurons_layer[1] = (*network1)->hid1;
    neurons_layer[2] = (*network1)->hid2;
    neurons_layer[3] = OUTPUT;
    
    NeuronsMut(duplicated_neurons, eliminated_neurons, neurons_layer, new_neurons_layer, protected_neurons, Dn, In);
    
    validmutations = 1;
    
    for (i=1;i<=LAYERS;i++)
        if (new_neurons_layer[i]>MAXSIZE)
            validmutations = 0;
    
    if (validmutations)
    {
        NeuronsMutMatricesRestablish((*network1)->inhid1, 0, 1, duplicated_neurons, eliminated_neurons, neurons_layer, new_neurons_layer);
        NeuronsMutMatricesRestablish((*network1)->hid1hid2, 1, 2, duplicated_neurons, eliminated_neurons, neurons_layer, new_neurons_layer);
        NeuronsMutMatricesRestablish((*network1)->hid2out, 2, 3, duplicated_neurons, eliminated_neurons, neurons_layer, new_neurons_layer);
        
        NeuronsMutActThreshRestablish(1, (*network1)->theta_hid1, neurons_layer, new_neurons_layer, duplicated_neurons, eliminated_neurons);
        NeuronsMutActThreshRestablish(2, (*network1)->theta_hid2, neurons_layer, new_neurons_layer, duplicated_neurons, eliminated_neurons);
        
        (*network1)->hid1 = new_neurons_layer[1];
        (*network1)->hid2 = new_neurons_layer[2];
    }
    
    
    free(neurons_layer);
    free(new_neurons_layer);
    free(duplicated_neurons);
    free(eliminated_neurons);
    free(protected_neurons);
};

void NeuronsMut(int *d_neurons, int *e_neurons, int *neurons_layer, int *new_neurons_layer, int *protected_neurons, double Dn, double In){
    int i,j;
    double rand1;
    
    for(i=1;i<=LAYERS;i++)
        for (j=0;j<neurons_layer[i];j++)
        {
            
            rand1 = (double)rand()/RAND_MAX;
            if (j==protected_neurons[i]&&rand1<Dn)
                d_neurons[i*MAXSIZE+j] = 1;
            if (j!=protected_neurons[i])
            {
                if (rand1<In) //Eliminate node
                    e_neurons[i*MAXSIZE+j] = 1;
                else if (rand1 < (In+Dn)) //Duplicate node
                    d_neurons[i*MAXSIZE+j] = 1;
            }
        }
    
    for (i=0;i<LAYERS+2;i++)
        new_neurons_layer[i] = neurons_layer[i];
    
    for (i=1;i<=LAYERS;i++)
        for (j=0;j<neurons_layer[i];j++)
        {
            if (e_neurons[i*MAXSIZE+j] == 1) new_neurons_layer[i]--;
            if (d_neurons[i*MAXSIZE+j] == 1) new_neurons_layer[i]++;
        }
    
};

void NeuronsMutMatricesRestablish(int *weights_matrix, int layer_from, int layer_to, int *d_neurons, int *e_neurons, int *neurons_layer, int *new_neurons_layer){
    int i,j,k,l;
    int nn, mm;
    nn=0;
    mm=0;
 
    int *aux_weights;
    aux_weights = (int*)malloc(sizeof(int)*MAXSIZE*MAXSIZE);
    for (i=0;i<MAXSIZE*MAXSIZE;i++) aux_weights[i]=0;
    
    for (i=0;i<neurons_layer[layer_from];i++)
    {
        if(d_neurons[layer_from*MAXSIZE+i] == 1)
            for (k=0;k<2;k++)
            {
                mm=0;
                for (j=0;j<neurons_layer[layer_to];j++)
                {
                    if (d_neurons[layer_to*MAXSIZE+j] == 1)
                    {
                        for (l=0;l<2;l++)
                        {
                            aux_weights[nn*new_neurons_layer[layer_to]+mm]=weights_matrix[i*neurons_layer[layer_to]+j];
                            mm++;
                        }
                    }
                    else if (e_neurons[layer_to*MAXSIZE+j] == 0)
                    {
                        aux_weights[nn*new_neurons_layer[layer_to]+mm] = weights_matrix[i*neurons_layer[layer_to]+j];
                        mm++;
                    }
                }
                nn++;
            }
        
        else if (e_neurons[layer_from*MAXSIZE+i] == 0)
        {
            mm=0;
            for (j=0;j<neurons_layer[layer_to];j++)
            {
                if (d_neurons[layer_to*MAXSIZE+j] ==1)
                {
                    for (l=0;l<2;l++)
                    {
                        aux_weights[nn*new_neurons_layer[layer_to]+mm]=weights_matrix[i*neurons_layer[layer_to]+j];
                        mm++;
                    }
                }
                else if (e_neurons[layer_to*MAXSIZE+j] == 0)
                {
                    aux_weights[nn*new_neurons_layer[layer_to]+mm] = weights_matrix[i*neurons_layer[layer_to]+j];
                    mm++;
                }
            }
            nn++;
        }
    }
    
    for (i=0;i<new_neurons_layer[layer_from]*new_neurons_layer[layer_to];i++)
        weights_matrix[i]=aux_weights[i];
    
    free(aux_weights);
};

void NeuronsMutActThreshRestablish(int layer, double *act_thresh_vec, int *neurons_layer, int *new_neurons_layer, int *d_neurons, int *e_neurons){
    int i,k,l;
    double *new_aux;
    new_aux = (double*)malloc(sizeof(double)*MAXSIZE);
    for (i=0;i<MAXSIZE;i++)
        new_aux[i] = 0;
    l=0;
    for (i=0;i<neurons_layer[layer];i++)
    {
        if (d_neurons[layer*MAXSIZE+i] == 1)
            for (k=0;k<2;k++)
            {
                new_aux[l] = act_thresh_vec[i]; l++;
            }
        else if (e_neurons[layer*MAXSIZE+i] == 0)
        {
            new_aux[l] = act_thresh_vec[i]; l++;
        }
    }
    for (i=0;i<new_neurons_layer[layer];i++) act_thresh_vec[i] = new_aux[i];
    
    free(new_aux);
};

void ConnectionsMutation(int *weights_matrix, int n_from, int n_to, int nprotected_from, int nprotected_to ,double Cc, double Ec, double Ca){
    
    double rand1, rand2;
    int i,j;
    for (i=0;i<n_from;i++)
        for (j=0;j<n_to;j++)
        {
            if (weights_matrix[i*n_to+j] == 0 )
            {
                rand1 = (double)rand()/RAND_MAX;
                rand2 = (double)rand()/RAND_MAX;
                if (rand1 < Cc)
                {
                    if (rand2<0.5)
                        weights_matrix[i*n_to+j]=1;
                    else
                        weights_matrix[i*n_to+j]=-1;
                }
            }
            else
            {
                rand1 = (double)rand()/RAND_MAX;
                if (i==nprotected_from && j==nprotected_to)
                {
                    if (rand1 <Ca)
                        weights_matrix[i*n_to+j] *= -1;
                }
                else
                {
                    if (rand1< Ec) weights_matrix[i*n_to+j]=0; else if (rand1<(Ec+Ca)) weights_matrix[i*n_to+j] *= -1;
                }
            }
        }
};

void ActivationThresholdMutation(double *thresh_vec, int size_vec, double Ctheta, double tau){
    double rand1, rand3, rand2;
    int i;
    double eta;
    double low_boundary = 0.01; 
    
    for (i=0;i<size_vec;i++)
    {
        rand2 = (double) rand()/RAND_MAX;
        if (rand2 < Ctheta)
        {
            rand1 = (double) rand()/RAND_MAX;
            eta = rand1*tau;
            rand3 = (double) rand()/RAND_MAX;
            if (rand3 < 0.5) eta *= -1;
            thresh_vec[i] *= (1+eta);
            if (thresh_vec[i] > 1) thresh_vec[i] = 1;
            if (thresh_vec[i] < low_boundary) thresh_vec[i] = low_boundary;
        }
    }
};

void RegenerationRateMutation(double *reg_rate, double tau2, double Creg){
    double rand1, rand2, rand3;
    double eta;
    
    rand2 = (double) rand()/RAND_MAX;
    if (rand2 < Creg)
    {
        rand1 = (double) rand()/RAND_MAX;
        eta = rand1*tau2;
        rand3 = (double) rand()/RAND_MAX;
        if (rand3 < 0.5) eta *= -1;
        *reg_rate *= (1+eta);
    }
    if (*reg_rate > 1) *reg_rate = 1;
    if (*reg_rate < minimum_regeneration_rate) *reg_rate = minimum_regeneration_rate;
};
