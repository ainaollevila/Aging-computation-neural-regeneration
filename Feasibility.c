#include "Feasibility.h"
#include <stdlib.h>
#include "Constants.h"
#include <stdbool.h>

void InputsFeasibility(int hid1, int *weights_inhid1, int *feas_hid1, int totalfeas1)
{
    int i,j;
    int randhid1,count;
    double rand1;
    bool feas_input = 0;
    
    for (i=0;i<INPUT;i++)
    {
        feas_input = 0;
        for (j=0;j<hid1;j++)
        {
            if(weights_inhid1[i*hid1+j]!=0 && feas_hid1[j]==1)
                feas_input = 1;
        }
        count = 0;
        if (!feas_input)
        {
            randhid1 = (rand() % totalfeas1) + 1;
            for (j=0;j<hid1;j++)
            {
                if (feas_hid1[j] == 1)
                    count++;
                if (count == randhid1)
                {
                    rand1 = (double) rand()/RAND_MAX;
                    if (rand1<0.5)
                        weights_inhid1[i*hid1+j] = 1;
                    else
                        weights_inhid1[i*hid1+j] = -1;
                }
                
            }
        }
    }
};

void ValidNeurons(int layer, int n_pre, int n_actual, int n_post, int *weights_preactual, int *weights_actualpost, int *non_valid_neurons){
    int i,j;
    int aux_nv_in, aux_nv_out; //nv= non valid
    for (j=0;j<n_actual;j++)
    {
        aux_nv_in = 0;
        aux_nv_out = 0;
        for (i=0;i<n_pre;i++)
            if (weights_preactual[i*n_actual+j] != 0)
                aux_nv_in = 1;
        for (i=0;i<n_post;i++)
            if (weights_actualpost[j*n_post+i] != 0)
                aux_nv_out = 1;
        if (aux_nv_in == 0 || aux_nv_out==0)
            non_valid_neurons[layer*MAXSIZE+j] = 1;
    }
};

void CheckFeasibility(int twotoinput, network *network1){
    
    int i,j,k;
    int *non_valid_neurons; //1=non-feasible. 0 = feasible.
    non_valid_neurons = (int*)malloc(sizeof(int)*LAYERS*MAXSIZE);
    for (i=0;i<LAYERS*MAXSIZE;i++)
        non_valid_neurons[i]=0;
    
    int input1_nc,input2_nc,input3_nc;
    
    ValidNeurons(0, INPUT, network1->hid1, network1->hid2, network1->inhid1, network1->hid1hid2, non_valid_neurons);
    ValidNeurons(1, network1->hid1, network1->hid2, OUTPUT, network1->hid1hid2, network1->hid2out, non_valid_neurons);
    
    network1->truly_feas_hid1 = 0;
    network1->truly_feas_hid2 = 0;
    for (i=0; i<MAXSIZE; i++)
    {
        network1->vec_aux_truly_feas_hid1[i] = 0;
        network1->vec_aux_truly_feas_hid2[i] = 0;
    }
    
    network1->feas_paths = 0;
    for (k=0;k<INPUT;k++)
        for (i=0;i<network1->hid1;i++)
            if (non_valid_neurons[0*MAXSIZE+i] == 0 && network1->inhid1[k*network1->hid1+i]!=0)
                for (j=0;j<network1->hid2;j++)
                    if (non_valid_neurons[1*MAXSIZE+j] == 0 && network1->hid1hid2[i*network1->hid2+j]!=0)
                    {
                        network1->vec_aux_truly_feas_hid1[i] = 1;
                        network1->vec_aux_truly_feas_hid2[j] = 1;
                        network1->vec_feas_paths[network1->feas_paths].input = k;
                        network1->vec_feas_paths[network1->feas_paths].hid1 = i;
                        network1->vec_feas_paths[network1->feas_paths].hid2 = j;
                        network1->vec_feas_paths[network1->feas_paths].output = 0; //valid only if one Output unit
                        network1->feas_paths++;
                    }

    for (i=0; i<network1->hid1; i++) if (network1->vec_aux_truly_feas_hid1[i] == 1) network1->truly_feas_hid1++;
    for (i=0; i<network1->hid2; i++) if (network1->vec_aux_truly_feas_hid2[i] == 1) network1->truly_feas_hid2++;

    
    network1->inputs_nonfeas = 0;
    input1_nc=1; input2_nc=1; input3_nc=1; //nc=non-correct
    for (i=0;i<INPUT;i++)
        for(j=0;j<network1->hid1;j++)
        {
            if (i==0) if(network1->inhid1[i*network1->hid1+j] != 0 && network1->vec_aux_truly_feas_hid1[j]==1) input1_nc = 0;
            if (i==1) if(network1->inhid1[i*network1->hid1+j] != 0 && network1->vec_aux_truly_feas_hid1[j]==1) input2_nc = 0;
            if (i==2) if(network1->inhid1[i*network1->hid1+j] != 0 && network1->vec_aux_truly_feas_hid1[j]==1) input3_nc = 0;
        }
    
    if(input1_nc == 1 || input2_nc == 1 || input3_nc == 1) network1->inputs_nonfeas = 1;
    free(non_valid_neurons);
};

void CheckFeasibilityPTP(int twotoinput, network **network1){
    
    int i,j,k;
    int *non_valid_neurons; //1=non-feasible. 0 = feasible.
    non_valid_neurons = (int*)malloc(sizeof(int)*LAYERS*MAXSIZE);
    for (i=0;i<LAYERS*MAXSIZE;i++)
        non_valid_neurons[i]=0;
    
    int input1_nc,input2_nc,input3_nc;
    
    ValidNeurons(0, INPUT, (*network1)->hid1, (*network1)->hid2, (*network1)->inhid1, (*network1)->hid1hid2, non_valid_neurons);
    ValidNeurons(1, (*network1)->hid1, (*network1)->hid2, OUTPUT, (*network1)->hid1hid2, (*network1)->hid2out, non_valid_neurons);
    
    (*network1)->truly_feas_hid1 = 0;
    (*network1)->truly_feas_hid2 = 0;
    for (i=0; i<MAXSIZE; i++)
    {
        (*network1)->vec_aux_truly_feas_hid1[i] = 0;
        (*network1)->vec_aux_truly_feas_hid2[i] = 0;
    }
    
    (*network1)->feas_paths = 0;
    for (k=0;k<INPUT;k++)
        for (i=0;i<(*network1)->hid1;i++)
            if (non_valid_neurons[0*MAXSIZE+i] == 0 && (*network1)->inhid1[k*(*network1)->hid1+i]!=0)
                for (j=0;j<(*network1)->hid2;j++)
                    if (non_valid_neurons[1*MAXSIZE+j] == 0 && (*network1)->hid1hid2[i*(*network1)->hid2+j]!=0)
                    {
                        (*network1)->vec_aux_truly_feas_hid1[i] = 1;
                        (*network1)->vec_aux_truly_feas_hid2[j] = 1;
                        (*network1)->vec_feas_paths[(*network1)->feas_paths].input = k;
                        (*network1)->vec_feas_paths[(*network1)->feas_paths].hid1 = i;
                        (*network1)->vec_feas_paths[(*network1)->feas_paths].hid2 = j;
                        (*network1)->vec_feas_paths[(*network1)->feas_paths].output = 0; //valid only if one Output unit
                        (*network1)->feas_paths++;
                    }
    
    for (i=0; i<(*network1)->hid1; i++) if ((*network1)->vec_aux_truly_feas_hid1[i] == 1) (*network1)->truly_feas_hid1++;
    for (i=0; i<(*network1)->hid2; i++) if ((*network1)->vec_aux_truly_feas_hid2[i] == 1) (*network1)->truly_feas_hid2++;
    
    
    (*network1)->inputs_nonfeas = 0;
    input1_nc=1; input2_nc=1; input3_nc=1; //nc=non-correct
    
    for (i=0;i<INPUT;i++)
        for(j=0;j<(*network1)->hid1;j++)
        {
            if (i==0) if((*network1)->inhid1[i*(*network1)->hid1+j] != 0 && (*network1)->vec_aux_truly_feas_hid1[j]==1) input1_nc = 0;
            if (i==1) if((*network1)->inhid1[i*(*network1)->hid1+j] != 0 && (*network1)->vec_aux_truly_feas_hid1[j]==1) input2_nc = 0;
            if (i==2) if((*network1)->inhid1[i*(*network1)->hid1+j] != 0 && (*network1)->vec_aux_truly_feas_hid1[j]==1) input3_nc = 0;
        }
    
    if(input1_nc == 1 || input2_nc == 1 || input3_nc == 1) (*network1)->inputs_nonfeas = 1;
    free(non_valid_neurons);
};

