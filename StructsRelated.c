#include "StructsRelated.h"
#include <stdlib.h>
#include "Constants.h"
#include "GlobalVariables.h"

void DeepCopyStruct (network *network1, network *network2){
    network2->hid1 = network1->hid1;
    network2->hid2 = network1->hid2;
    network2->regeneration_rate = network1->regeneration_rate;
    int i;
    for (i=0;i<network2->hid1*INPUT;i++) network2->inhid1[i] = network1->inhid1[i];
    for (i=0;i<network2->hid1*network1->hid2;i++) network2->hid1hid2[i] = network1->hid1hid2[i];
    for (i=0;i<network2->hid2*OUTPUT;i++) network2->hid2out[i] = network1->hid2out[i];
    for (i=0;i<network2->hid1;i++) network2->theta_hid1[i] = network1->theta_hid1[i];
    for (i=0;i<network2->hid2;i++) network2->theta_hid2[i] = network1->theta_hid2[i];
    for (i=0;i<OUTPUT;i++) network2->theta_output[i] = network1->theta_output[i];
};

void DeepCopyStructPTP(network **network1, network **network2){
    (*network2)->hid1 = (*network1)->hid1;
    (*network2)->hid2 = (*network1)->hid2;
    (*network2)->regeneration_rate = (*network1)->regeneration_rate;
    int i;
    for (i=0;i<(*network2)->hid1*INPUT;i++) (*network2)->inhid1[i] = (*network1)->inhid1[i];
    for (i=0;i<(*network2)->hid1*(*network1)->hid2;i++) (*network2)->hid1hid2[i] = (*network1)->hid1hid2[i];
    for (i=0;i<(*network2)->hid2*OUTPUT;i++) (*network2)->hid2out[i] = (*network1)->hid2out[i];
    for (i=0;i<(*network2)->hid1;i++) (*network2)->theta_hid1[i] = (*network1)->theta_hid1[i];
    for (i=0;i<(*network2)->hid2;i++) (*network2)->theta_hid2[i] = (*network1)->theta_hid2[i];
    for (i=0;i<OUTPUT;i++) (*network2)->theta_output[i] = (*network1)->theta_output[i];
};

void AllocateMemoryStruct(network *netw1, int twotoinput){
    
    netw1->inhid1 = (int*)malloc(sizeof(int)*INPUT*MAXSIZE);
    netw1->hid1hid2 = (int*)malloc(sizeof(int)*MAXSIZE*MAXSIZE);
    netw1->hid2out = (int*)malloc(sizeof(int)*MAXSIZE*OUTPUT);
    
    netw1->theta_hid1 = (double*)malloc(sizeof(double)*MAXSIZE);
    netw1->theta_hid2 = (double*)malloc(sizeof(double)*MAXSIZE);
    netw1->theta_output = (double*)malloc(sizeof(double)*OUTPUT);
    
    netw1->vec_aux_truly_feas_hid1 = (int*)malloc(sizeof(int)*MAXSIZE);
    netw1->vec_aux_truly_feas_hid2 = (int*)malloc(sizeof(int)*MAXSIZE);
    
    
    netw1->vec_feas_paths = (feasible_path*)malloc(sizeof(feasible_path)*INPUT*MAXSIZE*MAXSIZE);
    
};

void InitStructNetworks(network *netw1, int twotoinput){
    int j;
    
    for (j=0;j<MAXSIZE*INPUT;j++) netw1->inhid1[j] = 0;
    for (j=0;j<MAXSIZE*MAXSIZE;j++) netw1->hid1hid2[j] = 0;
    for (j=0;j<MAXSIZE*OUTPUT;j++) netw1->hid2out[j] = 0;
    
    for (j=0;j<MAXSIZE;j++)
    {
        netw1->theta_hid1[j] = 0;
        netw1->theta_hid2[j] = 0;
    }
    for (j=0;j<OUTPUT;j++) netw1->theta_output[j] = 0;
    for (j=0;j<MAXSIZE;j++)
    {
        netw1->vec_aux_truly_feas_hid1[j] = 0;
        netw1->vec_aux_truly_feas_hid2[j] = 0;
    }
    
    for (j=0;j<INPUT*MAXSIZE*MAXSIZE;j++)
    {
        netw1->vec_feas_paths[j].input = 0;
        netw1->vec_feas_paths[j].hid1 = 0;
        netw1->vec_feas_paths[j].hid2 = 0;
        netw1->vec_feas_paths[j].output = 0;
    }

    netw1->protected_path = 0;
    netw1->hid1 = 0;
    netw1->hid2 = 0;
    netw1->truly_feas_hid1 = 0;
    netw1->truly_feas_hid2 = 0;
    netw1->feas_paths = 0;
    netw1->inputs_nonfeas = 0;
};

void InitPop(network *netw1){
    
    int i,j;
    double rand1, rand2;
    
    int hid1ini = MAXSIZE, hid2ini = MAXSIZE;
    netw1->hid1 = hid1ini;
    netw1->hid2 = hid2ini;
  
    double connection_prob = 1;
    
    //INPUT-HIDDEN1
    for (i=0;i<INPUT;i++)
        for(j=0;j<hid1ini;j++)
        {
            rand1 = (double)rand() / RAND_MAX;
            rand2 = (double)rand() / RAND_MAX;
            if (rand1<connection_prob) { if (rand2<0.5) netw1->inhid1[i*hid1ini+j] = 1; else netw1->inhid1[i*hid1ini+j] = -1;}
        }
        
    //HIDDEN1-HIDDEN2
    for (i=0;i<hid1ini;i++)
        for(j=0;j<hid2ini;j++)
        {
            rand1 = (double)rand() / RAND_MAX;
            rand2 = (double)rand() / RAND_MAX;
            if (rand1<connection_prob) { if (rand2<0.5) netw1->hid1hid2[i*hid2ini+j] = 1; else netw1->hid1hid2[i*hid2ini+j] = -1;}
        }
        
    //HIDDEN2-OUTPUT
    for (i=0;i<hid2ini;i++)
        for(j=0;j<OUTPUT;j++)
        {
            rand1 = (double)rand() / RAND_MAX;
            rand2 = (double)rand() / RAND_MAX;
            if (rand1<connection_prob) { if (rand2<0.5) netw1->hid2out[i*OUTPUT+j] = 1; else netw1->hid2out[i*OUTPUT+j] = -1;}
        }
    
    //ACTIVATION THRESHOLD
    for (i=0;i<hid1ini;i++)
    {
        rand1 = (double) rand() /RAND_MAX;
        if (rand1>=0.01)
            netw1->theta_hid1[i] = rand1;
        else
        {
            while (rand1<0.01)
                rand1 = (double) rand() /RAND_MAX;
            netw1->theta_hid1[i] = rand1;
        }
    }
    for (i=0;i<hid2ini;i++)
    {
        rand1 = (double) rand()/RAND_MAX;
        if (rand1>=0.01)
            netw1->theta_hid2[i] = rand1;
        else
        {
            while (rand1<0.01)
                rand1 = (double) rand() /RAND_MAX;
            netw1->theta_hid2[i] = rand1;
        }
    }
    for (i=0;i<OUTPUT;i++)
    {
        rand1 = (double) rand() /RAND_MAX;
        if (rand1>=0.01)
            netw1->theta_output[i] = rand1;
        else
        {
            while (rand1<0.01)
                rand1 = (double) rand() /RAND_MAX;
            netw1->theta_output[i] = rand1;
        }
    }

    //REGENERATION RATES
    rand1 = (double)rand()/RAND_MAX;
    if (rand1>=minimum_regeneration_rate)
        netw1->regeneration_rate = rand1;
    else
    {
        while (rand1<minimum_regeneration_rate)
            rand1 = (double)rand()/RAND_MAX;
        netw1->regeneration_rate = rand1;
    }
    
};
