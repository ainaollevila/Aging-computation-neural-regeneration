#include "EA.h"
#include "Constants.h"
#include "MutationsCrossover.h"
#include "Feasibility.h"
#include "StructsRelated.h"
#include <stdlib.h>

void computeDominanceScore(int nSolutions, int nTarget, double *tFunctions, int *dominanceScore){
    int i,j,k;
    for (i=0;i<nSolutions;i++)
        dominanceScore[i] = 0;
    bool fDominates_i, fDominates_j;
    for (i=0;i<nSolutions;i++)
    {
        for(j=i+1;j<nSolutions;j++)
        {
            fDominates_i = true;
            fDominates_j = true;
            for (k=0;k<nTarget;k++)
            {
                if (tFunctions[j*nTarget+k]<tFunctions[i*nTarget+k]) fDominates_i=false;
                if (tFunctions[i*nTarget+k]<tFunctions[j*nTarget+k]) fDominates_j=false;
            }
            if (fDominates_i || fDominates_j)
            {
                if (fDominates_i) dominanceScore[j]++;
                else if (fDominates_j) dominanceScore[i]++;
            }
        }
    }
};

void MatingNetworks(bool crossover, network *netw1, network *netw2, network *aux1, network *aux2, int twotoinput, double Cc, double Ec, double Ca, double In, double Dn, double Ctheta, double tau, double tau2, double Creg){
    
    int rand4;
    double rand1;
    rand1= (double) rand()/RAND_MAX;
    
    CheckFeasibilityPTP(twotoinput, &netw1);
    CheckFeasibilityPTP(twotoinput, &netw2);
    
    if (crossover && rand1<0.3)
    {
        DeepCopyStructPTP(&netw1, &aux1);
        DeepCopyStructPTP(&netw2, &aux2);
        
        rand4 = rand() % 3;
        if (rand4 == 0)
            Crossover(aux1->inhid1, INPUT, aux1->hid1, aux2->inhid1, INPUT, aux2->hid1);
        if (rand4 == 1)
            Crossover(aux1->hid1hid2, aux1->hid1, aux1->hid2, aux2->hid1hid2, aux2->hid1, aux2->hid2);
        if (rand4 == 2)
            Crossover(aux1->hid2out, aux1->hid2, OUTPUT, aux2->hid2out, aux2->hid2, OUTPUT);
        
        
        CheckFeasibilityPTP(twotoinput, &aux1);
        CheckFeasibilityPTP(twotoinput, &aux2);
        
        if (aux1->feas_paths >=1)
            DeepCopyStructPTP(&aux1, &netw1);
        if (aux2->feas_paths >=1)
            DeepCopyStructPTP(&aux2, &netw2);
        
        
        CheckFeasibilityPTP(twotoinput, &netw1);
        CheckFeasibilityPTP(twotoinput, &netw2);
    }
    
    
    netw1->protected_path = rand() % netw1->feas_paths;
    netw2->protected_path = rand() % netw2->feas_paths;
    
    
    ConnectionsMutation(netw1->inhid1, INPUT, netw1->hid1, netw1->vec_feas_paths[netw1->protected_path].input, netw1->vec_feas_paths[netw1->protected_path].hid1, Cc, Ec, Ca);
    ConnectionsMutation(netw1->hid1hid2, netw1->hid1, netw1->hid2, netw1->vec_feas_paths[netw1->protected_path].hid1, netw1->vec_feas_paths[netw1->protected_path].hid2, Cc, Ec, Ca);
    ConnectionsMutation(netw1->hid2out, netw1->hid2, OUTPUT, netw1->vec_feas_paths[netw1->protected_path].hid2, 0, Cc, Ec, Ca);
    
    ConnectionsMutation(netw2->inhid1, INPUT, netw2->hid1, netw2->vec_feas_paths[netw2->protected_path].input, netw2->vec_feas_paths[netw2->protected_path].hid1, Cc, Ec, Ca);
    ConnectionsMutation(netw2->hid1hid2, netw2->hid1, netw2->hid2, netw2->vec_feas_paths[netw2->protected_path].hid1, netw2->vec_feas_paths[netw2->protected_path].hid2, Cc, Ec, Ca);
    ConnectionsMutation(netw2->hid2out, netw2->hid2, OUTPUT, netw2->vec_feas_paths[netw2->protected_path].hid2, 0, Cc, Ec, Ca);
    
    
    NeuronsMutation(&netw1, In, Dn);
    NeuronsMutation(&netw2, In, Dn);
    
    CheckFeasibilityPTP(twotoinput, &netw1);
    CheckFeasibilityPTP(twotoinput, &netw2);
    
    InputsFeasibility(netw1->hid1, netw1->inhid1, netw1->vec_aux_truly_feas_hid1, netw1->truly_feas_hid1);
    InputsFeasibility(netw2->hid1, netw2->inhid1, netw2->vec_aux_truly_feas_hid1, netw2->truly_feas_hid1);
    
    
    ActivationThresholdMutation(netw1->theta_hid1, netw1->hid1, Ctheta, tau);
    ActivationThresholdMutation(netw1->theta_hid2, netw1->hid2, Ctheta, tau);
    ActivationThresholdMutation(netw1->theta_output, OUTPUT, Ctheta, tau);
    
    ActivationThresholdMutation(netw2->theta_hid1, netw2->hid1, Ctheta, tau);
    ActivationThresholdMutation(netw2->theta_hid2, netw2->hid2, Ctheta, tau);
    ActivationThresholdMutation(netw2->theta_output, OUTPUT, Ctheta, tau);
    
    RegenerationRateMutation(&netw1->regeneration_rate, tau2, Creg);
    RegenerationRateMutation(&netw2->regeneration_rate, tau2, Creg);
};

