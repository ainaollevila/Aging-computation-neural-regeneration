#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <limits.h>
#include <time.h>
#include <assert.h>
#include "Constants.h"
#include "GlobalVariables.h"
#include "Structs.h"
#include "StructsRelated.h"
#include "TargetsComputation.h"
#include "Feasibility.h"
#include "MutationsCrossover.h"
#include "Outputs.h"
#include "GeneralUsageFunctions.h"
#include "Aging.h"
#include "EA.h"

int main(int argc, const char * argv[]) {
    
    if (argc < 6)
    {
        printf("ERROR. Wrong number of argv. Use: rho_min delta tau BooleanFunction(0 for Mux, 1 for Maj) replicaID\n");
        exit(-1);
    }
    int BooleanFunction=0;
    char *Bool[2]={"Mux","Maj"};

    minimum_regeneration_rate = atof(argv[1]);
    damage_rate = atof(argv[2]);
    TIME_AGEING = atoi(argv[3]);
    BooleanFunction = atoi(argv[4]);
    REPLICA = atoi(argv[5]);
    
    srand(REPLICA*10000);


    printf("Running EA of %s function: rho_min = %lf; damage rate =%lf; lifetime=%d; REPLICA: %03d\n", Bool[BooleanFunction], minimum_regeneration_rate, damage_rate, TIME_AGEING, REPLICA);
    

    bool Crossover=1;
    bool ElitistOperator = 1;
    
    //Mutational probabilites
    double Ec, Cc, In, Dn, Ca, Ctheta;
    Ec = 0.1; //Eliminate existing connection
    Cc = 0.1; //Originate new connection
    In = 0.1; //Eliminate a node
    Dn = 0.1; //Duplicate a node and its connections
    Ca = 0.1; //Change the weight sign
    Ctheta = 0.1; //Change the activation threshold
    double tau = 0.01; //used to weight mutational change of activation threshold
    double Creg=1; //Change the regeneration rate
    double tau2=0.1;//used to weight mutational change of regeneration rate
    
    //Output files definition
    FILE *data=NULL,*file_netw=NULL,*file_all=NULL,*conn=NULL;
    
    
    //Variables definition
    int i=0,k=0,g=0,in=0,a=0,t=0,sel=0;
    int twotoinput = pow(2, INPUT);
    
    //Mating variables
    int k1, k2;
    int newind;
    int pareto_solutions;
    int remain;
    
    //Truth table definition
    int *input_comb, *exp_out;
    input_comb = (int*)malloc(twotoinput*INPUT*sizeof(int));
    exp_out = (int*)malloc(twotoinput*sizeof(int));
    
    
    input_comb[0]= 0;
    input_comb[1]= 0;
    input_comb[2]=0;
    
    input_comb[3]=0;
    input_comb[4]=0;
    input_comb[5]=1;
    
    input_comb[6]=0;
    input_comb[7]=1;
    input_comb[8]=0;
    
    input_comb[9]= 0;
    input_comb[10]= 1;
    input_comb[11]=1;
    
    input_comb[12]=1;
    input_comb[13]=0;
    input_comb[14]=0;
    
    input_comb[15]=1;
    input_comb[16]=0;
    input_comb[17]=1;
    
    input_comb[18]= 1;
    input_comb[19]= 1;
    input_comb[20]=0;
    
    input_comb[21]=1;
    input_comb[22]=1;
    input_comb[23]=1;
    

    if (BooleanFunction == 0)
    {
        //Multiplexer (Mux)
        exp_out[0] = 0;
        exp_out[1]= 0;
        exp_out[2] = 1;
        exp_out[3] = 1;
        exp_out[4]= 0;
        exp_out[5] = 1;
        exp_out[6] = 0;
        exp_out[7]= 1;
    }
    else if (BooleanFunction == 1)
    {
        //Majority (Maj)
        exp_out[0] = 0;
        exp_out[1]= 0;
        exp_out[2] = 0;
        exp_out[3] = 1;
        exp_out[4]= 0;
        exp_out[5] = 1;
        exp_out[6] = 1;
        exp_out[7]= 1;
    }

    //target1: computation error;; target2: number of connections;; target3:network regeneration rate
    double *target1, *target2, *target3;
    target1 = (double*)malloc(sizeof(double)*SIZEPOP);
    target2 = (double*)malloc(sizeof(double)*SIZEPOP);
    target3 = (double*)malloc(sizeof(double)*SIZEPOP);
    
    double *Ff; //Feasibility function
    Ff = (double*)malloc(sizeof(double)*SIZEPOP);

    double aux_t1=0, aux_Ff=0;

    for (i=0;i<SIZEPOP;i++)
    {
        target1[i]=0;
        target2[i]=0;
        target3[i]=0;
        Ff[i]=0;
    }
    
    int *dominanceScore;
    dominanceScore = (int*)malloc(sizeof(int)*SIZEPOP);
    for (i=0;i<SIZEPOP;i++) dominanceScore[i] = 0;
    
    double *tFunctions;
    tFunctions = (double*)malloc(sizeof(double)*nTarget*SIZEPOP);
    for (i=0;i<nTarget*SIZEPOP;i++) tFunctions[i] = 0;
    
    int *DomScoreOrdered,*DomScorePosition;
    DomScoreOrdered = (int*)malloc(sizeof(int)*SIZEPOP);
    DomScorePosition = (int*)malloc(sizeof(int)*SIZEPOP);
    for (i=0;i<SIZEPOP;i++) {DomScoreOrdered[i]=0; DomScorePosition[i]=0;}
    
    //Computation-related
    int *obs_out_ageing; //Observed output during AGING
    double *int_hid1, *int_hid2, *int_out; //Internal States
    int *ext_in, *ext_hid1, *ext_hid2, *ext_out; //External States
    
    obs_out_ageing = (int*)malloc(sizeof(int)*twotoinput*TIME_AGEING);
    
    int_hid1 = (double*)malloc(sizeof(double)*MAXSIZE);
    int_hid2 = (double*)malloc(sizeof(double)*MAXSIZE);
    int_out = (double*)malloc(sizeof(double)*OUTPUT);
    
    ext_in = (int*)malloc(sizeof(int)*INPUT);
    ext_hid1 = (int*)malloc(sizeof(int)*MAXSIZE);
    ext_hid2 = (int*)malloc(sizeof(int)*MAXSIZE);
    ext_out = (int*)malloc(sizeof(int)*OUTPUT);
    
    for (i=0;i<twotoinput*TIME_AGEING;i++) obs_out_ageing[i] = 0;

    for (i=0;i<MAXSIZE;i++)
    {
        int_hid1[i] = 0;
        int_hid2[i] = 0;
    }
    for (i=0;i<OUTPUT;i++) int_out[i] = 0;
    
    for (i=0;i<INPUT;i++) ext_in[i] = 0;
    for (i=0;i<MAXSIZE;i++)
    {
        ext_hid1[i] = 0;
        ext_hid2[i] = 0;
    }
    for (i=0;i<OUTPUT;i++) ext_out[i] = 0;
    
    
    //"Memory vectors", to restablish the former connections' weights once regenerated.
    int *inhid1_damaged, *hid1hid2_damaged, *hid2out_damaged;
    inhid1_damaged = (int*)malloc(sizeof(int)*INPUT*MAXSIZE);
    hid1hid2_damaged = (int*)malloc(sizeof(int)*MAXSIZE*MAXSIZE);
    hid2out_damaged = (int*)malloc(sizeof(int)*OUTPUT*MAXSIZE);
    for (i=0;i<INPUT*MAXSIZE;i++) inhid1_damaged[i] = 0;
    for (i=0;i<MAXSIZE*MAXSIZE;i++) hid1hid2_damaged[i] = 0;
    for (i=0;i<MAXSIZE*OUTPUT;i++) hid2out_damaged[i] = 0;
    
    //Used to keep a copy of the network state at t=0 (before aging)
    int *inhid1_aux,*hid1hid2_aux, *hid2out_aux;
    inhid1_aux = (int*)malloc(sizeof(int)*INPUT*MAXSIZE);
    hid1hid2_aux = (int*)malloc(sizeof(int)*MAXSIZE*MAXSIZE);
    hid2out_aux = (int*)malloc(sizeof(int)*OUTPUT*MAXSIZE);
    for (i=0;i<INPUT*MAXSIZE;i++) inhid1_aux[i] = 0;
    for (i=0;i<MAXSIZE*MAXSIZE;i++) hid1hid2_aux[i] = 0;
    for (i=0;i<MAXSIZE*OUTPUT;i++) hid2out_aux[i] = 0;
    
    //Used to keep track of the units' "feasibility".
    int *aux_ageing_hid1_feas, *aux_ageing_hid2_feas;
    aux_ageing_hid1_feas = (int*)malloc(sizeof(int)*MAXSIZE);
    aux_ageing_hid2_feas = (int*)malloc(sizeof(int)*MAXSIZE);
    for (i=0;i<MAXSIZE;i++) {aux_ageing_hid1_feas[i]=0; aux_ageing_hid2_feas[i]=0;}
    
    int *aux_ageing_input_feas, *aux_ageing_output_feas;
    aux_ageing_input_feas = (int*)malloc(sizeof(int)*INPUT);
    aux_ageing_output_feas = (int*)malloc(sizeof(int)*OUTPUT);
    for (i=0;i<INPUT;i++) aux_ageing_input_feas[i] = 1;
    for (i=0;i<OUTPUT;i++) aux_ageing_output_feas[i] = 1;
    

    //Structs' definition
    network *networks, *next_networks, *temp, *aux, mating_network;
    networks = (network*)malloc(sizeof(network)*SIZEPOP);
    next_networks = (network*)malloc(sizeof(network)*SIZEPOP);
    aux = (network*)malloc(sizeof(network)*2);
    
    //Allocating memory
    for (i=0;i<SIZEPOP;i++) AllocateMemoryStruct(&networks[i],twotoinput);
    for (i=0;i<SIZEPOP;i++) AllocateMemoryStruct(&next_networks[i], twotoinput);
    for (i=0;i<2;i++) AllocateMemoryStruct(&aux[i],twotoinput);
    AllocateMemoryStruct(&mating_network,twotoinput);
    
    //Initalisation of structs' variables
    for (i=0;i<SIZEPOP;i++) InitStructNetworks(&networks[i],twotoinput);
    for (i=0;i<SIZEPOP;i++) InitStructNetworks(&next_networks[i],twotoinput);
    for (i=0;i<2;i++) InitStructNetworks(&aux[i],twotoinput);
    InitStructNetworks(&mating_network,twotoinput);
    
    //Initialisation of networks' population undergoing evolution
    for (i=0;i<SIZEPOP;i++) InitPop(&networks[i]);
    
    while(g<GENERATIONS)
    {
        newind=0;
        if (g%100==0)
            printf("\nGeneration %d; REPLICA %03d; lifetime = %d; damagerate = %lf; rho_min %lf;\n",g,REPLICA, TIME_AGEING, damage_rate, minimum_regeneration_rate);
        
        for (k=0;k<SIZEPOP;k++)
        {
            for (i=0;i<INPUT*MAXSIZE;i++)
                inhid1_aux[i] = 0;
            for (i=0;i<MAXSIZE*MAXSIZE;i++)
                hid1hid2_aux[i] = 0;
            for (i=0;i<OUTPUT*MAXSIZE;i++)
                hid2out_aux[i] = 0;
            
            CheckFeasibility(twotoinput, &networks[k]);
            
            //Computing target 2
            target2[k] = 0;
            target2[k] += CountConnections(INPUT, networks[k].hid1, aux_ageing_input_feas, networks[k].vec_aux_truly_feas_hid1, networks[k].inhid1);
            target2[k] += CountConnections(networks[k].hid1, networks[k].hid2, networks[k].vec_aux_truly_feas_hid1, networks[k].vec_aux_truly_feas_hid2, networks[k].hid1hid2);
            target2[k] += CountConnections(networks[k].hid2, OUTPUT, networks[k].vec_aux_truly_feas_hid2, aux_ageing_output_feas, networks[k].hid2out);
            
            
            //Storing target 3
            target3[k] = networks[k].regeneration_rate;

            
            //Network copied to auxiliary arrays at t=0
            for (i=0;i<INPUT*networks[k].hid1;i++)
                inhid1_aux[i] = networks[k].inhid1[i];
            for (i=0;i<networks[k].hid1*networks[k].hid2;i++)
                hid1hid2_aux[i] = networks[k].hid1hid2[i];
            for (i=0;i<OUTPUT*networks[k].hid2;i++)
                hid2out_aux[i] = networks[k].hid2out[i];
            
            //Assess computation error at t=0 (before aging)
            for (in=0;in<twotoinput;in++)
            {
                //All units' internal state set to zero
                for (i=0;i<networks[k].hid1;i++) int_hid1[i] = 0;
                for (i=0;i<networks[k].hid2;i++) int_hid2[i] = 0;
                for (i=0;i<OUTPUT;i++) int_out[i] = 0;
                
                //Setting external states of input neurons
                ext_in[0] = input_comb[in*INPUT+0];
                ext_in[1] = input_comb[in*INPUT+1];
                ext_in[2] = input_comb[in*INPUT+2];
                
                //Computing external states of the remaining network units
                StepFunction(INPUT, networks[k].hid1, int_hid1, ext_hid1, ext_in, networks[k].inhid1, networks[k].theta_hid1);
                StepFunction(networks[k].hid1, networks[k].hid2, int_hid2, ext_hid2, ext_hid1, networks[k].hid1hid2, networks[k].theta_hid2);
                StepFunction(networks[k].hid2, OUTPUT, int_out, ext_out, ext_hid2, networks[k].hid2out, networks[k].theta_output);
                
                //Store the Observed Output
                obs_out_ageing[in] = ext_out[0];
            }
            
            for(i=0;i<MAXSIZE;i++)
            {
                aux_ageing_hid1_feas[i] = networks[k].vec_aux_truly_feas_hid1[i];
                aux_ageing_hid2_feas[i] = networks[k].vec_aux_truly_feas_hid2[i];
            }
            
            aux_t1=0;
            aux_Ff=0;
            
            for (a=0;a<AVERAGE_OVER;a++)
            {
                Ff[k]=0;
                
                for (i=0;i<INPUT*MAXSIZE;i++)
                    inhid1_damaged[i] = 0;
                for (i=0;i<MAXSIZE*MAXSIZE;i++)
                    hid1hid2_damaged[i] = 0;
                for (i=0;i<OUTPUT*MAXSIZE;i++)
                    hid2out_damaged[i] = 0;
                
                
                for (t=1;t<TIME_AGEING;t++)
                {
                    
                    Damage(INPUT, networks[k].hid1, networks[k].inhid1, inhid1_damaged, aux_ageing_input_feas, aux_ageing_hid1_feas);
                    Damage(networks[k].hid1, networks[k].hid2, networks[k].hid1hid2, hid1hid2_damaged, aux_ageing_hid1_feas, aux_ageing_hid2_feas);
                    Damage(networks[k].hid2, OUTPUT, networks[k].hid2out, hid2out_damaged, aux_ageing_hid2_feas, aux_ageing_output_feas);
                    
                    Regeneration(INPUT, networks[k].hid1, networks[k].inhid1, inhid1_damaged, networks[k].regeneration_rate);
                    Regeneration(networks[k].hid1, networks[k].hid2, networks[k].hid1hid2, hid1hid2_damaged, networks[k].regeneration_rate);
                    Regeneration(networks[k].hid2, OUTPUT, networks[k].hid2out, hid2out_damaged, networks[k].regeneration_rate);
                    
                    
                    //Assess computation error at t=0 (before aging)
                    for (in=0;in<twotoinput;in++)
                    {
                        //All units' internal state set to zero
                        for (i=0;i<networks[k].hid1;i++) int_hid1[i] = 0;
                        for (i=0;i<networks[k].hid2;i++) int_hid2[i] = 0;
                        for (i=0;i<OUTPUT;i++) int_out[i] = 0;
                        
                        //Setting external states of input neurons
                        ext_in[0] = input_comb[in*INPUT+0];
                        ext_in[1] = input_comb[in*INPUT+1];
                        ext_in[2] = input_comb[in*INPUT+2];
                        
                        //Computing external states of the remaining network units
                        StepFunction(INPUT, networks[k].hid1, int_hid1, ext_hid1, ext_in, networks[k].inhid1, networks[k].theta_hid1);
                        StepFunction(networks[k].hid1, networks[k].hid2, int_hid2, ext_hid2, ext_hid1, networks[k].hid1hid2, networks[k].theta_hid2);
                        StepFunction(networks[k].hid2, OUTPUT, int_out, ext_out, ext_hid2, networks[k].hid2out, networks[k].theta_output);
                        
                        //Store the Observed Output
                        obs_out_ageing[t*twotoinput+in] = ext_out[0];
                    }
                    
                    CheckFeasibility(twotoinput, &networks[k]);
                    if (networks[k].feas_paths == 0 || networks[k].inputs_nonfeas == 1)
                        Ff[k] += 1;
                }
                target1[k] = ExpectedObservedOutputAccuracy(twotoinput, TIME_AGEING, obs_out_ageing, exp_out);
                
                Ff[k] /= (double) TIME_AGEING;
                aux_Ff += Ff[k];
                aux_t1 += target1[k];
                
                //Network topology restablished to initial state (before undergoing aging)
                for (i=0;i<INPUT*networks[k].hid1;i++)
                    networks[k].inhid1[i]=inhid1_aux[i];
                for (i=0;i<networks[k].hid1*networks[k].hid2;i++)
                    networks[k].hid1hid2[i] = hid1hid2_aux[i];
                for (i=0;i<OUTPUT*networks[k].hid2;i++)
                    networks[k].hid2out[i]=hid2out_aux[i];
                CheckFeasibility(twotoinput, &networks[k]);
            }
            target1[k] = aux_t1/AVERAGE_OVER;
            Ff[k] = aux_Ff/AVERAGE_OVER;
        }
        
        for (i=0;i<SIZEPOP;i++)
        {
            tFunctions[i*nTarget] = target1[i];
            tFunctions[i*nTarget+1] = target2[i];
            tFunctions[i*nTarget+2] = target3[i];
        }
        
        computeDominanceScore(SIZEPOP, nTarget, tFunctions, dominanceScore);
        
        if (g%500==0 || (g<500 && g%50==0))
            TargetsData(data, file_all, g, dominanceScore, target1, target2, target3, Ff, Bool[BooleanFunction]);
            
        
        if (g==GENERATIONS-1)
        {
            TargetsData(data, file_all, g, dominanceScore, target1, target2, target3, Ff, Bool[BooleanFunction]);
            NetworksDOT(file_netw, networks, dominanceScore, g, target1, target2, target3, Ff, Bool[BooleanFunction]);
            ConnectivityNetworksStorage(conn,networks,g,dominanceScore,Bool[BooleanFunction]);
        }
        
        //Breeding operations to create next generation
        
        //Elitist operator: Pareto solutions are copied to the next generation
        pareto_solutions=0;
        if (ElitistOperator)
        {
            for (i=0;i<SIZEPOP;i++)
            {
                if (dominanceScore[i]==0)
                {
                    DeepCopyStruct(&networks[i], &next_networks[pareto_solutions]);
                    pareto_solutions++;
                }
            }
        }
        
        //Breeding operations to complete next generation
        remain = (int)(SIZEPOP - pareto_solutions)/2;
        newind = pareto_solutions;
        
        SortSmalltoBig(DomScorePosition, dominanceScore, DomScoreOrdered, SIZEPOP);
        
        for (sel=0;sel<remain;sel++)
        {
            //Only the best half of networks' population can be selected (randomly) for breeding
            k1 = rand() % (SIZEPOP/2);
            k2 = rand() % (SIZEPOP/2);
            DeepCopyStruct(&networks[DomScorePosition[k1]],&next_networks[newind]);
            DeepCopyStruct(&networks[DomScorePosition[k2]],&next_networks[newind+1]);
            
            
            MatingNetworks(Crossover, &next_networks[newind], &next_networks[newind+1], &aux[0], &aux[1], twotoinput, Cc, Ec, Ca, In, Dn, Ctheta, tau, tau2, Creg);
            
            newind += 2; 
        }
        

        if ((SIZEPOP - pareto_solutions)%2 !=0)
        {
            k1 = rand() % (SIZEPOP/2);
            k2 = rand() % (SIZEPOP/2);
            DeepCopyStruct(&networks[DomScorePosition[k1]],&next_networks[newind]);
            DeepCopyStruct(&networks[DomScorePosition[k2]],&mating_network);
            
            
            MatingNetworks(Crossover, &next_networks[newind], &mating_network, &aux[0], &aux[1], twotoinput, Cc, Ec, Ca, In, Dn, Ctheta, tau, tau2, Creg);
        }
        
        //Pointers' swap
        temp = networks;
        networks = next_networks;
        next_networks = temp;
        g++;
    }
    

    free(networks);
    free(next_networks);
    
    free(aux);
    
    free(input_comb);
    free(exp_out);

    free(target1);
    free(target2);
    free(target3);
    free(Ff);
    
    free(dominanceScore);
    free(tFunctions);
    free(DomScoreOrdered);
    free(DomScorePosition);
    
    free(inhid1_damaged);
    free(hid1hid2_damaged);
    free(hid2out_damaged);

    free(inhid1_aux);
    free(hid1hid2_aux);
    free(hid2out_aux);
    
    free(aux_ageing_hid1_feas);
    free(aux_ageing_hid2_feas);
    free(aux_ageing_input_feas);
    free(aux_ageing_output_feas);
    
    return 0;
}



