#ifndef Outputs_h
#define Outputs_h

#include <stdio.h>
#include "Structs.h"

void NetworksDOT (FILE *file_netw, network *networks, int *dominanceScore, int g, double *target1, double *target2, double *target3, double *Ff, const char *argvCOND);
void TargetsData (FILE *pareto_data, FILE *nonpareto_data, int g, int *dominanceScore, double *target1, double *target2, double *target3, double *Ff, const char *argvCOND);
void ConnectivityNetworksStorage(FILE *connectivity, network *networks, int g, int *domscore,const char *argvCOND);

#endif 
