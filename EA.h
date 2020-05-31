#ifndef EA_h
#define EA_h

#include <stdio.h>
#include <stdbool.h>
#include "Structs.h"

void computeDominanceScore(int nSolutions, int nTarget, double *tFunctions, int *dominanceScore);
void MatingNetworks(bool crossover, network *netw1, network *netw2, network *aux1, network *aux2, int twotoinput, double Cc, double Ec, double Ca, double In, double Dn, double Ctheta, double tau, double tau2, double Creg);
#endif 
