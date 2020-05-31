#ifndef MutationsCrossover_h
#define MutationsCrossover_h

#include <stdio.h>
#include "Structs.h"

void Crossover (int *matrix1, int i1, int j1, int *matrix2, int i2, int j2);
void NeuronsMutation(network **network1, double In, double Dn);
void NeuronsMut(int *d_neurons, int *e_neurons, int *neurons_layer, int *new_neurons_layer, int *protected_neurons, double Dn, double In);
void NeuronsMutMatricesRestablish(int *weights_matrix, int layer_from, int layer_to, int *d_neurons, int *e_neurons, int *neurons_layer, int *new_neurons_layer);
void NeuronsMutActThreshRestablish(int layer, double *act_thresh_vec, int *neurons_layer, int *new_neurons_layer, int *d_neurons, int *e_neurons);
void ConnectionsMutation(int *weights_matrix, int n_from, int n_to, int nprotected_from, int nprotected_to ,double Cc, double Ec, double Ca);
void ActivationThresholdMutation(double *thresh_vec, int size_vec, double Ctheta, double tau);
void RegenerationRateMutation(double *reg_rate, double tau2, double Creg);

#endif
