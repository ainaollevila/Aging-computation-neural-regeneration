#ifndef StructsRelated_h
#define StructsRelated_h

#include <stdio.h>
#include "Structs.h"

void DeepCopyStruct (network *network1, network *network2);
void DeepCopyStructPTP(network **network1, network **network2);
void AllocateMemoryStruct(network *netw1, int twotoinput);
void InitStructNetworks(network *netw1, int twotoinput);
void InitPop(network *netw1);

#endif 
