
#ifndef GENALG_H
#define GENALG_H

#include "params.h"
#include "casim.h"


typedef unsigned int UINT;
typedef int BOOL;

// definice typu chromozomu GA - UPRAVTE SI DLE POTREBY
typedef struct {
    UINT chromosome[RULES*RULEN];  // vlastni napln chromozomu
    UINT fitness;   // fitness daneho jedince
    BOOL evaluate;  // zda je treba znovu vyhodnotit fitness
} GA_chromosome;

// prototypy funkci pro GA
UINT urandom(UINT low, UINT high);

void initialize(GA_chromosome *g);
void crossover(GA_chromosome *parent1, GA_chromosome *parent2, 
                GA_chromosome *child1, GA_chromosome *child2);
BOOL mutator(GA_chromosome *genome, UINT _pmut);
UINT fitness(GA_chromosome *genome, CA & ca);
void gprint(GA_chromosome *genome);
void evolve(CA & ca);
BOOL stop();

#endif
