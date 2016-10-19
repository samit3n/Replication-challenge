//////////////////////////////////////////////////////////////
// VZOROVA IMPLEMENTACE JEDNODUCHEHO GENETICKEHO ALGORITMU  //
// JAZYK C, PREKLAD S OHLEDEM NA NORMU C99                  //
// (c) MICHAL BIDLO, 2011                                   //
// UVEDENY KOD SMI BYT POUZIT VYHRADNE PRO POTREBY PROJEKTU //
// V KURZU BIOLOGII INSPIROVANE POCITACE NA FIT VUT V BRNE  //
//////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "genalg.h"
#include <signal.h>
#include <unistd.h>


using namespace std;

// jednotky pro specifikaci pravdepodobnosti genetickych operatoru:
// 100 - procenta, 1000 - promile. NENI TREBA MENIT
const UINT unit = 100;
// delka chromozomu - s touto promennou pracuje GA
//
const UINT glength = RULES * RULEN;
//
// maximalni fitness - zde odpovida delce jedince, protoze se pouze snazime
// nalezt retezec se samymi 1 geny
//
//const UINT max_fit = CHLEN -1 

const UINT MIN_FIT = MFIT;

#ifdef DEBUG
const UINT generations = 0; // 0 - pocet generaci neni pri ladeni omezen
#else
const UINT generations = GENERATIONS;   // po tomto poctu je GA zastaven
#endif

// !!!! nas GA pracuje pouze se sudymi pocty jedincu v populaci !!!!
UINT _popsize = (POPSIZE & 1) ? POPSIZE + 1 : POPSIZE;

// ------------------- implementace genetickeho algoritmu ----------------------
// *****************************************************************************

GA_chromosome best;	// njelepsi dosud nalezene reseni
UINT best_ever; // fitness dosud nejlepsiho jedince

UINT generation;    // pocitadlo generaci
GA_chromosome *population;
GA_chromosome *next_population;
// pracovni populace - pouze sudy pocet jedincu !!!!
GA_chromosome pool1[(POPSIZE & 1) ? POPSIZE + 1 : POPSIZE];
GA_chromosome pool2[(POPSIZE & 1) ? POPSIZE + 1 : POPSIZE];


// evolucni cyklus SGA: nova populace o _popsize jedincich je tvorena krizenimi
// a mutaci turnajem vybranych jedincu predchazejici generace. 

void my_handler(int s){
    
    cout << "EA stopeed at " << generation << " generations by SIGINT" << endl;
    cout << "Best chrom" << endl;
    gprint(&best);
    cout << "Best fit" << best_ever << endl;
           
    exit(1);
	   
}

void evolve(CA & ca)
{
    //ctrl-c handl
    
    struct sigaction sigIntHandler;

   sigIntHandler.sa_handler = my_handler;
   sigemptyset(&sigIntHandler.sa_mask);
   sigIntHandler.sa_flags = 0;

   sigaction(SIGINT, &sigIntHandler, NULL);
   
   
   
    // inicializace promennych
    generation = 0;
    best.fitness = 0;
    best_ever = 0;
    GA_chromosome ind1_new, ind2_new;
    UINT _tour = (TOUR >= 2 ? TOUR : 2);
    UINT i1;
    
    // inicializace populace
    for (UINT i = 0; i < _popsize; i++)
    {
    	initialize(&pool1[i]);
    	pool1[i].evaluate = 1;
    }
    // evolucni cyklus
    do
    {
    	generation++;
    	if (generation & 1)
    	{
    		population = pool1;
    		next_population = pool2;
    	}
    	else
    	{
    		population = pool2;
    		next_population = pool1;
    	}
    	// ohodnoceni populace
    	for (UINT i = 0; i < _popsize; i++)
    	{
    		if (population[i].evaluate)
    		{
    			population[i].fitness = fitness(&population[i], ca);
    			if (population[i].fitness >= best.fitness)
    			    best = population[i];
    			population[i].evaluate = 0;
    		}
    	}
    	// elitizmus
    	next_population[0] = best;  // dosud nejlepsi nalezeny jedinec...
        GA_chromosome mutant = best; mutator(&mutant, unit);
    	next_population[1] = mutant;    // ...a mutant nejlepsiho
    	// tvorba nove populace
    	for (UINT i = 2; i < _popsize; i += 2)
    	{
    		GA_chromosome *ind1 = NULL, *ind2 = NULL;
    		// turnajovy vyber jedincu
    		for (UINT t = 0; t < _tour; t++)
    		{
    			i1 = urandom(0, _popsize - 1);
    			if (ind1 == NULL) ind1 = &population[i1];
    			else if (ind2 == NULL) ind2 = &population[i1];
    			else if (population[i1].fitness > ind1->fitness)
    				ind1 = &population[i1];
    			else if (population[i1].fitness > ind2->fitness)
    				ind2 = &population[i1];
    		}
    		// krizeni
    		if (urandom(0, unit) < PCROSS)
    		{
    			crossover(ind1, ind2, &ind1_new, &ind2_new);
    			ind1_new.evaluate = 1;
    			ind2_new.evaluate = 1;
    		}
    		else	// prechod jedincu bez krizeni
    		{
    			ind1_new = *ind1;
    			ind2_new = *ind2;
    		}
    		// mutace
    		if (mutator(&ind1_new, PMUT)) ind1_new.evaluate = 1;
    		if (mutator(&ind2_new, PMUT)) ind2_new.evaluate = 1;
    		// umisteni potomku do nove populace
    		next_population[i] = ind1_new;
    		next_population[i + 1] = ind2_new;
    	}
    } while (!stop());
}

// *****************************************************************************
// --------------- geneticke operatory a podpurne funkce pro GA ----------------

// generuje cele cislo v rozsahu low-high vcetne
UINT urandom(UINT low, UINT high)
{
    return rand() % (high - low + 1) + low;
}

// vypis chromozomu - ZMENTE SI DLE POTREBY
void gprint(GA_chromosome *genome)
{
//    for (UINT i = 0; i < glength; i++)
//        printf("%d ", genome->chromosome[i]);
//    putchar('\n');

    printf("%u\n", CA_STATES);
    for (UINT j = 0; j < RULES*RULEN; j++){
        if (j % RULEN  == RULEN -1)
            printf("-> %d\n", genome->chromosome[j]);
        else
            printf("%d ", genome->chromosome[j]);
    }
}

// inicializace populace nahodne - ZMENTE SI DLE POTREBY
void initialize(GA_chromosome *genome)
{
    //CMR rules COND:VAL

    for(UINT i = 0; i < RULES; i++){
        for (UINT j = 0; j < RULEN-1; j++){
            
            if( j & 1)
                genome->chromosome[i*RULEN + j] = urandom(0, STATES-1);
            else
                genome->chromosome[i*RULEN + j] = urandom(0, 4);
        }

        //target value        
        genome->chromosome[i*RULEN + (RULEN-1)] = urandom(0, STATES-1);
    }
//    gprint(genome);
//    printf("--------------------------------\n");

}

// krizeni - ZMENTE SI DLE POTREBY
void crossover(GA_chromosome *parent1, GA_chromosome *parent2,
                GA_chromosome *child1, GA_chromosome *child2)
{
    // zde standardni jednobodove krizeni

    UINT cpoint;

    for (UINT i = 0; i < glength; i++)
    {
        if (i % RULEN == 0)
            cpoint = urandom(0, RULEN - 1);

        if ( i % RULEN < cpoint)
        {
            child1->chromosome[i] = parent1->chromosome[i];
            child2->chromosome[i] = parent2->chromosome[i];
        }
        else
        {
            child1->chromosome[i] = parent2->chromosome[i];
            child2->chromosome[i] = parent1->chromosome[i];
        }
    }
}

// mutace - ZMENTE SI DLE POTREBY. je vsak treba zachovat navratovou hodnotu!
BOOL mutator(GA_chromosome *genome, UINT _pmut)
{
    if (urandom(0, unit) <= _pmut)  // mutace s pravdepodobnosti _pmut
    {
        for (UINT i = 0; i < MUTAGENES; i++)
        {
            UINT g = urandom(0, glength - 1);
            int mod = g % RULEN; //pozice v ramci pravidla

            if(mod & 1) //licha pozice - vzdy stav
                genome->chromosome[g] = (STATES-1) - genome->chromosome[g]; // inverze stavu
            else{
                 if(mod == RULEN-1) 
                    genome->chromosome[g] = (STATES-1) - genome->chromosome[g]; // inverze stavu
                 else
                    genome->chromosome[g] = 4 - genome->chromosome[g]; // inverze operace
            }
        
        return 1;   // probehla-li mutace, vratim true...
        }
    }

    return 0;   // ...jinak vracim false
}

// test na zastaveni evoluce - V PRIPADE POTREBY ZMENTE

BOOL stop()
{
    if (generations > 0 && generation == generations)
        return 1;
    if(best.fitness > MFIT)
        return 1;
   
    return 0; 
}

// evaluace fitness pro zadaneho jedince - ZMENTE PRO RESENI SVEHO PROBLEMU

UINT fitness(GA_chromosome *genome, CA & ca)
{
    ca.init();

    vector<UINT *> rules;

    for (int i = 0; i < RULES;i++){
        rules.push_back(&(genome->chromosome[i*RULEN]));
    }
//    gprint(genome);

    ca.setRules(rules);
  //  ca.printRulesCMR();
    double fitness = ca.run(0);
  //  ca.print();
//    cout << "gen: " << generation << " fit: " << fitness << endl;
    cout << "gen: " << generation <<  " UINT fit " <<  (UINT) (fitness * SCALE) << endl;
    //exit(0);
    return (UINT) (fitness * SCALE);
}

// *****************************************************************************
// ------------------------------ hlavni program -------------------------------


int main(int argc, char *argv[])
{
    if (argc < 2){
        cerr << "Not enough arguments, pass me cas file" << endl;
        return 1;
    }

    CA ca;
    ca.parseInit(argv[1]);
    ca.init();
    
    ca.print();

    srand(time(0)); // random seed - NUTNE
    evolve(ca);     //ruleCnt 10 -> X
    
    cout << "EA stopeed at " << generation << " generations";
    cout << "Best chrom" << endl;
    gprint(&best);


    ca.printRulesCMR();

        
    return 0;
}
