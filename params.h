//////////////////////////////////////////////////////////////
// VZOROVA IMPLEMENTACE JEDNODUCHEHO GENETICKEHO ALGORITMU  //
// JAZYK C, PREKLAD S OHLEDEM NA NORMU C99                  //
// (c) MICHAL BIDLO, 2011                                   //
// UVEDENY KOD SMI BYT POUZIT VYHRADNE PRO POTREBY PROJEKTU //
// V KURZU BIOLOGII INSPIROVANE POCITACE NA FIT VUT V BRNE  //
//////////////////////////////////////////////////////////////



#ifndef PARAMS_H
#define PARAMS_H

/*
 * Genetic Algorithm params
 *
 */

// pro vypis prubehu evoluce na stdout, jinak zakomentujte
//#define DEBUG

//----------------------- parametry genetickeho algoritmu ----------------------
// pravdepodobnost mutace
#define PMUT 65
// pocet mutovanych genu v chromozomu
#define MUTAGENES 25
// pravdepodobnost krizeni
#define PCROSS 80
// pocet jedincu v turnajove selekce
#define TOUR 15
// velikost populace
#define POPSIZE 30
// maximalni pocet generaci
#define GENERATIONS 200000
// delka chromozomu

#define RULES 60 //number of evolving rules
#define RULEN 11//len of 1 CMR rules
#define STATES 3 //CA states
#define CA_STATES 3 //CA states again
#define MFIT 800 //minimal fitness to end
#define SCALE 1000 //fitness  = SCALE * <0,1>
// CASIM params
#define CA_WIDTH  30
#define CA_HEIGHT 30
#define COND_COUNT 5 //number of CMR conditions


//CA genereations
const unsigned MES_RNDS = 15;

//fitness sigmoid inflex param
const double THETA  = 4.0;

//weights for fitness proportion
//they must sum to 1
//fitness = Fg * w1 + Fp * w2 + Fr * w3

const double w1 = 0.3;
const double w2 = 0.4;
const double w3 = 0.3;



#endif
