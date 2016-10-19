

/*
 *  2D Cellular Automata simulator for BIN 2015 project purpose
 *  
 *  Author: Vojtech Dvoracek
 *  email: xdvora0y@stud.fit.vutbr.cz
 *  
 */

#ifndef CASIM_H
#define CASIM_H

#include <cstdio>
#include <vector>
#include <array>
#include <string>
#include <iostream>
#include <fstream>
#include <numeric>
#include "params.h"
#include <algorithm>
#include <cmath>
#include <stdlib.h>

/**
* Some terminal coloring escape sequences
*/

#define RESET   "\033[0m"
#define BLACK   "\033[30m"      /* Black */
#define RED     "\033[31m"      /* Red */
#define GREEN   "\033[32m"      /* Green */
#define YELLOW  "\033[33m"      /* Yellow */
#define BLUE    "\033[34m"      /* Blue */
#define MAGENTA "\033[35m"      /* Magenta */
#define CYAN    "\033[36m"      /* Cyan */
#define WHITE   "\033[37m"      /* White */
#define BOLDBLACK   "\033[1m\033[30m"      /* Bold Black */
#define BOLDRED     "\033[1m\033[31m"      /* Bold Red */
#define BOLDGREEN   "\033[1m\033[32m"      /* Bold Green */
#define BOLDYELLOW  "\033[1m\033[33m"      /* Bold Yellow */
#define BOLDBLUE    "\033[1m\033[34m"      /* Bold Blue */
#define BOLDMAGENTA "\033[1m\033[35m"      /* Bold Magenta */
#define BOLDCYAN    "\033[1m\033[36m"      /* Bold Cyan */
#define BOLDWHITE   "\033[1m\033[37m"      /* Bold White */

using namespace std;


class CA_STATE{
    
    public:

    CA_STATE(int cx, int cy, unsigned cstate);
    void print();
    int x,y;
    unsigned state;
};

/**
* 2D Cellular automata simulator
*
*/

class CA{


    public:

    CA();
    ~CA();

    void init(void);
    double run(unsigned stop = 0);          //run simulation
    void print(void);

    unsigned fitness(unsigned round);       //returns CA fitness function
    int setRules(vector<unsigned *> &rls);
    int parseInit(string filename);
    void printInit(void);
    int parseRules(string filename);
    void printRulesCMR();
    void getCoords(int * x, int *y);

    private:

    void getShift(int cellNum, int * x, int * y);
    int stable(void);
    void cpLattice(bool last=false);

    unsigned getVal(int x, int y);
    void setVal(int x, int y, unsigned val);
    void setValTmp(int x, int y, unsigned val);
    bool isNeigh(int x1, int y1, int x2,int y2);

    unsigned int lattice[CA_WIDTH][CA_HEIGHT];
    unsigned int tmpLat[CA_WIDTH][CA_HEIGHT];
    unsigned int lastLat[CA_WIDTH][CA_HEIGHT];

    unsigned INIT_DIMX, INIT_DIMY, SQ, midX, midY;

    bool parsedRules;

    vector<unsigned int *> rules;
    vector<CA_STATE> initVect;

    //fitness stuff
    vector<unsigned> FgMes[MES_RNDS];
    vector<unsigned> FpMes[MES_RNDS];
    vector<unsigned> Replics;

    double cntMaxNeigh(unsigned state, unsigned Mvt);
    vector< vector<CA_STATE> > neigh;
    vector<unsigned> neighCnt;
    unsigned neighCntSum;
    unsigned initComps[CA_STATES];

    
};


typedef enum  {DC = 0, EQ, NEQ, MOREQ, LESEQ} RULE_CONDS;

#endif
