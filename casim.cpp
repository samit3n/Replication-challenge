
/*
 *
 *  2D Cellular Automata simulator for BIN 2015 project purpose
 *  
 *  Author: Vojtech Dvoracek
 *  email: xdvora0y@stud.fit.vutbr.cz
 *  
 */


#include "casim.h"


using namespace std;

CA::CA()
{
    parsedRules = false;

}

void CA::print(void)
{

        for(unsigned j = 0; j < CA_HEIGHT;j++){
            for(unsigned i = 0; i < CA_WIDTH;i++){

            if(lattice[i][j] == 0)
                cout << BLUE <<  lattice[i][j] << RESET << " ";
            if(lattice[i][j] == 1)
                cout << RED <<  lattice[i][j] << RESET << " ";
            if(lattice[i][j] == 2)
                cout << CYAN <<  lattice[i][j] << RESET << " ";
            if(lattice[i][j] == 3)
                cout << GREEN <<  lattice[i][j] << RESET << " ";
            if(lattice[i][j] == 4)
                cout << YELLOW <<  lattice[i][j] << RESET << " ";
        }
        cout << "\n";
    }
}

CA::~CA()
{
    if(parsedRules == true){
        //freeing rules
        for(vector<unsigned *>::iterator it = rules.begin(); it != rules.end();it++){
            delete[] *it;
        }
    }
        
}

CA_STATE::CA_STATE(int cx, int cy, unsigned cstate)
{
    x = cx;
    y = cy;
    state = cstate;   
}

void CA_STATE::print(void)
{
    cout << "x: " << x << " y: " << y << " st: " << state << endl;
}

int CA::stable (void)
{
   int flag = 1;

    for(unsigned i = 0; i < CA_WIDTH; i++){
        for(unsigned j = 0; j < CA_HEIGHT; j++){

            if(lattice[i][j] != lastLat[i][j])
                flag = 0;
        }
    }

    return flag;
}

bool CA::isNeigh(int x1, int y1, int x2,int y2)
{
    if(abs(x1-x2) == 1 && abs(y1-y2) == 0)
        return true;
    else if(abs(x1-x2) == 0  && abs(y1-y2) == 1)
        return true;
    else
        return false;
}



void CA::init(void)
{

    //init lattice to zero
    for(unsigned i = 0; i < CA_WIDTH;i++){
        for(unsigned j = 0; j < CA_HEIGHT;j++){
            lattice[i][j] = 0;
            tmpLat[i][j] = 0;
        }   
    }

    //odd / even check
    midX = (CA_WIDTH & 1)? (CA_WIDTH+1) /2 : CA_WIDTH /2;
    midY = (CA_HEIGHT & 1)? (CA_HEIGHT+1) /2 : CA_HEIGHT /2;
    //centering mid
    midX -= (INIT_DIMX & 1)? (INIT_DIMX+1) /2 : INIT_DIMX /2;
    midY -= (INIT_DIMY & 1)? (INIT_DIMY+1) /2 : INIT_DIMY /2;

    //CA initial state 
    for (vector<CA_STATE>::iterator it = initVect.begin(); it != initVect.end(); it++){

        unsigned posX = midX +  (*it).x;
        unsigned posY = midY + (*it).y;
        
      //  lattice[posX][posY] = (*it).state;
        lattice[posX][posY] = (*it).state;
    }
    
    for(unsigned i = 0; i < MES_RNDS;i++)
    {
        FgMes[i].resize( CA_STATES);
        fill(FgMes[i].begin(),FgMes[i].end(),0);

        FpMes[i].resize( neigh.size() );
        fill(FpMes[i].begin(),FpMes[i].end(),0);

    }

    //create vector of component neighbors

    Replics.clear();

}

void CA::getCoords(int * x, int *y)
{

    // lattice overflow
    
    if (*x < 0) 
        *x = CA_WIDTH + *x;
    else if (*x >= CA_WIDTH)
        *x = *x - CA_WIDTH;

    if (*y < 0)
        *y = CA_HEIGHT + *y;
    else if (*y >= CA_HEIGHT)
        *y = *y - CA_HEIGHT;

}
unsigned CA::getVal(int x, int y)
{   

    getCoords(&x,&y);
    return lattice[x][y];
        
}
void CA::setVal(int x, int y, unsigned val)
{

    lattice[x][y] = val;

}
void CA::setValTmp(int x, int y, unsigned val)
{

    tmpLat[x][y] = val;

}
int CA::setRules(vector<unsigned *> &rls)
{
    rules.clear();
    rules = rls;   
}


//return appropriate coords in neighborhood
//for purposes of next state computing
//cellNum coresponds to neighborhood
//counted from top-left to bottom-right

void CA::getShift(int cellNum, int * x, int * y)
{

    switch(cellNum){
        case 0: *y -= 1;break;
        case 2: *x -= 1;break;
//        case 4: no change
        case 6: *x += 1;break;
        case 8: *y += 1;break;      
    }
}

void CA::cpLattice(bool last)
{
   
    if(last){

        for(unsigned i = 0; i < CA_WIDTH;i++){
            for(unsigned j = 0; j < CA_HEIGHT;j++){
                lastLat[i][j] = lattice[i][j];
            }   
        }

    }else{
        for(unsigned i = 0; i < CA_WIDTH;i++){
            for(unsigned j = 0; j < CA_HEIGHT;j++){
                lattice[i][j] = tmpLat[i][j];
            }   
        }
    }

}

//CA fitness
unsigned CA::fitness(unsigned round)
{

    unsigned score = 0; //score of intiail pattern
    unsigned copies = 0;

    unsigned fitness = 0, act;
    int posX, posY;
    bool match;

    for(unsigned i = 0; i < CA_WIDTH - INIT_DIMX; i++){ //x axes
        for(unsigned j = 0; j < CA_HEIGHT - INIT_DIMY;j++){ //y axes

            score = 0;
            act = lattice[i][j];
            //counting states for tendency fitness (Fg)

            FgMes[round][act] += 1;
        
            if( act != 0){
                //searching neighbors mv(t)
                for(auto Nit = neigh.begin(); Nit != neigh.end(); Nit++){
                    if ((*Nit)[0].state != act)
                        continue;
                    
                    match = true;
                    for(auto Vit = (*Nit).begin()+1; Vit != (*Nit).end(); Vit++){
                        //iterate from 2nd item, first is source state

                        posX = i + (*Vit).x;
                        posY = j + (*Vit).y;

                        getCoords(&posX, &posY);

                        if(lattice[posX][posY] != (*Vit).state)
                            match = false;

                    }
                    if(match == true)
                        FpMes[round][Nit - (neigh.begin())] +=  neighCnt[(*Nit)[0].state];
                        
                }
            }

            //finding full pattern
            for (unsigned it = 0; it < initVect.size();it++){
//               initVect[it].print();

                posX = i +  (initVect[it]).x;
                posY = j + (initVect[it]).y;


                if( lattice[posX][posY] == (initVect[it]).state)
                    score++;
                    
            }

            if(score == SQ){
                copies++;
            }
        }
    
    }
    Replics.push_back(copies);

    return 0; 
    
}

//counts potential neighborhood for given state
double CA::cntMaxNeigh(unsigned state, unsigned Mvt)
{
    vector<unsigned> cnts;

    for(int i = 0; i < neigh.size();i++)
    {
        if( (neigh[i][0]).state == state){
            cnts.push_back(neighCnt[i]);
        }
    }
    double frac; 

    frac = Mvt / cnts.size();

    double res = 0;
    
    for(auto it = cnts.begin(); it != cnts.end();it++)
        res += frac *(*it);

    return res;
}


double CA::run(unsigned stop){

    bool stopFlag = stop > 0 ? true : false; 
    unsigned round = 0;
    bool change;
    unsigned latX, score;
    unsigned op;

    double Fg = 0, Fp = 0, Fr = 0, Fs = 0, fit;
    double tmpFg = 0, tmpSig = 0;

    int flag;
    int tmpX, tmpY;
    

    
    while(1){
        
        fitness(round);

        if(round == MES_RNDS-1){ //fitness computation and return
         
    
            for(unsigned j = 1; j < FgMes[1].size();j++){

                tmpFg = 0; //sum fo Mvt at time <1,t>
                tmpSig = 0;

                for(unsigned t = 1; t < MES_RNDS;t++){

                    if (FgMes[t][j] > FgMes[t-1][j])
                        tmpFg += 1;
                    else if (FgMes[t][j] ==  FgMes[t-1][j])
                        tmpFg += 0.5;

                //sigma(t)
                    if(FgMes[t][j] <= 1)
                        tmpSig  += 0;
                    else{
                        unsigned mvtSum = 0;

                        for (int i =0 ;i < neigh.size(); i++){
                            if( neigh[i][0].state == j )
                                mvtSum += FpMes[t][i];
                        }
                        
                        tmpSig += mvtSum / (FgMes[t][j] * neighCnt[j]);
                    }

                }
                
                Fp += tmpSig * neighCnt[j];
                Fg += tmpFg; //sum rounds
            }

            //normalize to <0,1>
            Fg = ( 1/( (double) MES_RNDS*(CA_STATES-1)) ) * Fg;
            Fp = ( 1/( (double) MES_RNDS*neighCntSum) ) * Fp;
            
            vector<unsigned>::iterator repMax = max_element(Replics.begin(), Replics.end());
          
            cout << "replics " << *repMax << endl;
    
            Fr =  1 / (1 + exp( -((*repMax)-THETA) ));
            
            cout << "Fg: " << Fg << " Fp: "  << Fp << " Fr: " << Fr << endl;


            fit = Fg * w1 + Fp * w2 + Fr * w3;

            return fit;
        }

        cpLattice(true);
        

        //lattice iteration
        for(unsigned i = 0; i < CA_WIDTH;i++){ 
            for(unsigned j = 0; j < CA_HEIGHT;j++){
                
                change = false; 
                //iterate over rules
                for(vector<unsigned int *>::iterator it = rules.begin(); it != rules.end();  it++){

                    flag = 1;
                    for(unsigned x = 0; x < COND_COUNT*2;x++){             //iterate over conditions
                       
                        tmpX = (int) i;   
                        tmpY = (int) j;

                        if( x % 2 != 0){
                                
                            op = (*it)[x-1];
                            getShift(x-1, &tmpX, &tmpY);
                    

                            unsigned val = getVal(tmpX,tmpY);
                            unsigned condVal = (*it)[x];

                            switch(op){
                                    
                                case EQ:
                                    if(val != condVal)
                                        flag = 0;
                                    break;
                                case NEQ:
                                    if(val == condVal)
                                        flag = 0;
                                    break;
                                case MOREQ:
                                    if(val < condVal)
                                        flag = 0;
                                    break;
                                case LESEQ:
                                    if (val > condVal)
                                        flag = 0;
                                    break;
                            }
                        }
                    }

                    if (flag){
                        setValTmp(i,j,(*it)[2*COND_COUNT] ); //last elem is targe value
                        change = true;
                    }
                }
                if(!change)
                    setValTmp(i,j, lattice[i][j]);
            }   
        }
        cpLattice(false);
        round++;
   }
}


/*
 * initial pattern parsing
 */

int CA::parseInit(string filename)
{
    ifstream in;

    in.open(filename.c_str(), ios::in);

    if (!in){
        cerr << "Cannot open file " << filename << endl;
        return 1;
    }


    string str ( (istreambuf_iterator<char>(in)),istreambuf_iterator<char>() );
    string token;

    unsigned spcs = count(str.begin(), str.end(), ' ');
    unsigned nls = count(str.begin(), str.end(), '\n');
    
    INIT_DIMY = nls;
    INIT_DIMX = (spcs / nls) +1;

    unsigned row = 0;
    unsigned col = 0;
    
    //tokenizing
    for (string::iterator it = str.begin(); it != str.end();it++){

        if(*it == ' '){
            initVect.push_back( CA_STATE(col, row, stoi(token)) );
            token.erase();
            col++;
        }else if(*it == '\n'){
  //          INIT_DIMX = col+1;
            initVect.push_back( CA_STATE(col ,row,  stoi(token)) );
            token.erase();
            col = 0;
            row++;
        }
        else{
            token.push_back(*it);
        }
    }
    SQ = INIT_DIMX * INIT_DIMY;

    in.close();
    
    vector<CA_STATE> tmp;
   
    for (int i = 0; i < CA_STATES; i++)
        initComps[i] = 0;

    for(int i = 0; i < initVect.size();i++){
        //distribution of states in initial pattern
        initComps[initVect[i].state]++;

        if(initVect[i].state != 0){
            tmp.push_back(initVect[i]);

            for(int j = 0; j < initVect.size();j++){

                if(initVect[j].state != 0 && isNeigh(initVect[i].x, initVect[i].y, initVect[j].x, initVect[j].y)) 
                {

                    CA_STATE tmpSt(initVect[j].x - initVect[i].x, initVect[j].y - initVect[i].y, initVect[j].state);

                    //tmp.print();
                    tmp.push_back(tmpSt);
                            
                }
            }
            neigh.push_back(tmp);
          //  neighCnt.push_back(tmp.size()-1);
            tmp.clear();
        }
    }


    neighCnt.resize(CA_STATES);
    fill(neighCnt.begin(), neighCnt.end(), 0);


    for(auto it = neigh.begin(); it != neigh.end(); it++){
        if (neighCnt[(*it)[0].state] < ((*it).size()-1 ))
            neighCnt[(*it)[0].state] = (*it).size()-1; 
    }
    neighCntSum = accumulate(neighCnt.begin(), neighCnt.end(), 0);

    return 0;
}

int CA::parseRules(string filename)
{

    parsedRules = true;
    ifstream in;

    in.open(filename.c_str(), ios::in);

    if (!in){
        cerr << "Cannot open file " << filename << endl;
        return 1;
    }


    string str ( (istreambuf_iterator<char>(in)),istreambuf_iterator<char>() );

    //simple rule check
    unsigned spcs = count(str.begin(), str.end(), ' ');
    unsigned nls = count(str.begin(), str.end(), '\n')+1;

    if ( (spcs % nls) != 0){
        cerr << "spaces: " << spcs << " nls: " << nls << endl;
        cerr << "Rule len mismatch in file " << filename << endl;
        return 1;
    }
        

    string token;

    vector<unsigned> tmp;
    
    //tokenizing
    for (string::iterator it = str.begin(); it != str.end();it++){

        if(*it == ' '){

            tmp.push_back(stoi(token));
            token.erase();

        }else if(*it == '\n'){

            tmp.push_back( stoi(token));
            token.erase();
            unsigned * arr = new unsigned[tmp.size()];
            //ugly but sigsegv proof
            copy(tmp.begin(), tmp.end(), arr);      
            tmp.clear();
            rules.push_back(arr);
        }
        else{
            token.push_back(*it);
        }
    }

    in.close();
    return 0;

}

void CA::printInit(void)
{

    for(vector<CA_STATE>::iterator it = initVect.begin(); it != initVect.end();it++)
        it->print();


}

void CA::printRulesCMR()
{
    cout << STATES << endl;

    for(vector<unsigned *>::iterator it = rules.begin(); it != rules.end(); it++)
    {
        for (unsigned i = 0; i < 2*COND_COUNT; i++){
            cout << (*it)[i] << " ";
        }
        cout << "-> " << (*it)[2*COND_COUNT] << "\n";
    }

}


