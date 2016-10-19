/*
 *
 *  main function for casim standalone running
 *
 */

#include "casim.h"


int main(int argc, char * argv[])
{
    if(argc < 3){
        cerr << "Not enough arguments\n";
        cerr << "Usage: ./casim <pattern.cas> <rules.cmr>\n";
        return 1;

    }

    CA ca;

    ca.parseInit(argv[1]);
    ca.parseRules(argv[2]);
    ca.init();
    ca.print();
    ca.run();
    ca.print();


}
