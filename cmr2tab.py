#!/usr/bin/env python2

import sys
from itertools import izip, product

# rules {DC = 0, EQ, NEQ, MOREQ, LESEQ}

def main(inf,outf ):

    input = open(inf, "r")
    output = open(outf, 'w')

    
    states = int(input.readline()[:-1])
    tabRules = list()

    for line in input:
	if line[-1] == '\n':
	  line = line[:-1]
	  
	sp = line.split(' ')
	
        sp.remove('->')
        sp =  map(int, sp)
        rules = list()
        pairs = izip( sp[::2], sp[1::2])

        for tup in pairs:

            if tup[0] == 0:
                rules.append( range(0,states) )
            if tup[0] == 1:
                rules.append( [tup[1]] )
            if tup[0] == 2:
                x = range(0,states)
                x.remove(tup[1])
                rules.append( x )
            if tup[0] == 3:
                rules.append( range(tup[1], states) )
            if tup[0] == 4:
                rules.append( range(0, tup[1]+1) )

        rules.append([sp[-1]])
        resRules = product(*rules)

        for r in resRules:
           tabRules.append(r)
            

    output.write(str(states) + "\n")
    for rule in tabRules:
        o = map(str, rule)
        output.write(" ".join(o) + "\n")
        
    input.close()
    output.close()



if __name__ == '__main__':
    
    main(sys.argv[1], sys.argv[2])
