#!/sw/bin/python2.7

import getopt
import math
import numpy as np
import os
from scipy.stats import chi2
from sympy.utilities.iterables import multiset_permutations
import sys

#calculating elasticity tensor from fabric tensor
#assumed to be fabric tensor of second kind, F0, F2 or F4

def read_file(filename):
    base = os.path.splitext(filename)[0]
    format = base[-2:]
    print format
    assert(format in ['F2','F4'])
    data = np.load(filename)
    return data, format

def write_elasticity(filename, elasticity, order):
    if order == 1:
        key = 'i'
    elif order == 2:
        key = 'ij'
    elif order == 3:
        key = 'ijk'
    elif order == 4:
        key = 'ijkl'
    elif order == 6:
        key = 'ijklmn'
    C_ijkl = '1111 1112 1113 1121 1122 1123 1131 1132 1133   1112 1212 1213 1221 1222 1223 1231 1232 1233    1113 1213 1313 1321 1322 1323 1331 1332 1333     1121 1221 1321 2121 2122 2123 2131 2132 2133     1122 1222 1322 2122 2222 2223 2231 2232 2233     1123 1223 1323 2123 2223 2323 2331 2332 2333     1131 1231 1331 2131 2231 2331 3131 3132 3133     1132 1232 1332 2132 2232 2332 3132 3232 3233     1133 1233 1333 2133 2233 2333 3133 3233 3333'
    dim = 3
    entries = []

    #generate entry lists
    if (order == 2):
        for i in range(0, dim):
            for j in range(0, dim):
                entry = '{}{}'.format(i,j) 
                entries.append(entry)

    elif (order == 4):
        for i in range(0, dim):
            for j in range(0, dim):
                for k in range(0, dim):
                    for l in range(0, dim):
                        entry = '{}{}{}{}'.format(i,j,k,l) 
                        entries.append(entry)
    else:
        assert(1 == 0)
    
    assert(len(entries) == math.pow(dim, order)) 

    #write data to file
    eFile = filename+'.elas'

    data = []
    for entry in entries:
        if (order == 2):
            data.append( elasticity[entry[0], entry[1]] )
        elif (order == 4):
            print "entry = ",entry
            print "len(elast) = ",len(elasticity)
            print "shape = ",elasticity.shape
            data.append( elasticity[entry] )
        else:
            assert(1 == 0)
    print data

    return    

def calc_rahmoun(data, order):
    dim = 3
    shape = (dim,3)*order

    elasticity = np.zeros( shape )
    return elasticity

def calc_shertzer(data, order):
    dim = 3
    shape = (dim,3)*order

    elasticity = np.zeros( shape )
    return elasticity

def calc_constant():
    dim = 3
    shape = (dim,3)*order

    elasticity = np.zeros( shape )
    return elasticity

def calc_elasticity(files, style): 
    for filename in files:
        print filename
        data,format = read_file(filename)

        order = len(data.shape)

        assert(order in [0,2,4])

        if (order == 0):
            elasticity = calc_constant()

        elif (style == 'shertzer'):
            assert(order == 2)
            elasticity = calc_shertzer(data, order)       

        elif (style == 'rahmoun'):
            elasticity = calc_rahmoun(data, order)       

        write_elasticity(filename, elasticity, order)
  
    return elasticity 


if __name__ == "__main__":
    optlist,args = getopt.getopt(sys.argv[1:],'',longopts=['style='])
    style = None
    for item in optlist:
        if (item[0].lower() == '--style'):
            style = item[1].lower()
    calc_elasticity(args, style) 
