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

def kronecker(delta_dimension, data_dimension):

    if (delta_dimension == 1):
        delta = 1

    elif (delta_dimension == 2):
        delta = np.zeros( (data_dimension, data_dimension) )
        for i in range(0, data_dimension):
            delta[i,i] = 1
            
    elif (delta_dimension == 3):
        delta = np.array( (data_dimension, data_dimension, data_dimension) )
        for i in range(0, data_dimension):
            delta[i,i,i] = 1

    elif (delta_dimension == 4):
        delta = np.ndarray( (data_dimension, data_dimension, data_dimension, data_dimension) )
        delta[:,:,:,:] = 0
        for i in range(0, data_dimension):
            delta[i,i,i,i] = 1

    return delta

def read_file(filename):
    base = os.path.splitext(filename)[0]
    format = base[-2:]
    assert(format in ['F2','F4', 'N2', 'N4'])
    data = np.load(filename)
    return data, format

def write_elasticity(filename, elasticity):

    elas_order = len( elasticity.shape )

    if elas_order == 1:
        key = 'i'
    elif elas_order == 2:
        key = 'ij'
    elif elas_order == 3:
        key = 'ijk'
    elif elas_order == 4:
        key = 'ijkl'
    elif elas_order == 6:
        key = 'ijklmn'
#    C_ijkl = '1111 1112 1113 1121 1122 1123 1131 1132 1133   1112 1212 1213 1221 1222 1223 1231 1232 1233    1113 1213 1313 1321 1322 1323 1331 1332 1333     1121 1221 1321 2121 2122 2123 2131 2132 2133     1122 1222 1322 2122 2222 2223 2231 2232 2233     1123 1223 1323 2123 2223 2323 2331 2332 2333     1131 1231 1331 2131 2231 2331 3131 3132 3133     1132 1232 1332 2132 2232 2332 3132 3232 3233     1133 1233 1333 2133 2233 2333 3133 3233 3333'
    dim = 3
    entries = []

    #generate entry lists
    if (elas_order == 2):
        for i in range(0, dim):
            for j in range(0, dim):
                entry = '{}{}'.format(i,j) 
                entries.append(entry)

    elif (elas_order == 4):
        for i in range(0, dim):
            for j in range(0, dim):
                for k in range(0, dim):
                    for l in range(0, dim):
                        entry = '{}{}{}{}'.format(i,j,k,l) 
                        entries.append(entry)
    else:
        assert(1 == 0)
    
    assert(len(entries) == math.pow(dim, elas_order)) 

    #write data to file
    base = os.path.splitext(filename)[0]
    eFile = base+'.elas'

    data = []
    for entry in entries:
        entryList = []
        for letter in entry:
            entryList.append(int(letter))
        entryList = tuple(entryList)
        data.append( elasticity[entryList] )

    #convert data to string
    dataStr = ' '.join([str(x) for x in data])
    print dataStr

    #write string to file
    with open(eFile,'w') as f:
        f.write(dataStr)

    return    

def calc_rahmoun(data, order):
    dim = 3
    shape = (dim,)*order

    elasticity = np.zeros( shape )

    Ncp = 1.0 
    R = 1.0 
    V = 1.0
    rho = 1.0
    E = 1.0
    G = 1.0
    kn = math.pi * rho*rho * E / (2.0 * R) 
    ks = math.pi * rho*rho * G / (2.0 * R)

    if (order == 2):
        d2 = kronecker(2)
        dij_dkl = np.einsum('ij,kl->ijkl', d2, d2)
        dik_djl = np.einsum('ik,jl->ijkl', d2, d2)
        Fij_dkl = np.einsum('ij,kl->ijkl', data, d2)
        Fik_djl = np.einsum('ik,jl->ijkl', data, d2)
        dij_Fkl = np.einsum('ij,kl->ijkl', d2, data)
        dik_Fjl = np.einsum('ik,jl->ijkl', d2, data)
        elasticity = 4.0*Ncp*R*R/(7.0 * V) * (kn * (-1./5.*(dij_dkl + 2.0 * dik_djl) + (Fij_dkl + dij_Fkl) + 2.0 * (Fik_djl + dik_Fjl) ) +
            - ks * ( -1./5. * (dij_dkl + 2.*dik_djl) + (Fij_dkl + dij_Fkl) - 1.5 * (Fik_djl + dik_Fjl) ) )
    elif (order == 4):
        #needs actual definition
        elasticity = data 
        

    return elasticity

def calc_cowin(data, order):
    dim = 3
    shape = (dim,)*order

    elasticity = np.zeros( shape )

    Ncp = 1.0 
    R = 1.0 
    V = 1.0
    rho = 1.0
    E = 1.0
    G = 1.0
    kn = math.pi * rho*rho * E / (2.0 * R) 
    ks = math.pi * rho*rho * G / (2.0 * R)

    d2 = kronecker(2)
    dij_dkl = np.einsum('ij,kl->ijkl',d2,d2)
    Fij_dkl = np.einsum('ij,kl->ijkl',data,d2)
    dij_Fkl = np.einsum('ij,kl->ijkl',d2,data)
    dij_Fkq_Fql = np.einsum('ij,kq,ql->ijkl',d2,data,data)
    dkl_Fiq_Fqj = np.einsum('kl,iq,qj->ijkl',d2,data,data)
    Fij_Fkl = np.einsum('ij,kl->ijkl',data,data)
    Fij_Fkq_Fql = np.einsum('ij,kq,ql->ijkl',data,data,data)
    Fis_Fsj_Fkl = np.einsum('is,sj,kl->ijkl',data,data,data)
    Fis_Fsj_Fkq_Fql = np.einsum('is,sj,kq,ql->ijkl',data,data,data,data)
    dki_dlj = np.einsum('ki,lj->ijkl',d2,d2)
    dli_dkj = np.einsum('li,kj->ijkl',d2,d2)
    Fik_dlj = np.einsum('ik,lj->ijkl',data,d2)
    Fkj_dli = np.einsum('kj,li->ijkl',data,d2)
    Fil_dkj = np.einsum('il,kj->ijkl',data,d2)
    Flj_dki = np.einsum('lj,ki->ijkl',data,d2)
    Fir_Frk_dlj = np.einsum('ir,rk,lj->ijkl',data,data,d2)
    Fkr_Frj_dli = np.einsum('kr,rj,li->ijkl',data,data,d2)
    Fir_Frl_dkj = np.einsum('ir,rl,kj->ijkl',data,data,d2)
    Flr_Frj_dik = np.einsum('lr,rj,ik->ijkl',data,data,d2)

    dij_Fkl = np.einsum('ij,kl->ijkl',d2,data)
    elasticity = a1 * dij_dkl + a2 * (Fij_dkl + dij_Fkl) + a3 * (dij_Fkq_Fql + dkl_Fiq_Fqj) + \
        b1 * Fij_Fkl + b2 * (Fij_Fkq_Fql + Fis_Fsj_Fkl) + b3 * (Fis_Fsj_Fkq_Fql) + \
        c1 * (dki_dlj + dli_dkj) + c2 * (Fik_dlj + Fkj_dli + Fil_dkj + Flj_dki) + \
        c3 * (Fir_Frk_dlj + Fkr_Frk_dli + Fir_Frl_dkj + Flr_Frj_dik)
    return elasticity

def calc_shertzer(data, order, scaling):
    dim = 3
    shape = (dim,)*order
    assert( len(data.shape) == 2)
    assert( order == 2)

    elasticity = np.zeros( shape )

    #Shertzer Eq. 5.16
    rho = 1.0 #mean contact circle radius
    R = 10.0 #mean particle radius
    phi = 0.5; #volume fraction 
    N = 3; #average number of contacts per particle
    E = 20.0; #elastic modulus of solid
    G = 1.0; #shear modulus of solid
    zeta = G/E; 
    lambda_star = 1./20. * math.pow(rho/R, 2) * phi * N * E * (1.0 - zeta) 
    G_star = 1./40. * math.pow(rho/R, 2) * phi * N * E * (2.0 + 3.0 * zeta) 
    elasticity = math.pow(scaling, -2) * (lambda_star * np.einsum('ij,kl->ijkl',data,data) + 2.0 * G_star * np.einsum('ik,jl->ijkl',data,data) )

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

        #determine shertzer's scaling coefficient, f
        if (format[0] == 'N'):
            scalingF = 1./3.
        elif (format[0] == 'F'):
            scalingF = 1.0

        if (order == 0):
            elasticity = calc_constant()

        elif (style == 'shertzer'):
            assert(order == 2)
            elasticity = calc_shertzer(data, order, scalingF)

        elif (style == 'rahmoun'):
            elasticity = calc_rahmoun(data, order)       

        elif (style == 'cowin'):
            elasticity = calc_cowin(data, order)       

        write_elasticity(filename, elasticity)
  
    return elasticity 


if __name__ == "__main__":
    optlist,args = getopt.getopt(sys.argv[1:],'',longopts=['style='])
    style = None
    for item in optlist:
        if (item[0].lower() == '--style'):
            style = item[1].lower()
    calc_elasticity(args, style) 
