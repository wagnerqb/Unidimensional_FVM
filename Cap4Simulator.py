# -*- coding: iso-8859-1 -*-
"""
@author: Wagner Queiroz Barros
@date:   Fri Apr 11 13:24:24 2014
@email:  wagnerqb@gmail.com
@brief:  Simulator for Navier-Stokes problem, using COUPLE method
"""


from __future__ import division
import numpy as np
from GridFluid import *
from Model_COUPLE_CDS import *


def run():

    # Pipe properties
    dx = 0.1
    ncells = 20
    a_ = 1       # Area de saida
    b_ = 2       # Area de entrada

    # Fluid properties
    rho = 100                 # Fluid density
    mu = 100                 # Fluid viscosity
    msrc = 0.               # Mass Source term per volume unity

    # Initial properties
    v_ini = 1               # Initial Condition for v
    p_ini = 0               # Initial Condition for p

    # Boundary Condition
    # Left Boundary Condtion
    lbc_t = 1                # LBC Type (0 - Pressure / 1 - Velocity)
    lbc = 1                 # LBC Value

    # Right Boundary Condtion
    rbc_t = 0               # LBC Type (0 - Pressure / 1 - Velocity)
    rbc = 0                # LBC Value

    # Creating Grid
    grid = GridFluid(lbc_t, lbc, rbc_t, rbc, msrc)

    # Creating Model
    model = Model_COUPLE_CDS()

    v_real = []
    p_real = []
    for i in range(ncells):
        A = b_- (b_-a_)*(.5+i)/ncells
        A_ = b_- (b_-a_)*(1+i)/ncells
        v_real.append(b_*lbc/A_)
        p_real.append(rbc+0.5*rho*(lbc*b_)**2*(1/a_**2-1/A**2))
        grid.add_cell(A, dx, rho, mu, v_ini, p_ini)
 
    for j in range(3):
        print '#'*40
        print "Iteration ", (j+1)
        #Construindo Matrix Inicial do Sistema
        A = np.matrix(model.build_Jacobian_matrix(grid))
        b = - np.matrix(model.build_Residual_Vector(grid))
        x = (A.I*b.T).A1
        
##        print "JACOBIAN MATRIX"    
#        print A
#        print
##        
##        print "RESIDUAL VECTOR"
#        print b
#        print
##        
##        print "SOLUTION"
#        print x
#        print
    
        #Atribuindo os resultados
        new_p = []
        new_v = []
        for i in range(ncells):
            newp = (grid.p(i) + x[2*i])
            new_p.append(newp)
            grid.set_p(i, newp)
            
            newv = (grid.v_rh(i) + x[2*i + 1])
            new_v.append(newv)
            grid.set_v_rh(i, newv)
            
        print 'New p:\n', new_p
        print 'New v:\n', new_v
        print

    import matplotlib.pyplot as plt
    
    print '*'*30
    print 'L1_v: ', max(np.array(new_v)-np.array(v_real))
    print 'L1_p: ', max(np.array(new_p)-np.array(p_real))
    #Velocidade
    plt.subplot(211)
    plt.plot(new_v, 'ko', label='Numerico')
    plt.plot(v_real, 'b', label='Real')
    plt.title('Velocidade')
    plt.grid()
    
    #Pressao
    plt.subplot(212)
    plt.plot(new_p, 'ko', label='Numerico')
    plt.plot(p_real, 'b', label='Real')
    plt.title(u'Pressão')
    plt.grid()
    plt.show()
    

#    for i in range(len(new_v)):
#        print grid.v(i)*grid.A(i)
if __name__ == '__main__':

    run()
