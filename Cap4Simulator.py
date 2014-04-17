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
    r = 0.1                 # Pipe radius
    dx = 0.2                # Discretization lenght delta_x
    ncells = 3              # Number of cells in domai

    # Fluid properties
    rho = 1                 # Fluid density
    mu = 0.01               # Fluid viscosity
    msrc = 0                # Mass Source term per volume unity

    # Initial properties
    v_ini = 0               # Initial Condition for v
    p_ini = 0               # Initial Condition for p

    # Boundary Condition
    # Left Boundary Condtion
    lbc_t = 1                # LBC Type (0 - Pressure / 1 - Velocity)
    lbc = 10                 # LBC Value

    # Right Boundary Condtion
    rbc_t = 0               # LBC Type (0 - Pressure / 1 - Velocity)
    rbc = 0                 # LBC Value

    # Calculated Values
    Ar = np.pi*r*r           # Pipe Area
    Re = rho*v_ini*2*r/mu       # Reynolds Number

    # Creating Grid
    grid = GridFluid(lbc_t, lbc, rbc_t, rbc, msrc)

    # Creating Model
    model = Model_COUPLE_CDS()

    for i in range(ncells):
        grid.add_cell(Ar, dx, rho, mu, v_ini, p_ini)

    print "Initial Reynolds Number -->", Re
    print "Pipe Area  -->", Ar   
    print

    for j in range(50):
        print "Iteration ", (j+1)
        #Construindo Matrix Inicial do Sistema
        A = np.matrix(model.build_Jacobian_matrix(grid))
        b = - np.matrix(model.build_Residual_Vector(grid))
        x = (A.I*b.T).A1
        
#        print "JACOBIAN MATRIX"    
        print A
        print
#        
#        print "RESIDUAL VECTOR"
        print b
        print
#        
#        print "SOLUTION"
        print x
        print
    
        #Atribuindo os resultados
        for i in range(ncells):
            newp = (grid.p(i) + x[2*i])
            print "NEW P --> ", newp
            grid.set_p(i, newp)
            
            newv = (grid.v_rh(i) + x[2*i + 1])
            print "NEW v --> ", newv
            grid.set_v_rh(i, newv)
        print
#    grid.print_phi()
#
#    print
#    print A
#    print
#    print b
##    print grid.get_all_x()
##    print grid.get_all_T()
#    grid.plot_T()


if __name__ == '__main__':

    run()
