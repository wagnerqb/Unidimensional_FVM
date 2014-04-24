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
    r = 150E-3/2                 # Pipe radius
    dx = 2               # Discretization lenght delta_x
    ncells = 5              # Number of cells in domai

    # Fluid properties
    rho = 999                 # Fluid density
    mu = 1E-3*rho               # Fluid viscosity
    msrc = 0.               # Mass Source term per volume unity

    # Initial properties
    v_ini = 5               # Initial Condition for v
    p_ini = 15000               # Initial Condition for p

    # Boundary Condition
    # Left Boundary Condtion
    lbc_t = 1                # LBC Type (0 - Pressure / 1 - Velocity)
    lbc = 5.66                 # LBC Value

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
        grid.add_cell(Ar, dx, rho, mu, v_ini, p_ini-p_ini/ncells*i)

    print "Initial Reynolds Number -->", Re
    print "Pipe Area  -->", Ar
    print

    for j in range(4):
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
        Re = rho*lbc*2*r/mu       # Reynolds Number

        print "Initial Reynolds Number -->", Re
        dp = 0.0149*dx*ncells/(2*r)*grid.v_r(0)**2/2*rho
        print "dp:", dp
        print "dp/dx:", dp/ncells
        print ((grid.p(1)+x[2])-(grid.p(0)+x[0]))/(grid.dx_rh(0))
        for i in range(ncells):
            newp = (grid.p(i) + x[2*i])
            print "NEW P --> ", newp
            grid.set_p(i, newp)
            
            newv = (grid.v_rh(i) + x[2*i + 1])
            print "\t\tNEW v --> ", newv
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
