# -*- coding: iso-8859-1 -*-
"""
@author: Wagner Queiroz Barros
@date:   Fri Apr 11 13:24:24 2014
@email:  wagnerqb@gmail.com
@brief:  Simulator for Navier-Stokes problem, using SIMPLE UDS method
"""


from __future__ import division
import numpy as np
from GridFluid_SIMPLE import *
from Model_SIMPLE_UDS import *


def run():

    # Pipe properties
    r = 0.1                 # Pipe radius
    dx = 0.1                # Discretization lenght delta_x
    ncells = 10              # Number of cells in domai

    # Fluid properties
    rho = 1                 # Fluid density
    mu = 0.                # Fluid viscosity
    msrc = 0                # Mass Source term per volume unity

    # Initial properties
    v_ini = 14.5               # Initial Condition for v
    p_ini = 1               # Initial Condition for p

    # Boundary Condition
    # Left Boundary Condtion
    lbc_t = 1               # LBC Type (0 - Pressure / 1 - Velocity)
    lbc = 15                # LBC Value

    # Right Boundary Condtion
    rbc_t = 0               # LBC Type (0 - Pressure / 1 - Velocity)
    rbc = 10                # LBC Value

    # Numerical Properties
    p_relax = 0.4           # Pressure Relaxation Factor
    v_relax = 0.4           # Velocity Relaxation Factor

    # Calculated Values
    Ar = np.pi*r*r           # Pipe Area
#    Re = rho*v_ini*2*r/mu       # Reynolds Number

    # Creating Grid
    grid = GridFluid_SIMPLE(lbc_t, lbc, rbc_t, rbc, msrc)

    # Creating Model
    model = Model_SIMPLE_UDS(p_relax, v_relax)

    for i in range(ncells):
        grid.add_cell(Ar, dx, rho, mu, v_ini, p_ini)

#    print "Initial Reynolds Number -->", Re
    print "Pipe Area  -->", Ar
    print

    for j in range(50):
        print "Iteration ", (j+1)

        #Construindo Matrix Previsão de Velocidade do Sistema
        Av = np.matrix(model.build_matrix_v(grid))
        bv = np.matrix(model.build_coef_vector_v(grid))
        xv = (Av.I*bv.T).A1

#        print "VELOCITY PREDICTION MATRIX"
#        print Av
#        print
#
#        print "VELOCITY PREDICTION VECTOR"
#        print bv
#        print
##
#        print "VELOCITY PREDICTION"
#        print xv
#        print

#        Atribuindo os resultados de velocidade
        model.set_guess_v(grid, xv)

#        for i in range(ncells):
#            print "NEW v --> ", grid.v_rh(i)

        #Construindo a Matriz de correção de pressão do sistema
        Ap = np.matrix(model.build_matrix_p(grid))
        Bp = np.matrix(model.build_coef_vector_p(grid))
        Rp = (np.matrix(Ap).I*np.matrix(Bp).T).A1

        print "PRESSURE CORRECTION MATRIX"
        print Ap
        print

        print "PRESSURE CORRECTION VECTOR"
        print Bp
        print

        print "PRESSURE CORRECTION"
        print Rp
        print

        model.correct_p(grid, Rp)
        model.correct_v_rh(grid, Rp)

        print " --- NEW p --- "
        for i in range(len(grid.cells)):
            print grid.p(i)
        print

        print " --- NEW v_r --- "
        for i in range(len(grid.cells)):
            print grid.v_r(i)
        print


if __name__ == '__main__':

    run()
