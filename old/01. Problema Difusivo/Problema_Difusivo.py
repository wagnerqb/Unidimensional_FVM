# -*- coding: iso-8859-1 -*-
"""
@author: Wagner Queiroz Barros
@date:   Wed Apr 09 15:06:58 2014
@email:  wagnerqb@gmail.com
@brief:  Simulator for difusive problems, described in Chapter 1 of report
"""

from __future__ import division
import numpy as np
from GridD import *
from Model_D_CDS import *


def run():

    A = 1                   # Pipe Area
    k = 0.5                 # Conductive factor kappa
    T = 0                   # Inicial Condition for Temperature
    dx = 0.004              # Discretization lenght delta_x

    LBC = 100               # Left Boundary Condition
    RBC = 200               # Right Boundary Condition

    ncells = 5              # Number of cells in domain
    source = 1000000        # Source term per volume unity

    grid = GridD(LBC, RBC, source)
    model = Model_D_CDS()

    for i in range(ncells):
        grid.add_cell(A, k, dx, T)

    #Resolvendo o Sistema
    A = np.matrix(model.build_matrix(grid))
    b = np.matrix(model.build_coef_vector(grid))
    x = (A.I*b.T).A1

    #Atribuindo os resultados
    for i in range(ncells):
        grid.set_phi(i, x[i])

    grid.print_phi()

    print A
    print b
#    print grid.get_all_x()
#    print grid.get_all_T()
    grid.plot_T()


if __name__ == '__main__':

    run()
