# -*- coding: iso-8859-1 -*-
"""
@author: Wagner Queiroz Barros
@date:   Wed Apr 09 15:06:58 2014
@email:  wagnerqb@gmail.com
@brief:  Main class, used to run simulation
"""

from __future__ import division
import numpy as np
from Grid import *
from Interpolation import *
from Model import *



def run():

    temp_A = 1
    temp_kappa = 0.5
    temp_T = 0
    temp_deltax = 0.004

    temp_LBC = 100
    temp_RBC = 200
    
    ncells = 5
    temp_source = 1000000

    grid = Grid(Interpolation(), temp_LBC, temp_RBC, temp_source)
    model = Model()

    for i in range(ncells):
        grid.add_cell(temp_A, temp_kappa, temp_T, temp_deltax)

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
