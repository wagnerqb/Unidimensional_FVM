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

def run():
        
    temp_A = 10e-3
    temp_kappa = 1000
    temp_T = 100
    temp_deltax = 0.1
#    ncells = input("Digite o número de Células")
    ncells = 5
    
    grid = Grid(Model())
    model = Model()
    
    for i in range(ncells) :
        grid.Addcell(temp_A, temp_kappa, temp_T, temp_deltax)
 
    A = model.BuildMatrix(grid)    
    print A
    
    
if __name__ == '__main__':
    
    run()
    
        
    
