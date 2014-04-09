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
        
    A = 1
    kappa = 2
    T = 3
    deltax = 4
#    ncells = input("Digite o número de Células")
    ncells = 5
    
    grid = []
    
    for i in range(ncells) :
        grid.addcell(A, kappa, T, deltax)
        
#    for g in grid :    
#        print g.A, g.kappa, g.T, g.deltax
#        print g.kappa

            
        
            
        
if __name__ == '__main__':
    
    run()
    
