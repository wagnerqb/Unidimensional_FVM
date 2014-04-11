# -*- coding: iso-8859-1 -*-
"""
@author: Wagner Queiroz Barros
@date:   Fri Apr 11 11:14:14 2014
@email:  wagnerqb@gmail.com
@brief:  GridCD class, used to store cell data in Difusive/Convective Problems
"""
from __future__ import division

from GridD import *

class GridCD(GridD):
    
    def __init__(self, lbc, rbc, Source):
        
        GridD.__init__(self, lbc, rbc, Source)
        
    #Métodos
    def add_cell(self, A, k, dx, phi, rho, v):
        "Adiciona uma celula no grid."
        self.cells.append(CellCD(A, k, dx, phi, rho, v))
        
    def rho(self, index):
        "densidade no Centro da célula"
        if (index < 0):
            return self[0].rho

        if (index > (len(self.cells)-1)):
            return self[(len(self.cells)-1)].rho

        return self[index].rho
        
    def v(self, index):
        "velocidade no Centro da célula"
        if (index < 0):
            return self[0].v

        if (index > (len(self.cells)-1)):
            return self[(len(self.cells)-1)].v

        return self[index].v


if __name__ == '__main__':

    #from Model_D_CDS import *

    ncells = 5

    grid = GridCD(100, 500, 0)
    #model = Model_D_CDS()

    for i in range(ncells):
        grid.add_cell(10e-3, 1000, 0.1, 100, 250, 0.1)

    print grid.A(-1)
    print grid.A(2)
    print grid.A(5)
    print
    
    print grid.k(-1)
    print grid.k(1)
    print grid.k(6)
    print
    
    print grid.dx(-1)
    print grid.dx(2)
    print grid.dx(5)
    print
    
    print grid.phi(-1)
    print grid.phi(3)
    print grid.phi(5)
    print
    
    print grid.rho(-1)
    print grid.rho(2)
    print grid.rho(5)
    print
    
    print grid.v(-1)
    print grid.v(2)
    print grid.v(6)
