# -*- coding: iso-8859-1 -*-
"""
@author: Wagner Queiroz Barros
@date:   Wed Apr 09 14:53:01 2014
@email:  wagnerqb@gmail.com
@brief:  Grid class, used to store cell data
"""
from __future__ import division
from Model import *
from Cell import *


class Grid():
    "Classe Celula."

    def __init__(self, model):
        #Atributos
        self.grid = []
        self.model = model

    #Métodos
    def addcell(self, A, kappa, T, deltax):
        self.grid.append(Cell(A, kappa, T, deltax))

    def rightdT_dx(self, index):
        cell = self.grid[index]
        cellRight = self.grid[index+1]
        return self.model.RightDerivative(cell.T, cellRight.T, cell.deltax,
                                          cellRight.deltax)

    def leftdT_dx(self, index):
        cell = self.grid[index]
        leftcell = self.grid[index-1]
        return self.model.RightDerivative(leftcell.T, cell.T, leftcell.deltax,
                                          cell.deltax)                 

if __name__ == '__main__':

    gridteste = Grid(Model())

    gridteste.addcell(1, 2, 3, 4)
    gridteste.addcell(5, 6, 7, 15)
    gridteste.addcell(9,10,11,12)

    print gridteste.grid[0].T
    print gridteste.grid[1].T

    print gridteste.grid[0].deltax
    print gridteste.grid[1].deltax

    print gridteste.rightdT_dx(0)
    
    print gridteste.grid[1].T
    print gridteste.grid[2].T

    print gridteste.grid[1].deltax
    print gridteste.grid[2].deltax
    
    print gridteste.rightdT_dx(1)
    print gridteste.leftdT_dx(2)
    