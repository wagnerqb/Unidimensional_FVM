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
    def Addcell(self, A, kappa, T, deltax):
        self.grid.append(Cell(A, kappa, T, deltax))

    def LeftdT_dx(self, index):
        cell = self.grid[index]
        leftcell = self.grid[index-1]
        return self.model.RightDerivative(leftcell.T, cell.T, leftcell.deltax,
                                          cell.deltax)
                
    def RightdT_dx(self, index):
        cell = self.grid[index]
        rightcell = self.grid[index+1]
        return self.model.RightDerivative(cell.T, rightcell.T, cell.deltax,
                                          rightcell.deltax)
    
    def LeftKappa(self, index):
        cell = self.grid[index]
        leftcell = self.grid[index-1]
        return self.model.CenterScheme(leftcell.kappa, cell.kappa)
        
    def RightKappa(self, index):
        cell = self.grid[index]
        rightcell = self.grid[index+1]
        return self.model.CenterScheme(cell.kappa, rightcell.kappa)
        
    def LeftArea(self, index):
        cell = self.grid[index]
        leftcell = self.grid[index-1]
        return self.model.CenterScheme(leftcell.A, cell.A)
        
    def RightArea(self, index):
        cell = self.grid[index]
        rightcell = self.grid[index+1]
        return self.model.CenterScheme(cell.A, rightcell.A)

if __name__ == '__main__':

    gridteste = Grid(Model())

    gridteste.Addcell(1, 2, 3, 4)
    gridteste.Addcell(5, 6, 7, 15)
    gridteste.Addcell(9,10,11,12)

    print gridteste.grid[0].T
    print gridteste.grid[1].T

    print gridteste.grid[0].deltax
    print gridteste.grid[1].deltax

    print gridteste.RightdT_dx(0)
    
    print gridteste.grid[1].T
    print gridteste.grid[2].T

    print gridteste.grid[1].deltax
    print gridteste.grid[2].deltax
    
    print gridteste.RightdT_dx(1)
    print gridteste.LeftdT_dx(2)
    