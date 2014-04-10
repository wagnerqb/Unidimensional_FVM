# -*- coding: iso-8859-1 -*-
"""
@author: Wagner Queiroz Barros
@date:   Wed Apr 09 15:35:41 2014
@email:  wagnerqb@gmail.com
@brief:  Class Storing all models using to discretize PDEs
"""

from __future__ import division
from Grid import *
import numpy as np


class Model():
    "Classe de Modelo"
    
    #Funções de Derivada
    def RightDerivative(self, leftprop, rightprop, leftx, rightx):
        return (rightprop - leftprop)/(rightx/2 + leftx/2)
        
    #Funções de Interpolação
    def CenterScheme(self, leftprop, rightprop):
        #Classe interpola utilizando diferencas centrais
        return (rightprop + leftprop)/2
        
    def BuildMatrix(self, grid):
        
        cpoints = len(grid.grid)
        A = np.zeros((cpoints,cpoints))
        
        AK_right = grid.RightKappa(0)*grid.RightArea(0)
        r_k = 2*AK_right/(grid.grid[1].deltax + grid.grid[0].deltax)
        c_k = - r_k -2*grid.grid[0].A*grid.grid[0].kappa/grid.grid[0].deltax 
        
        A[0][0] = c_k
        A[0][1] = r_k
        for  i in range(cpoints-2):
            
            AK_left = grid.LeftKappa(i+1)*grid.LeftArea(i+1)
            l_k = 2*AK_left/(grid.grid[i+1].deltax + grid.grid[i].deltax)
            AK_right = grid.RightKappa(i+1)*grid.RightArea(i+1)
            r_k = 2*AK_right/(grid.grid[i+2].deltax + grid.grid[i+1].deltax)
            
            A[i+1][i] = l_k
            A[i+1][i+1] = - l_k - r_k
            A[i+1][i+2] = r_k
            
        
        AK_left = grid.LeftKappa(cpoints-1)*grid.LeftArea(cpoints-1)
        l_k = 2*AK_left/(grid.grid[cpoints-1].deltax + grid.grid[cpoints-2].deltax)         
        c_k = - l_k -2*grid.grid[cpoints-1].A*grid.grid[cpoints-1].kappa/grid.grid[cpoints-1].deltax
        
        A[cpoints-1][cpoints-2] = l_k
        A[cpoints-1][cpoints-1] = c_k
        
        return A

if __name__ == '__main__':
    
    modelteste = Model()
    gridteste = Grid(Model())
    
    for i in range(5) :
        gridteste.Addcell(10e-3, 1000, 100, 0.1)
    
    B = modelteste.BuildMatrix(gridteste)
    print B