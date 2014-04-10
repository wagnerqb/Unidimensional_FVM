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
    "Classe de Modelo."

    def build_matrix(self, grid):
        "Classe que constrói a matrz que será resolvida pelo solver"
        cpoints = len(grid.cells)
        A = np.zeros((cpoints, cpoints))

        for i in range(cpoints):
#            k_lh = center_scheme()
        CONTINUAR DAQUI
            AK_left = grid.k_lh(i)*grid.A_lh(i)
            l_k = 2*AK_left/(grid.dx(i) + grid.dx(i-1))

            AK_right = grid.k_rh(i) * grid.A_rh(i)
            r_k = 2*AK_right/(grid.dx(i+1) + grid.dx(i))

            if i == 0:
                A[i][i] = - l_k - r_k
                A[i][i+1] = r_k
            else:
                if i == cpoints-1:
                    A[i][i-1] = l_k
                    A[i][i] = - l_k - r_k
                else:
                    A[i][i-1] = l_k
                    A[i][i] = - l_k - r_k
                    A[i][i+1] = r_k
        return A

    def build_coef_vector(self, grid):
        #classe que constrói o vetor que será resolvido pelo solver
        cpoints = len(grid.cells)
        B = np.zeros(cpoints)

        for i in range(cpoints):
            #Termos Fonte
            B[i] = - grid.dx(i)*grid.A(i)*grid.Source(i)

        #Condição de Contorno Esquerda
        AK_left = grid.k_lh(0)*grid.A_lh(0)
        l_k = 2*AK_left/(grid.dx(0) + grid.dx(-1))
        B[0] = B[0] - l_k*grid.get_phi(-1)

        #Condição de Contorno Direita
        AK_right = grid.k_rh(cpoints - 1) * grid.A_rh(cpoints - 1)
        r_k = 2*AK_right/(grid.dx(cpoints) + grid.dx(cpoints - 1))
        B[cpoints - 1] = B[cpoints - 1] - r_k*grid.get_phi(cpoints)

        return B
        
    #Funções de Interpolação
    def center_scheme(self, leftprop, rightprop):
        "Classe interpola utilizando diferencas centrais." 
        return (rightprop + leftprop)/2

    #Funções para Calcular Derivadas
    def right_derivative(self, leftprop, rightprop, leftx, rightx):
        return (rightprop - leftprop)/(rightx/2 + leftx/2)
        
    def dphi_dx_lh(self, index, grid):
        "dphi/dx left half"
        if index == 0:
            return 0
        cell = grid[index]
        l_cell = grid[index-1]
        return self.right_derivative(l_cell.phi, cell.phi, l_cell.dx,
                                                  cell.dx)
                                                  

    def dphi_dx_rh(self, index, grid):
        "dT/dx right half"
        if (index == len(self.cells) - 1):
            return 0
        cell = grid[index]
        rightcell = grid[index+1]
        return self.right_derivative(cell.phi, rightcell.phi,
                                                  cell.dx, rightcell.dx)


if __name__ == '__main__':

    from Interpolation import *
    modelteste = Model()
    gridteste = Grid(Interpolation(), 100, 500, 0)

    for i in range(5):
        gridteste.add_cell(10e-3, 1000, 100, 0.1)

    A = modelteste.build_matrix(gridteste)
    print A

    b = modelteste.build_coef_vector(gridteste)
    print b
