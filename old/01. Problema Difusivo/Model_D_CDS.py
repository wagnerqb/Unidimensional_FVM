# -*- coding: iso-8859-1 -*-
"""
@author: Wagner Queiroz Barros
@date:   Wed Apr 09 15:35:41 2014
@email:  wagnerqb@gmail.com
@brief:  Class used to discretize Difuse PDEs using CDS Scheme
"""

from __future__ import division
from GridD import *
import numpy as np


class Model_D_CDS():
    """Classe de Modelo para Equa��es Difusivas, utilizando o esquema
    CDS Central Discretization Scheme ."""

    def build_matrix(self, grid):
        "Classe que constr�i a matrz que ser� resolvida pelo solver"
        cpoints = len(grid.cells)
        A = np.zeros((cpoints, cpoints))

        for i in range(cpoints):
            k_lh = self.center_scheme(grid.k(i-1), grid.k(i))
            k_rh = self.center_scheme(grid.k(i), grid.k(i+1))
            A_lh = self.center_scheme(grid.A(i-1), grid.A(i))
            A_rh = self.center_scheme(grid.A(i), grid.A(i+1))

            l_k = - 2*k_lh*A_lh/(grid.dx(i) + grid.dx(i-1))
            r_k = - 2*k_rh*A_rh/(grid.dx(i+1) + grid.dx(i))

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
        #classe que constr�i o vetor que ser� resolvido pelo solver
        cpoints = len(grid.cells)
        B = np.zeros(cpoints)

        for i in range(cpoints):
            #Termos Fonte
            B[i] = grid.dx(i)*grid.A(i)*grid.Source(i)

        #Condi��o de Contorno Esquerda
        k_lh_lbc = self.center_scheme(grid.k(-1), grid.k(0))
        A_lh_lbc = self.center_scheme(grid.A(-1), grid.A(0))

        l_k = - 2*k_lh_lbc*A_lh_lbc/(grid.dx(0) + grid.dx(-1))
        B[0] = B[0] - l_k*grid.phi(-1)

        #Condi��o de Contorno Direita
        k_rh_rbc = self.center_scheme(grid.k(cpoints - 1), grid.k(cpoints))
        A_rh_rbc = self.center_scheme(grid.A(cpoints - 1), grid.A(cpoints))

        r_k = - 2*k_rh_rbc*A_rh_rbc/(grid.dx(cpoints) + grid.dx(cpoints - 1))
        B[cpoints - 1] = B[cpoints - 1] - r_k*grid.phi(cpoints)

        return B

    #Fun��es de Interpola��o
    def center_scheme(self, leftprop, rightprop):
        "Classe interpola utilizando diferencas centrais."
        return (rightprop + leftprop)/2

    #Fun��es para Calcular Derivadas
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

    modelteste = Model_D_CDS()
    gridteste = GridD(100, 500, 0)

    for i in range(5):
        gridteste.add_cell(10e-3, 1000, 0.1, 100)

    A = modelteste.build_matrix(gridteste)
    print A
    print

    b = modelteste.build_coef_vector(gridteste)
    print b
