# -*- coding: iso-8859-1 -*-
"""
@author: Wagner Queiroz Barros
@date:   Tue Apr 15 10:55:33 2014
@email:  wagnerqb@gmail.com
@brief:  Class used to discretize Navier-Stokes using SIMPLE Method, CDS for P
and CDS for v. 
"""

from __future__ import division
import numpy as np

class Model_SIMPLE_CDS():
    """Classe de Modelo para Equações de Navier Stokes, utilizando
    o modelo SIMPLE para o acoplamento pressão velocidade, e o
    esquema CDS Central Discretization Scheme para interpolar a velocidade."""

    def build_matrix_v(self, grid):
        "Classe que constrói a matrz de v* (previsão de velocidade)"
        cpoints = len(grid.cells)
        v_matrix = np.zeros((cpoints, cpoints))

        for i in range(cpoints):

            A = grid.A(i)
            A_r = grid.A(i+1)
            rho = grid.rho(i)
            rho_r = grid.rho(i+1)
            mu = grid.mu(i)
            mu_r = grid.mu(i+1)
            dx = grid.dx(i)
            dx_r = grid.dx(i+1)
            v = self.center_scheme(grid.v_r(i-1), grid.v_r(i))
            v_r = self.center_scheme(grid.v_r(i), grid.v_r(i+1))

            # Constructing matrix coefficients
            ac_vl = 0.5*A*rho*v
            ad_vl = (A*mu)/dx

            ac_vr = 0.5*A_r*rho_r*v_r
            ad_vr = (A_r*mu_r)/dx_r

            a_vl = - ac_vl - ad_vl
            a_vc = ac_vr + ad_vr - ac_vl + ad_vl
            a_vr = ac_vr - ad_vr

            A_lh = self.center_scheme(grid.A(i-1), grid.A(i))
            A_rh = self.center_scheme(grid.A(i), grid.A(i+1))

            #Filling Matrix
            if i == 0:
                v_matrix[i][i] = a_vc - a_vl*A_lh/A_rh
                v_matrix[i][i+1] = a_vr
            else:
                if i == cpoints-1:
                    v_matrix[i][i-1] = a_vl
                    v_matrix[i][i] = a_vc - a_vr*A_rh/A_lh
                else:
                    v_matrix[i][i-1] = a_vl
                    v_matrix[i][i] = a_vc
                    v_matrix[i][i+1] = a_vr

        return v_matrix

    def build_coef_vector_v(self, grid):
        #classe que constrói o vetor que será resolvido pelo solver
        n = len(grid.cells)
        b = np.zeros(n)

        for i in range(n-1):
            A_c = self.center_scheme(grid.A(i), grid.A(i+1))
            b[i] = -A_c*(grid.p(i+1)-grid.p(i))

        return b

    #Funções de Interpolação
    def center_scheme(self, leftprop, rightprop):
        "Classe interpola utilizando diferencas centrais."
        return (rightprop + leftprop)/2

if __name__ == '__main__':

    from GridFluid import *

    gridteste = GridFluid(10, 0, 1)
    for i in range(5):
        A = 0.5-0.05*(i+1)
        gridteste.add_cell(A=A, dx=0.5, rho=1, mu=0, v_r=1/A, p=10-2.5*(i+1))

    modelCD = Model_SIMPLE_CDS()

    A = modelCD.build_matrix_v(gridteste)
    b = modelCD.build_coef_vector_v(gridteste)
    print "\nA:\n", A
    print "\nb:\n", b
    print "\nv:\n", (np.matrix(A).I*np.matrix(b).T).A1
    print

