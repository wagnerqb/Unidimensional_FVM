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

    def build_v_matrix(self, grid):
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

            A_vl = - ac_vl - ad_vl
            A_vc = ac_vr + ad_vr - ac_vl + ad_vl
            A_vr = ac_vr - ad_vr

            #Filling Matrix
            if i == 0:
                v_matrix[i][i] = A_vc
                v_matrix[i][i+1] = A_vr
            else:
                if i == cpoints-1:
                    v_matrix[i][i-1] = A_vl
                    v_matrix[i][i] = A_vc
                else:
                    v_matrix[i][i-1] = A_vl
                    v_matrix[i][i] = A_vc
                    v_matrix[i][i+1] = A_vr
        return v_matrix

    #Funções de Interpolação
    def center_scheme(self, leftprop, rightprop):
        "Classe interpola utilizando diferencas centrais."
        return (rightprop + leftprop)/2

if __name__ == '__main__':

    from GridFluid import *

    gridteste = GridFluid(100, 500, 0)
    for i in range(5):
        gridteste.add_cell(10e-3, 1000, 0.1, 100, 1, 2)

    modelCD = Model_SIMPLE_CDS()

    A = modelCD.build_v_matrix(gridteste)
    print A
    print

    print "TESTE"
