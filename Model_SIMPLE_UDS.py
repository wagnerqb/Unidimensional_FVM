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


class Model_SIMPLE_UDS():
    """Classe de Modelo para Equações de Navier Stokes, utilizando
    o modelo SIMPLE para o acoplamento pressão velocidade, e o
    esquema UDS Central Discretization Scheme para interpolar a velocidade."""

    def build_matrix_v(self, grid):
        "Classe que constrói a matrz de v* (previsão de velocidade)"
        cpoints = len(grid.cells)
        v_matrix = np.zeros((cpoints, cpoints))

        for i in range(cpoints):
            A = grid.A(i)
            A_r = grid.A_r(i)
            rho = grid.rho(i)
            rho_r = grid.rho_r(i)
            mu = grid.mu(i)
            mu_r = grid.mu_r(i)
            dx = grid.dx(i)
            dx_r = grid.dx_r(i)
            v = grid.v(i)
            v_r = grid.v_r(i)

            ac_k = A*rho*v
            ad_k = (A*mu)/dx
            ac_kr = A_r*rho_r*v_r
            ad_kr = A_r*mu_r/dx_r

            ##  Velocidades Positivas  ##
            if (grid.v_rh(i) >= 0):
                a_vl = - ac_k - ad_k
                a_vc = ac_kr + ad_kr + ad_k
                a_vr = -ad_kr

            ##  Velocidades Negativas  ##
            else:
                a_vl = - ad_k
                a_vc = - ac_k + ad_k + ad_kr
                a_vr = ac_kr - ad_kr

            #Filling Matrix
            if i == 0:
                v_matrix[i][i] = a_vc
                v_matrix[i][i+1] = a_vr

            elif i == (cpoints - 1):
                v_matrix[i][i-1] = a_vl
                v_matrix[i][i] = a_vc
                print

            else:
                v_matrix[i][i-1] = a_vl
                v_matrix[i][i] = a_vc
                v_matrix[i][i+1] = a_vr

        return v_matrix

    def build_coef_vector_v(self, grid):
        #classe que constroi o vetor que será resolvido pelo solver
        n = len(grid.cells)
        b = np.zeros(n)

        for i in range(n):
            A_rh = grid.A_rh(i)
            p_r = grid.p_r(i)
            p = grid.p(i)

            b[i] = - A_rh*(p_r - p)

        #Condicao de Contorno a Esquerda
        A0 = grid.A(0)
        rho0 = grid.rho(0)
        dx0 = grid.dx(0)
        mu0 = grid.mu(0)
        v0 = grid.v(0)
        v_lh0 = grid.v_lh(0)

        if (grid.v_r(0) >= 0):
            a_vl = - (A0*rho0*v0) - (A0*mu0)/dx0
        else:
            a_vl = - (A0*mu0)/dx0
        b[0] += - a_vl*v_lh0

        #Condicao de Contorno a Direita
        A_rl = grid.A_r(n-1)
        v_rl = grid.v_r(n-1)
        mu_rl = grid.mu_r(n-1)
        rho_rh = grid.rho_rh(n-1)
        dx_rl = grid.dx_r(n-1)
        v_rrl = grid.v_r(n)

        if (grid.v_r(n-1) >= 0):
            a_vr = - (A_rl*mu_rl)/dx_rl
        else:
            a_vr = (A_rl*rho_rh*v_rl) - (A_rl*mu_rl)/dx_rl
        b[n-1] += - a_vr*v_rrl

        return b

    def build_matrix_p(self, grid):
        "Classe que constrói a matrz de p' (previsão de velocidade)"
        n = len(grid.cells)
        p_matrix = np.zeros((n, n))

        for i in range(n):
            v_r = grid.v_r(i)
            v = grid.v(i)
#            v_rh = grid.v_rh(i)
#            v_lh = grid.v_lh(i)
#            A_l = grid.A_l(i)
            A_rh = grid.A_rh(i)
            A_rrh = grid.A_rh(i+1)
            A_lh = grid.A_lh(i)
            A_r = grid.A_r(i)
            A = grid.A(i)
            rho = grid.rho(i)
            rho_r = grid.rho_r(i)
            rho_lh = grid.rho_lh(i)
            rho_rh = grid.rho_rh(i)
            rho_rrh = grid.rho_rh(i+1)

            if (grid.v_r(i) >= 0):
                al = A_lh*rho_lh*A_lh/(A*rho*v)
                ar = A_rh*rho_rh*A_rh/(A_r*rho_r*v_r)
                ac = - al - ar
            else:
                al = A_rh*rho_rh*A_lh/(A*rho*v)
                ar = A_rrh*rho_rrh*A_rh/(A_r*rho_r*v_r)
                ac = - al - ar

            if (i == 0):
                p_matrix[i][i] = ac
                p_matrix[i][i+1] = ar

            elif (i == (n - 1)):

                p_matrix[n-1][n-2] = al
                p_matrix[n-1][n-1] = ac

            else:
                p_matrix[i][i-1] = al
                p_matrix[i][i] = ac
                p_matrix[i][i+1] = ar

        return p_matrix

    def build_coef_vector_p(self, grid):
        "Classe que constroi o vetor de coeficientes da correcao da pressao."
        #classe que constroi o vetor que será resolvido pelo solver
        n = len(grid.cells)
        bp = np.zeros(n)

        for i in range(n):
            A = grid.A(i)
            A_lh = grid.A_lh(i)
            A_rh = grid.A_rh(i)
            A_rrh = grid.A_rh(i+1)
            A_r = grid.A_r(i)
            rho_lh = grid.rho_lh(i)
            rho_rh = grid.rho_rh(i)
            rho_rrh = grid.rho_rh(i+1)
            dx = grid.dx(i)
            dx_r = grid.dx_r(i)
            v = grid.v(i)
            v_r = grid.v_r(i)
            
            msrc = grid.msrc
            
            if (grid.v_r(i) >= 0):
                RHS = - (A*dx*msrc) + A_rh*rho_rh*v_r
                RHS +=  - A_lh*rho_lh*v
            else:
                RHS = - (A_r*dx_r*msrc) + A_rrh*rho_rrh*v_r
                RHS +=  - A_rh*rho_rh*v
                

            bp[i] = RHS

        #Condicao de Contorno a Esquerda
        #Automatica (REVER)

        #Condicao de Contorno a Direita
        #Automatica (REVER)

        return bp
        
    def correct_p(self, grid, dp):
        "Funcao que corrige os valores de p na celula de acordo com p'."
        n = len(grid.cells)

        for i in range(n):
            grid.set_p(i, grid.p(i) + dp[i])
            
    def correct_v_rh(self, grid, dp):
        "Funcao que corrige as velocidas na celula de acordo com p'."
        n = len(grid.cells)
        
        for i in range(n):
            A_lh = grid.A_lh(i)
            A = grid.A(i)
            rho = grid.rho(i)
            v = grid.v(i)
            
            coef = A_lh/(A*rho*v)
            if (i == 0):
                grid.set_v_rh(i, - coef*dp[i])
            else:
                grid.set_v_rh(i, - coef*(dp[i] - dp[i-1]))

    #Funções de Interpolacao
    def center_scheme(self, leftprop, rightprop):
        "Classe interpola utilizando diferencas centrais."
        return (rightprop + leftprop)/2

if __name__ == '__main__':

    from GridFluid import *

    gridteste = GridFluid(0, 10, 0, 0, 0)
    for i in range(5):
        A = 0.5 - 0.1*(i)
        Ai = 0.45 - 0.1*i
        gridteste.add_cell(A=A, dx=0.5, rho=1, mu=0, v_r=1/Ai, p=10-2.5*(i))

    modelCD = Model_SIMPLE_UDS()

    A = modelCD.build_matrix_v(gridteste)
    b = modelCD.build_coef_vector_v(gridteste)
    v = (np.matrix(A).I*np.matrix(b).T).A1

    print "\n A: \n", A
    print "\n b: \n", b
    print "\n v: \n", v
    print

    Ap = modelCD.build_matrix_p(gridteste)
    Bp = modelCD.build_coef_vector_p(gridteste)
    Rp = (np.matrix(Ap).I*np.matrix(Bp).T).A1

    print "\n Ap: \n", Ap
    print

    print "\n Bp: \n", Bp
    print

    print "\n Rp: \n", Rp
