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
        v_matrix = np.zeros((cpoints-1, cpoints-1))
        
#        for i in range(5):
#                print grid.v_r(i),

        for i in range(cpoints-1):
            A = grid.A(i)
            A_r = grid.A(i+1)
            A_kph = self.center_scheme(grid.A(i), grid.A(i+1))
            rho_kph = self.center_scheme(grid.rho(i), grid.rho(i+1))
            rho = grid.rho(i)
            rho_r = grid.rho(i+1)
            mu = grid.mu(i)
            mu_r = grid.mu(i+1)
            dx = grid.dx(i)
            dx_r = grid.dx(i+1)
            v = self.center_scheme(grid.v_r(i-1), grid.v_r(i))
            v_r = self.center_scheme(grid.v_r(i), grid.v_r(i+1))

            ac_k = A*rho*v
            ad_k = (A*mu)/dx
            ac_kp1 = A_r*rho_r*v_r
            ad_kp1 = A_r*mu_r/dx_r

            ##  Velocidades Positivas  ##
            if (grid.v_r(i) >= 0):
                a_vl = -ac_k - ad_k
                a_vc = ac_kp1 + ad_kp1 + ad_k
                a_vr = -ad_kp1

            ##  Velocidades Negativas  ##
            else:
                a_vl = -ad_k
                a_vc = ad_kp1 - ac_k + ad_k
                a_vr = ac_kp1 - ad_kp1


            #Filling Matrix
            if i == 0:
                v_BC = v*A_kph/A
#                v_matrix[i][i] =  -a_vl*A_kph/A + 0.5*v_BC*rho*A*(A_kph/A)**2 + ac_kp1 + ad_k +ad_kp1
                v_matrix[i][i] =  0.5*v_BC*rho*A*(A_kph/A)**2 + ac_kp1 + ad_k +ad_kp1
                v_matrix[i][i+1] = a_vr
            elif i == cpoints-2:
                v_matrix[i][i-1] = a_vl
#                v_matrix[i][i] = a_vc #- a_vr*A_rh/A_lh
                v_matrix[i][i] = ad_k + A_kph*rho_kph*grid.v_r(i)
                print 
            else:
                v_matrix[i][i-1] = a_vl
                v_matrix[i][i] = a_vc
                v_matrix[i][i+1] = a_vr

        return v_matrix

    def build_coef_vector_v(self, grid):
        #classe que constrói o vetor que será resolvido pelo solver
        n = len(grid.cells)
        b = np.zeros(n-1)

        for i in range(n-1):
            A_c = self.center_scheme(grid.A(i), grid.A(i+1))
            b[i] = -A_c*(grid.p(i+1)-grid.p(i))

        #Condicao de Contorno a Esquerda
        A = grid.A(0)
        A_kph = self.center_scheme(grid.A(0), grid.A(1))
        mu = grid.mu(0)
        rho = grid.rho(0)
        dx = grid.dx(0)
        v = self.center_scheme(grid.v_r(-1), grid.v_r(0))
        v_BC = v*A_kph/A
        ac_k = A*rho*v_BC
        ad_k = (A*mu)/dx
        b[0] += (ac_k + ad_k)*v*A_kph/A

        #Condicao de Contorno a Direita
        #Automatica

        return b
        
    def build_matrix_p(self, grid, v):
        "Classe que constrói a matrz de p' (previsão de velocidade)"
        n = len(v)
        A = np.zeros((n-1, n-1))
        
        for i in range(1, n):
            v_kp1 = self.center_scheme(v[i+2], v[i+1])
            A_kph = self.center_scheme(grid.A(i), grid.A(i+1))
            v_km1 = self.center_scheme(v[i], v[i+1])
            A_kmh = self.center_scheme(grid.A(i), grid.A(i-1))
            print A_kmh/v_k
            print 
            
        

    #Funções de Interpolação
    def center_scheme(self, leftprop, rightprop):
        "Classe interpola utilizando diferencas centrais."
        return (rightprop + leftprop)/2

if __name__ == '__main__':

    from GridFluid import *

    gridteste = GridFluid(0, 10, 0, 0, 1)
    for i in range(5):
        A = 0.5-0.1*(i)
        Ai = 0.45 -0.1*i
        gridteste.add_cell(A=A, dx=0.5, rho=1, mu=0, v_r=1/Ai, p=10-2.5*(i))

    modelCD = Model_SIMPLE_UDS()

    A = modelCD.build_matrix_v(gridteste)
    b = modelCD.build_coef_vector_v(gridteste)
    v = (np.matrix(A).I*np.matrix(b).T).A1
    
    A = modelCD.build_matrix_p(gridteste, v)    
    b = []
    
    print "\nA:\n", A
    print "\nb:\n", 
    print "\np:\n",

