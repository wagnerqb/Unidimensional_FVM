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

        for i in range(cpoints-1):
            A = grid.A(i)
            A_r = grid.A_r(i)
            A_rh = grid.A_rh(i)
            rho = grid.rho(i)
            rho_r = grid.rho_r(i)
            rho_rh = grid.rho_rh(i)
            mu = grid.mu(i)
            mu_r = grid.mu_r(i)
            dx = grid.dx(i)
            dx_r = grid.dx_r(i)
            v = grid.v(i)
            v_r = grid.v_r(i)

            ac_k = A*rho*v
            ad_k = (A*mu)/dx
            ac_kp1 = A_r*rho_r*v_r
            ad_kp1 = A_r*mu_r/dx_r

            ##  Velocidades Positivas  ##
            if (grid.v_rh(i) >= 0):
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
                v_BC = v*A_rh/A
#                v_matrix[i][i] =  -a_vl*A_kph/A + 0.5*v_BC*rho*A*(A_kph/A)**2 + ac_kp1 + ad_k +ad_kp1
                v_matrix[i][i] =  0.5*v_BC*rho*A*(A_rh/A)**2 + ac_kp1 + ad_k +ad_kp1
                v_matrix[i][i+1] = a_vr
            elif i == cpoints-2:
                v_matrix[i][i-1] = a_vl
#                v_matrix[i][i] = a_vc #- a_vr*A_rh/A_lh
                v_matrix[i][i] = ad_k + A_rh*rho_rh*grid.v_rh(i)
                print 
            else:
                v_matrix[i][i-1] = a_vl
                v_matrix[i][i] = a_vc
                v_matrix[i][i+1] = a_vr

        return v_matrix

    def build_coef_vector_v(self, grid):
        #classe que constroi o vetor que será resolvido pelo solver
        n = len(grid.cells)
        b = np.zeros(n-1)

        for i in range(n-1):
            A_c = self.center_scheme(grid.A(i), grid.A(i+1))
            b[i] = -A_c*(grid.p(i+1)-grid.p(i))

        #Condicao de Contorno a Esquerda
        A = grid.A(0)
        A_rh = grid.A_rh(0)
        mu = grid.mu(0)
        rho = grid.rho(0)
        dx = grid.dx(0)
        v = grid.v(0)
        v_BC = v*A_rh/A
        ac_k = A*rho*v_BC
        ad_k = (A*mu)/dx
        b[0] += (ac_k + ad_k)*v*A_rh/A

        #Condicao de Contorno a Direita
        #Automatica

        return b
        
    def build_matrix_p(self, grid, v_prev):
        "Classe que constrói a matrz de p' (previsão de velocidade)"
        n = len(v_prev)
        p_matrix = np.zeros((n-1, n-1))
        d = [0]*4
        
        for i in range(n):
            v_r = grid.v_r(i)
            v_rh = grid.v_rh(i)
            v_lh = grid.v_lh(i)
            A_rh = grid.A_rh(i)
            A_lh = grid.A_lh(i)
            A_r = grid.A_r(i)
            A = grid.A(i)
            rho = grid.rho(i)
            rho_r = grid.rho_r(i)
            rho_lh = grid.rho_lh(i)
            rho_rh = grid.rho_rh(i)
            
            if (i == 0):          
#                Termo da pressão de estagnaçao. RETIRAR FUTURO
                t = rho_lh*(v_rh*A_rh/A)*A_lh*0.5*(A_rh*A_rh)/(A*A)
                d[i] = A_rh/(A_r*rho_r*v_r + t) 
            elif (i == 3):
#                Termo de Direita. REVISAR NO FUTURO
                d[i] = A_rh/(A_rh*rho_rh*v_rh)
                
            else:
                d[i] = A_rh/(A_r*rho_r*v_r)
            
        print "d --> ", d
        print
            
        for i in range(1, n):
            
            v_r = grid.v_r(i)
            v_rh = grid.v_rh(i)
            v_lh = grid.v_lh(i)
            A_rh = grid.A_rh(i)
            A_lh = grid.A_lh(i)
            A_r = grid.A_r(i)
            A = grid.A(i)
            rho = grid.rho(i)
            rho_r = grid.rho_r(i)
            rho_lh = grid.rho_lh(i)
            rho_rh = grid.rho_rh(i)
            
             #Filling Matrix
            if i == 1:
                al = A_lh*rho_lh*d[i-1]
                ar = A_rh*rho_rh*d[i]
                ac = - al - ar
                
                p_matrix[i-1][i-1] = ac
                p_matrix[i-1][i] = ar
            elif i == (n-1):
                al = A_lh*rho_lh*d[i-1]
                ar = A_rh*rho_rh*d[i]
                ac = - al - ar
                
                p_matrix[i-1][i-2] = al
                p_matrix[i-1][i-1] = ac
            else:
                al = A_lh*rho_lh*d[i-1]
                ar = A_rh*rho_rh*d[i]
                ac = - al - ar
                
                p_matrix[i-1][i-2] = al
                p_matrix[i-1][i-1] = ac
                p_matrix[i-1][i] = ar
        
        return p_matrix

            
        

    #Funções de Interpolacao
    def center_scheme(self, leftprop, rightprop):
        "Classe interpola utilizando diferencas centrais."
        return (rightprop + leftprop)/2

if __name__ == '__main__':

    from GridFluid import *

    gridteste = GridFluid(0, 10, 0, 0, 0)
    for i in range(5):
        A = 0.5-0.1*(i)
        Ai = 0.45 -0.1*i
        gridteste.add_cell(A=A, dx=0.5, rho=1, mu=0, v_r=1/Ai, p=10-2.5*(i))

    modelCD = Model_SIMPLE_UDS()

    A = modelCD.build_matrix_v(gridteste)
    b = modelCD.build_coef_vector_v(gridteste)
    v = (np.matrix(A).I*np.matrix(b).T).A1
    
    print "\n A: \n", A
    print "\n b: \n", b
    print "\n p: \n", v
    print
    
    Ap = modelCD.build_matrix_p(gridteste, v)    
#    b = []
    
    print "\n Ap: \n", Ap
    print
    

