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


class Model_COUPLE_CDS():
    """Classe de Modelo para Equações de Navier Stokes, utilizando
    o modelo COUPLE para o acoplamento pressão velocidade, e o
    esquema CDS Central Discretization Scheme para interpolar a velocidade."""

    def build_Jacobian_matrix(self, grid):
        "Classe que constrói a matrz de v* (previsão de velocidade)"
        cpoints = len(grid.cells)
        v_matrix = np.zeros((2*cpoints, 2*cpoints))

        for i in range(cpoints):

            A = grid.A(i)
            A_r = grid.A_r(i)
            A_lh = grid.A_lh(i)
            A_rh = grid.A_rh(i)
            rho = grid.rho(i)
            rho_r = grid.rho_r(i)
            rho_lh = grid.rho_lh(i)
            rho_rh = grid.rho_rh(i)
            mu = grid.mu(i)
            mu_r = grid.mu_r(i)
            dx = grid.dx(i)
            dx_r = grid.dx_r(i)
#            dx_rh = self.center_scheme(grid.dx(i), grid.dx(i+1))
            v = grid.v(i)
            v_r = grid.v_r(i)

            # Calculando as Derivadas dos Resíduos
            #Derivada do resíduo da massa em relação a pressão esquerda
            dRm_dpl = 0

            #Derivada do resíduo da massa em relação a pressão central
            dRm_dpc = 0

            #Derivada do resíduo da massa em relação a pressão direita
            dRm_dpr = 0

            #Derivada do resíduo do momentum em relação a pressão esquerda
            dRp_dpl = 0

            #Derivada do resíduo do momentum em relação a pressão central
            dRp_dpc = -A_rh  # - (2*A_rh*dx_rh)/(dx_r + dx)

            #Derivada do resíduo do momentum em relação a pressão direita
            dRp_dpr = A_rh   # (2*A_rh*dx_rh)/(dx_r + dx)

            #Derivada do resíduo da massa em relação a velocidade face k-1/2
            dRm_dvl = -A_lh*rho_lh

            #Derivada do resíduo da massa em relação a velocidade face k+1/2
            dRm_dvc = A_rh*rho_rh

            #Derivada do resíduo da massa em relação a velocidade face k+3/2
            dRm_dvr = 0

            #Derivada do resíduo do momentum em relação a velocidade face k-1/2
            dvvr_dvlh = 0  # v_r         # CDS Method
            dvv_dvlh = v            # CDS Method
            dRp_dvl = (A_rh*rho_rh)*dvvr_dvlh - (A*rho)*dvv_dvlh - (A*mu)/dx

            #Derivada do resíduo do momentum em relação a velocidade face k+1/2
            dvvr_dvrh = v_r  # 1           # CDS Method
            dvv_dvrh = v  # 1            # CDS Method
            dRp_dvc = (A_r*rho_r)*dvvr_dvrh - (A*rho)*dvv_dvrh + (A_r*mu_r)/dx_r + (A*mu)/dx

            #Derivada do resíduo do momentum em relação a velocidade face k+3/2
            dvvr_dvrrh = 0  #  v_r          # CDS Method
            dvv_dvrrh = v_r  #  v       # CDS Method
            dRp_dvr = (A_r*rho_r)*dvvr_dvrrh - (A*rho)*dvv_dvrrh - (A_r*mu_r)/dx_r

            #Filling Matrix
            if i == 0:
                #Resíduos do momentum
                v_matrix[2*i][(2*i)] = dRp_dpc
                v_matrix[2*i][(2*i)+1] = dRp_dvc
                v_matrix[2*i][(2*i)+2] = dRp_dpr
                v_matrix[2*i][(2*i)+3] = dRp_dvr

                #Resíduos da massa
                v_matrix[2*i+1][(2*i+1)-1] = dRm_dpc
                v_matrix[2*i+1][(2*i+1)] = dRm_dvc
                v_matrix[2*i+1][(2*i+1)+1] = dRm_dpr
                v_matrix[2*i+1][(2*i+1)+2] = dRm_dvr

            elif i == cpoints-1:
                #Resíduos do momentum
                v_matrix[2*i][(2*i)-2] = dRp_dpl    # BC
                v_matrix[2*i][(2*i)-1] = dRp_dvl
                v_matrix[2*i][(2*i)] = dRp_dpc      # BC
                v_matrix[2*i][(2*i)+1] = dRp_dvc

                #Resíduos da massa
                v_matrix[2*i+1][(2*i+1)-3] = dRm_dpl
                v_matrix[2*i+1][(2*i+1)-2] = dRm_dvl
                v_matrix[2*i+1][(2*i+1)-1] = dRm_dpc
                v_matrix[2*i+1][(2*i+1)] = dRm_dvc

            else:
                #Resíduos do momentum
                v_matrix[2*i][(2*i)-2] = dRp_dpl
                v_matrix[2*i][(2*i)-1] = dRp_dvl
                v_matrix[2*i][(2*i)] = dRp_dpc
                v_matrix[2*i][(2*i)+1] = dRp_dvc
                v_matrix[2*i][(2*i)+2] = dRp_dpr
                v_matrix[2*i][(2*i)+3] = dRp_dvr

                #Resíduos da massa
                v_matrix[2*i+1][(2*i+1)-3] = dRm_dpl
                v_matrix[2*i+1][(2*i+1)-2] = dRm_dvl
                v_matrix[2*i+1][(2*i+1)-1] = dRm_dpc
                v_matrix[2*i+1][(2*i+1)] = dRm_dvc
                v_matrix[2*i+1][(2*i+1)+1] = dRm_dpr
                v_matrix[2*i+1][(2*i+1)+2] = dRm_dvr

        return v_matrix

    #Vetor de Resíduos
    def build_Residual_Vector(self, grid):
        #classe que constrói o vetor de resíduosque será resolvido pelo solver
        #Observa-se que o vetor de resíduos possui sinal positivo, devendo
        #ser multiplicado por -1 na hora de resolver o sistema

        cpoints = len(grid.cells)
        R = np.zeros(2*cpoints)

        for i in range(cpoints):
            #Termos Fonte
            A = grid.A(i)
            A_r = grid.A_r(i)
            A_lh = grid.A_lh(i)
            A_rh = grid.A_rh(i)
            rho = grid.rho(i)
            rho_r = grid.rho_r(i)
            rho_lh = grid.rho_lh(i)
            rho_rh = grid.rho_rh(i)
            mu = grid.mu(i)
            mu_r = grid.mu_r(i)
            dx = grid.dx(i)
            dx_r = grid.dx_r(i)
            dx_rh = self.center_scheme(grid.dx(i), grid.dx(i+1))
            v = grid.v(i)
            v_r = grid.v_r(i)
            v_lh = grid.v_lh(i)
            v_rh = grid.v_rh(i)
            v_rrh = grid.v_rh(i+1)
            p = grid.p(i)
            p_r = grid.p_r(i)
            msrc = grid.msrc
            
            # Mass Residual
            R_mass = (A_rh*rho_rh*v_rh) - (A_lh*rho_lh*v_lh) - A*msrc*dx
            
            # Momentum residual
            R_mom = (A_r*rho_r*v_r*v_r) - (A*rho*v*v)
            R_mom += A_rh*(p_r-p)  # 2*(A_rh*dx_rh)*(p_r-p)/(dx_r + dx)
            R_mom += -(A_r*mu_r)*(v_rrh - v_rh)/dx_r
            R_mom += (A*mu)*(v_rh - v_lh)/dx

            R[2*i] = R_mom
            R[2*i+1] = R_mass

        return R

    #Funções de Interpolação
    def center_scheme(self, leftprop, rightprop):
        "Classe interpola utilizando diferencas centrais."
        return (rightprop + leftprop)/2


if __name__ == '__main__':

    from GridFluid import *

    gridteste = GridFluid(0, 150, 1, 20, 0)
    for i in range(4):
        gridteste.add_cell(A=0.1, dx=0.01, rho=1, mu=1, v_r=10, p=100)

    modelCD = Model_COUPLE_CDS()

    A = modelCD.build_Jacobian_matrix(gridteste)
    print A
    print

    R = modelCD.build_Residual_Vector(gridteste)
    print R
    print
    
    x = (np.matrix(A).I*np.matrix(R).T).A1
    p = x[::2]
    v = x[1::2]
    print 'dv:', v
    print 'dp:', p
