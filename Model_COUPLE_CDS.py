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
    """Classe de Modelo para Equa��es de Navier Stokes, utilizando
    o modelo COUPLE para o acoplamento press�o velocidade, e o
    esquema CDS Central Discretization Scheme para interpolar a velocidade."""

    def build_v_matrix(self, grid):
        "Classe que constr�i a matrz de v* (previs�o de velocidade)"
        cpoints = len(grid.cells)
        v_matrix = np.zeros((2*cpoints, 2*cpoints))

        for i in range(cpoints):

            A = grid.A(i)
            A_r = grid.A(i+1)
            A_lh = self.center_scheme(grid.A(i-1), grid.A(i))
            A_rh = self.center_scheme(grid.A(i), grid.A(i+1))
            rho = grid.rho(i)
            rho_r = grid.rho(i+1)
            mu = grid.mu(i)
            mu_r = grid.mu(i+1)
            dx = grid.dx(i)
            dx_r = grid.dx(i+1)
            v = self.center_scheme(grid.v_r(i-1), grid.v_r(i))
            v_r = self.center_scheme(grid.v_r(i), grid.v_r(i+1))

            # Calculando as Derivadas dos Res�duos
            #Derivada do res�duo da massa em rela��o a press�o esquerda            
            dRm_dpl = 0
            
            #Derivada do res�duo da massa em rela��o a press�o central            
            dRm_dpc = 0
            
            #Derivada do res�duo da massa em rela��o a press�o direita            
            dRm_dpr = 0
            
            #Derivada do res�duo do momentum em rela��o a press�o esquerda            
            dRp_dpl = 0
            
            #Derivada do res�duo do momentum em rela��o a press�o central            
            dRp_dpc = 1
            
            #Derivada do res�duo do momentum em rela��o a press�o direita
            dRp_dpr = 1
            
            #Derivada do res�duo da massa em rela��o a velocidade face k-1/2            
            dRm_dvl = 1
            
            #Derivada do res�duo da massa em rela��o a velocidade face k+1/2
            dRm_dvc = 1
            
            #Derivada do res�duo da massa em rela��o a velocidade face k+3/2
            dRm_dvr = 0
            
            #Derivada do res�duo do momentum em rela��o a velocidade face k-1/2            
            dRp_dvl = 1
            
            #Derivada do res�duo do momentum em rela��o a velocidade face k+1/2            
            dRp_dvc = 1
            
            #Derivada do res�duo do momentum em rela��o a velocidade face k+3/2            
            dRp_dvr = 1

            #Filling Matrix
            if i == 0:
                #Res�duos da massa
                v_matrix[2*i][2*i] = dRm_dpc
                v_matrix[2*i][2*i+1] = dRm_dvc
                v_matrix[2*i][2*i+2] = dRm_dpr
                v_matrix[2*i][2*i+3] = dRm_dvr
                
                #Res�duos do momentum
                v_matrix[2*i+1][(2*i+1)-1] = dRp_dpc
                v_matrix[2*i+1][(2*i+1)] = dRp_dvc
                v_matrix[2*i+1][(2*i+1)+1] = dRp_dpr
                v_matrix[2*i+1][(2*i+1)+2] = dRp_dvr
                    
                
            else:
                if i == cpoints-1:
                    #Res�duos da massa
                    v_matrix[2*i][2*i-2] = dRm_dpl
                    v_matrix[2*i][2*i-1] = dRm_dvl
                    v_matrix[2*i][2*i] = dRm_dpc
                    v_matrix[2*i][2*i+1] = dRm_dvc
                    
                    #Res�duos do momentum
                    v_matrix[2*i+1][(2*i+1)-3] = dRp_dpl
                    v_matrix[2*i+1][(2*i+1)-2] = dRp_dvl
                    v_matrix[2*i+1][(2*i+1)-1] = dRp_dpc
                    v_matrix[2*i+1][(2*i+1)] = dRp_dvc
                                        
                else:
                    #Res�duos da massa
                    v_matrix[2*i][2*i-2] = dRm_dpl
                    v_matrix[2*i][2*i-1] = dRm_dvl
                    v_matrix[2*i][2*i] = dRm_dpc
                    v_matrix[2*i][2*i+1] = dRm_dvc
                    v_matrix[2*i][2*i+2] = dRm_dpr
                    v_matrix[2*i][2*i+3] = dRm_dvr
                    
                    #Res�duos do momentum
                    v_matrix[2*i+1][(2*i+1)-3] = dRp_dpl
                    v_matrix[2*i+1][(2*i+1)-2] = dRp_dvl
                    v_matrix[2*i+1][(2*i+1)-1] = dRp_dpc
                    v_matrix[2*i+1][(2*i+1)] = dRp_dvc
                    v_matrix[2*i+1][(2*i+1)+1] = dRp_dpr
                    v_matrix[2*i+1][(2*i+1)+2] = dRp_dvr
                    
        return v_matrix

    #Fun��es de Interpola��o
    def center_scheme(self, leftprop, rightprop):
        "Classe interpola utilizando diferencas centrais."
        return (rightprop + leftprop)/2

if __name__ == '__main__':

    from GridFluid import *

    gridteste = GridFluid(100, 500, 0)
    for i in range(5):
        gridteste.add_cell(10e-3, 1000, 0.1, 100, 1, 2)

    modelCD = Model_COUPLE_CDS()

    A = modelCD.build_v_matrix(gridteste)
    print A
    print

    print "TESTE"
