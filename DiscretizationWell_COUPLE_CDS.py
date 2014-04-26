# -*- coding: iso-8859-1 -*-
"""
@author: Bismarck Gomes Souza Jr e Wagner Queiroz Barros
@date:   Tue Apr 15 10:55:33 2014
@email:  bismarckgomes@gmail.com e wagnerqb@gmail.com
@brief:  Classe de discretização utilizando o modelo COUPLE CDS para interpolar
        a derivada da velocidade.
"""

from __future__ import division
import numpy as np


class DiscretizationWell_COUPLE_CDS():
    "Modelo de discretização COUPLE CDS."

    def iterate_x(self, grid, it_max=1E8, erro_max=1E-8):
        "Função que realiza iterações utilizando o modelo COUPLE."
        it = 0
        erro = erro_max
        new_p = new_v_rh = 0
        print '='*40
        print ' It |  E(v_rh)  |   E(p)    |  E_tot  '
        print '----+-----------+-----------+----------'
        while(it < it_max and erro >= erro_max):
            A = np.matrix(self.build_Jacobian_matrix(grid))
            b = np.matrix(self.build_Residual_Vector(grid))
            x = (A.I*b.T).A1

            #Valores novos
            new_p = x[::2]
            new_v_rh = x[1::2]

            #Calculando norma L2
            erro_p = (sum(new_p**2))**0.5
            erro_v_rh = (sum(new_v_rh**2))**0.5

            #Verificando alteração de informações
            if erro == erro_p + erro_v_rh:
                print '\n!!! Erro mínimo obtido: %.3e !!!\n' % erro
                raw_input('Digite enter para continuar...')
                break
            else:
                erro = erro_p + erro_v_rh

            #Atualizando valores
            it += 1
            for i in range(grid.n):
                grid.set_p(i, grid.p(i)+new_p[i])
                grid.set_v_rh(i, grid.v_rh(i)+new_v_rh[i])

            #Imprimindo resultados da iteração
            print '{0:03d} | {1:.3E} | {2:.3E} | {3:.3E}'.format(it,
                erro_v_rh, erro_p, erro_v_rh+erro_p)
        print '='*40
        return it, erro

    def build_Jacobian_matrix(self, grid):
        "Função que constrói o Jacobiano"
        n = grid.n
        J = np.zeros((2*n, 2*n))

        for i in range(n):
            A = grid.A(i)
            A_r = grid.A_r(i)
            A_lh = grid.A_lh(i)
            A_rh = grid.A_rh(i)
            rho = grid.rho(i)
            rho_r = grid.rho_r(i)
            rho_lh = grid.rho_lh(i)
            rho_rh = grid.rho_rh(i)
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
            dvvr_dvlh = 0                # CDS Method
            dvv_dvlh = v                 # CDS Method
            dRp_dvl = (A_rh*rho_rh)*dvvr_dvlh - (A*rho)*dvv_dvlh

            #Derivada do resíduo do momentum em relação a velocidade face k+1/2
            dvvr_dvrh = v_r              # CDS Method
            dvv_dvrh = v                 # CDS Method
            dRp_dvc = (A_r*rho_r)*dvvr_dvrh - (A*rho)*dvv_dvrh

            #Derivada do resíduo do momentum em relação a velocidade face k+3/2
            dvvr_dvrrh = 0               # CDS Method
            dvv_dvrrh = v_r              # CDS Method
            dRp_dvr = (A_r*rho_r)*dvvr_dvrrh - (A*rho)*dvv_dvrrh

            #Filling Matrix
            if i == 0:
                #Resíduos do momentum
                J[2*i][(2*i)] = dRp_dpc
                J[2*i][(2*i)+1] = dRp_dvc
                J[2*i][(2*i)+2] = dRp_dpr
                J[2*i][(2*i)+3] = dRp_dvr

                #Resíduos da massa
                J[2*i+1][(2*i+1)-1] = dRm_dpc
                J[2*i+1][(2*i+1)] = dRm_dvc
                J[2*i+1][(2*i+1)+1] = dRm_dpr
                J[2*i+1][(2*i+1)+2] = dRm_dvr

            elif i == n-1:
                #Resíduos do momentum
                J[2*i][(2*i)-2] = dRp_dpl
                J[2*i][(2*i)-1] = dRp_dvl
                J[2*i][(2*i)] = dRp_dpc
                J[2*i][(2*i)+1] = dRp_dvc

                #Resíduos da massa
                J[2*i+1][(2*i+1)-3] = dRm_dpl
                J[2*i+1][(2*i+1)-2] = dRm_dvl
                J[2*i+1][(2*i+1)-1] = dRm_dpc
                J[2*i+1][(2*i+1)] = dRm_dvc

            else:
                #Resíduos do momentum
                J[2*i][(2*i)-2] = dRp_dpl
                J[2*i][(2*i)-1] = dRp_dvl
                J[2*i][(2*i)] = dRp_dpc
                J[2*i][(2*i)+1] = dRp_dvc
                J[2*i][(2*i)+2] = dRp_dpr
                J[2*i][(2*i)+3] = dRp_dvr

                #Resíduos da massa
                J[2*i+1][(2*i+1)-3] = dRm_dpl
                J[2*i+1][(2*i+1)-2] = dRm_dvl
                J[2*i+1][(2*i+1)-1] = dRm_dpc
                J[2*i+1][(2*i+1)] = dRm_dvc
                J[2*i+1][(2*i+1)+1] = dRm_dpr
                J[2*i+1][(2*i+1)+2] = dRm_dvr

        return J

    def build_Residual_Vector(self, grid):
        """Função que constrói o vetor de resíduos. Observa-se que o vetor de
        resíduos já possui sinal negativo."""

        n = grid.n
        R = np.zeros(2*n)

        for i in range(n):
            #Termos Fonte
            A = grid.A(i)
            A_r = grid.A_r(i)
            A_lh = grid.A_lh(i)
            A_rh = grid.A_rh(i)
            rho = grid.rho(i)
            rho_r = grid.rho_r(i)
            rho_lh = grid.rho_lh(i)
            rho_rh = grid.rho_rh(i)
            dx = grid.dx(i)
            v = grid.v(i)
            v_r = grid.v_r(i)
            v_lh = grid.v_lh(i)
            v_rh = grid.v_rh(i)
            p = grid.p(i)
            p_r = grid.p_r(i)
            msrc = grid.msrc(i)

            # Mass Residual
            R_mass = (A_rh*rho_rh*v_rh) - (A_lh*rho_lh*v_lh) - A*msrc*dx

            # Momentum residual
            R_mom = (A_r*rho_r*v_r*v_r) - (A*rho*v*v)
            R_mom += A_rh*(p_r-p)  # 2*(A_rh*dx_rh)*(p_r-p)/(dx_r + dx)

            R[2*i] = R_mom
            R[2*i+1] = R_mass

        return -R


if __name__ == '__main__':

    from GridFluid import *

    gridteste = GridFluid(1, 15, 0, 20, 0)
    for i in range(4):
        gridteste.add_cell(A=0.3, dx=0.1, rho=1, mu=1, v_r=1, p=10)

    modelCD = Model_COUPLE_CDS()

    A = modelCD.build_Jacobian_matrix(gridteste)
    print A
    print

    R = modelCD.build_Residual_Vector(gridteste)
    print R
    print
    
    x = (np.matrix(A).I*np.matrix(R).T).A1
    dp = x[::2]
    dv = x[1::2]
    print 'dv:', dv
    print 'dp:', dp
