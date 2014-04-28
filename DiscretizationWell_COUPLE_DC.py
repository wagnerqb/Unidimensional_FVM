# -*- coding: iso-8859-1 -*-
"""
@author: Bismarck Gomes Souza Jr e Wagner Queiroz Barros
@date:   Tue Apr 15 10:55:33 2014
@email:  bismarckgomes@gmail.com e wagnerqb@gmail.com
@brief:  Classe de discretização utilizando o modelo COUPLE Donor-Cell
         para interpolar a derivada da velocidade.
"""

from __future__ import division
import numpy as np


class DiscretizationWell_COUPLE_DC():
    "Modelo de discretização COUPLE Donor-Cell."

    def iterate_x(self, grid, dt, it_max=1E8, erro_max=1E-8):
        "Função que realiza iterações no espaço utilizando o modelo COUPLE."
        it = 0
        erro = erro_max
        new_p = new_v_rh = 0
        print '='*40
        print ' It |  E(v_rh)  |   E(p)    |  E_tot  '
        print '----+-----------+-----------+----------'
        while(it < it_max and erro >= erro_max):
            A = np.matrix(self.build_Jacobian_matrix(grid, dt))
            b = - np.matrix(self.build_Residual_Vector(grid, dt))
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
        
##########################################################################
    def iterate_t(self, grid, dt, it_max=1E8, erro_max=1E-8):
        "Função que realiza iterações no tempo."
        
        pass
##########################################################################
        
    def build_Jacobian_matrix(self, grid, dt):
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
            dx = grid.dx(i)
            dx_rh = grid.dx_rh(i)
            v = grid.v(i)
            v_r = grid.v_r(i)
            v_lh = grid.v_lh(i)
            v_rh = grid.v_rh(i)
            p = grid.p(i)
            p_l = grid.p_l(i)
            p_r = grid.p_r(i)

            # Calculando as Derivadas dos Resíduos
            #Derivada do resíduo da massa em relação a pressão esquerda
            dRm_dpl = self.dRm_dpl(A_lh, v_lh, p_l, fluid)

            #Derivada do resíduo da massa em relação a pressão central
            dRm_dpc = self.dRm_dpc(dt, A, A_lh, A_rh, dx, v_lh, v_rh, p, fluid)

            #Derivada do resíduo da massa em relação a pressão direita
            dRm_dpr = self.dRm_dpr(A_rh, v_rh, p_r, fluid)

            #Derivada do resíduo do momentum em relação a pressão esquerda
            dRp_dpl = self.dRp_dpl()

            #Derivada do resíduo do momentum em relação a pressão central
            dRp_dpc = self.dRp_dpc(dt, A, A_rh, dx_rh, v, v_rh, p, fluid)

            #Derivada do resíduo do momentum em relação a pressão direita
            dRp_dpr = self.dRp_dpr(dt, A_r, A_rh, dx_rh, v_r, v_rh, p_r, fluid)

            #Derivada do resíduo da massa em relação a velocidade face k-1/2
            dRm_dvl = self.dRm_dvl(A_lh, rho_lh)

            #Derivada do resíduo da massa em relação a velocidade face k+1/2
            dRm_dvc = self.dRm_dvc(A_rh, rho_rh)

            #Derivada do resíduo da massa em relação a velocidade face k+3/2
            dRm_dvr = self.dRm_dvr()

            #Derivada do resíduo do momentum em relação a velocidade face k-1/2
            dRp_dvl = self.dRp_dvl(A, rho, v, v_rh)

            #Derivada do resíduo do momentum em relação a velocidade face k+1/2
            dRp_dvc = self.dRp_dvc(dt, A, A_r, A_rh, rho, rho_r, rho_rh, dx_rh,
                                   v, v_r, v_rh)

            #Derivada do resíduo do momentum em relação a velocidade face k+3/2
            dRp_dvr = self.dRp_dvr(A_r, rho_r, v_r, v_rh)

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

    def build_Residual_Vector(self, grid, dt):
        """Função que constrói o vetor de resíduos. Observa-se que o vetor de
        resíduos não possui sinal negativo."""

        n = grid.n
        R = np.zeros(2*n)

        for i in range(n):
            #Variáveis no passo de tempo atual
            A = grid.A(i)
            A_r = grid.A_r(i)
            A_lh = grid.A_lh(i)
            A_rh = grid.A_rh(i)
            rho = grid.rho(i)
            rho_r = grid.rho_r(i)
            rho_lh = grid.rho_lh(i)
            rho_rh = grid.rho_rh(i)
            dx = grid.dx(i)
            dx_rh = grid.dx_rh(i)
            v = grid.v(i)
            v_r = grid.v_r(i)
            v_lh = grid.v_lh(i)
            v_rh = grid.v_rh(i)
            p = grid.p(i)
            p_r = grid.p_r(i)
            msrc = grid.msrc(i)
            fsrc = 0

            #Variáveis no passo de tempo anterior
            v_rh_old = grid.v_rh_old(i)
            rho_old = grid.rho_old(i)
            rho_rh_old = grid.rho_rh_old(i)

            # Mass Residual
            R_mass = A*dx*(rho - rho_old)/dt
            R_mass += (A_rh*rho_rh*v_rh) - (A_lh*rho_lh*v_lh) - A*msrc*dx

            # Momentum residual
            R_mom = A_rh*dx_rh*(rho_rh*v_rh - rho_rh_old*v_rh_old)/dt
            R_mom += A_r*rho_r*v_r*v_r - A*rho*v*v + A_rh*(p_r-p) - A*dx*fsrc

            R[2*i] = R_mom
            R[2*i+1] = R_mass

        return R

    #======================= DERIVADAS DO RESÍDUO DA MASSA ===================#
    def dRm_dpl(self, A_lh, v_lh, p_l, fluid):
        "Derivada do resíduo da massa em relação a pressão esquerda."
        return -0.5*A_lh*v_lh*fluid.drho_dp(p_l)

    def dRm_dpc(self, dt, A, A_lh, A_rh, dx, v_lh, v_rh, p, fluid):
        "Derivada do resíduo da massa em relação a pressão central."
        return (A*dx/dt + 0.5*A_rh*v_rh - 0.5*A_lh*v_lh)*fluid.drho_dp(p)

    def dRm_dpr(self, A_rh, v_rh, p_r, fluid):
        "Derivada do resíduo da massa em relação a pressão direita."
        return 0.5*A_rh*v_rh*fluid.drho_dp(p_r)

    def dRm_dvl(self, A_lh, rho_lh):
        "Derivada do resíduo da massa em relação a velocidade esquerda."
        return -A_lh*rho_lh

    def dRm_dvc(self, A_rh, rho_rh):
        "Derivada do resíduo da massa em relação a velocidade central."
        return A_rh*rho_rh

    def dRm_dvr(self):
        "Derivada do resíduo da massa em relação a velocidade direita."
        return 0

    #============= DERIVADAS DO RESÍDUO DA QUANTIDADE DE MOVIMENTO ===========#
    def dRp_dpl(self):
        "Derivada do resíduo do momentum em relação a pressão esquerda"
        return 0

    def dRp_dpc(self, dt, A, A_rh, dx_rh, v, v_rh, p, fluid):
        "Derivada do resíduo do momentum em relação a pressão central"
        return ((A_rh*dx_rh*v_rh/(2*dt) - (A*v*v))*fluid.drho_dp(p) - A_rh)

    def dRp_dpr(self, dt, A_r, A_rh, dx_rh, v_r, v_rh, p_r, fluid):
        "Derivada do resíduo do momentum em relação a pressão direita"
        return ((A_rh*dx_rh*v_rh/(2*dt) + (A_r*v_r*v_r))*fluid.drho_dp(p_r)
                + A_rh)

    def dRp_dvl(self, A, rho, v, v_rh):
        "Derivada do resíduo do momentum em relação a velocidade face k-1/2"
        if (v_rh >= 0):
            return (-(A*rho)*2*v)
        elif (v_rh < 0):
            return 0

    def dRp_dvc(self, dt, A, A_r, A_rh, rho, rho_r, rho_rh, dx_rh, v, v_r,
                v_rh):
        "Derivada do resíduo do momentum em relação a velocidade face k+1/2"
        if (v_rh >= 0):
            return ((A_rh*dx_rh*rho_rh)/dt + (A_r*rho_r*2*v_r))
        elif (v_rh < 0):
            return ((A_rh*dx_rh*rho_rh)/dt - (A*rho*2*v))

    def dRp_dvr(self, A_r, rho_r, v_r, v_rh):
        "Derivada do resíduo do momentum em relação a velocidade face k+3/2"
        if (v_rh >= 0):
            return 0
        elif (v_rh < 0):
            return (A_r*rho_r*2*v_r)

    #=========================================================================#
    #####################   FUNÇÕES DE TESTE
    #========================= RESÍDUO DA MASSA ==============================#
    def mass_residual(self, index, v_lh, v_rh, dt):
        n = grid.n
        R = np.zeros(2*n)

        #Variáveis no passo de tempo atual
        A = grid.A(index)
        A_lh = grid.A_lh(index)
        A_rh = grid.A_rh(index)
        rho = grid.rho(index)
        rho_lh = grid.rho_lh(index)
        rho_rh = grid.rho_rh(index)
        dx_rh = grid.dx_rh(index)
#        v_lh = grid.v_lh(index)
#        v_rh = grid.v_rh(index)
        msrc = grid.msrc(index)

        #Variáveis no passo de tempo anterior
        rho_old = grid.rho_old(index)

        # Mass Residual
        R_mass = A*dx*(rho - rho_old)/dt
        R_mass += (A_rh*rho_rh*v_rh) - (A_lh*rho_lh*v_lh) - A*msrc*dx
        
        return R_mass

    def dmr_dv_rh(self, index, dt, eps):
        
        v_lh = grid.v_lh(index)
        v_rh = grid.v_rh(index)
        
        fpos = self.mass_residual(index, v_lh, (v_rh + eps), dt)
        fat = self.mass_residual(index, v_lh, (v_rh), dt)
        
        return (fpos - fat)/eps

    def dmr_dv_lh(self, index, dt, eps):
        
        v_lh = grid.v_lh(index)
        v_rh = grid.v_rh(index)
        
        fpos = self.mass_residual(index, (v_lh + eps), v_rh, dt)
        fat = self.mass_residual(index, (v_lh), v_rh, dt)
        
        return (fpos - fat)/eps

    def momentum_residual(self, index, v_rh, v_lh, p, p_r, dt):
#            #Variáveis no passo de tempo atual
            A = grid.A(index)
            A_r = grid.A_r(index)
            A_rh = grid.A_rh(index)
            rho = grid.rho(index)
            rho_r = grid.rho_r(index)
            rho_rh = grid.rho_rh(index)
            dx = grid.dx(index)
            dx_rh = grid.dx_rh(index)
            v = v_lh
            v_r = v_rh
            fsrc = 0
#
#            #Variáveis no passo de tempo anterior
            v_rh_old = grid.v_rh_old(index)
            rho_rh_old = grid.rho_rh_old(index)
    
            # Momentum residual
            R_mom = A_rh*dx_rh*(rho_rh*v_rh - rho_rh_old*v_rh_old)/dt
            R_mom += A_r*rho_r*v_r*v_r - A*rho*v*v + A_rh*(p_r-p) - A*dx*fsrc

            return R_mom
        
    def dmp_dp(self, index, dt, eps):
            
            v_rh = grid.v_rh(index)
            v_lh = grid.v_lh(index)
            p = grid.p(index)
            p_r = grid.p_r(index)
            
            fpos = self.momentum_residual(index, v_rh, v_lh, (p + eps), p_r, dt)
            fat = self.momentum_residual(index, v_rh, v_lh, (p), p_r, dt)
 
            return (fpos - fat)/eps
    
    def dmp_dv_rh(self, index, dt, eps):
            
            v_rh = grid.v_rh(index)
            v_lh = grid.v_lh(index)
            p = grid.p(index)
            p_r = grid.p_r(index)
            
            fpos = self.momentum_residual(index, (v_rh + eps), v_lh, p, p_r, dt)
            fat = self.momentum_residual(index, (v_rh), v_lh, p, p_r, dt)
 
            return (fpos - fat)/eps
            
    def dmp_dv_lh(self, index, dt, eps):
            
            v_rh = grid.v_rh(index)
            v_lh = grid.v_lh(index)
            p = grid.p(index)
            p_r = grid.p_r(index)
            
            fpos = self.momentum_residual(index, v_rh, (v_lh + eps), p, p_r, dt)
            fat = self.momentum_residual(index, v_rh, v_lh, p, p_r, dt)
 
            return (fpos - fat)/eps
        
if __name__ == '__main__':

    from GridWell import *
    from Fluid import *
    from Cell import *

    # Pipe properties
    A = 0.1
    dx = 100
    ncells = 5

    # Numerical Parameters
    dt = 1

    # Fluid properties
    rho = 1                       # Fluid density
    msrc = 0.                       # Mass Source term per volume unity
    fsrc = 0.                       # Termo fonte da QM

    # Initial properties
    v_ini = 0                       # Initial Condition for v
    p_ini = 5                       # Initial Condition for p

    # Boundary Condition
    # Left Boundary Condtion (Velocidade na entrada: v_ini)
    lbc_t = 1                       # LBC Type (0 - Pressure / 1 - Velocity)
    lbc = 5                         # LBC Value

    # Right Boundary Condtion (Pressão na saida: p_ini)
    rbc_t = 0                       # RBC Type (0 - Pressure / 1 - Velocity)
    rbc = 10                         # RBC Value

    #Fluido
    fluid = FluidIncompressible(rho)

    # Creating Grid
    grid = GridWell(fluid)
    grid.set_boundaries(lbc_t, lbc, rbc_t, rbc)

    for i in range(ncells):
        #Criando grid
        cell = CellWell(A, dx, fluid, v_ini, p_ini, msrc)
        grid.add_cell(cell)

    # Creating Model
    model = DiscretizationWell_COUPLE_DC()
    
    for tst in range(2):
        print "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
        J = model.build_Jacobian_matrix(grid, dt)
        J = np.matrix(J)
        r = - model.build_Residual_Vector(grid, dt)
        r = np.matrix(r)
        
        x = (J.I*r.T).A1
    
        print "JACOBIAN TEST"
        print J
        
        print "JACOBIAN Inverse"
        print J.I
    
        print "RESIDUAL TEST"
        print r.T
        
    #    print "RESULTS"
    #    print x
        
         #Valores novos
        new_p = x[::2]
        new_v_rh = x[1::2]
    
    #    model.iterate_x(grid, dt, 100, 1e-5)
    #    
        print 'Delta P', new_p
        print 'Delat v', new_v_rh
        
        #Atualizando valores
        for i in range(grid.n):
            grid.set_p(i, grid.p(i)+new_p[i])
            grid.set_v_rh(i, grid.v_rh(i)+new_v_rh[i])
         
         
        print 'NOVO P', grid.get_p()
        print 'Novo v', grid.get_v_rh()

    
#    index = 1
#    eps = 0.01
#    v_lh = grid.v_lh(index)
#    v_rh = grid.v_rh(index)
#    p = grid.p(index)
#    p_r = grid.p_r(index)
#    
##    print v_lh, v_rh
##    mass_residual(index, v_lh, v_rh, dt): 
#    print "MASS"
#    print "mass residual", model.mass_residual(index, v_lh, v_rh, dt)
#    print "dMr_dv_rh", model.dmr_dv_rh(index, dt, eps)
#    print "dMr_dv_lh", model.dmr_dv_lh(index, dt, eps)
#    
#    print "MOMENTUM"
#    print "Momentum Residual", model.momentum_residual(index, v_rh, v_lh, p, p_r, dt)
#    print "dMp_dp", model.dmp_dp(index, dt, eps)
#    print "dMp_dv_rh", model.dmp_dv_rh(index, dt, eps)
#    print "dMp_dv_lh", model.dmp_dv_lh(index, dt, eps)
    
