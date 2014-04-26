# -*- coding: iso-8859-1 -*-
"""
@author: Wagner Queiroz Barros
@date:   Fri Apr 11 16:29:03 2014
@email:  wagnerqb@gmail.com
@brief:  Class used to discretize Convective/Difuse PDEs using UDS Scheme
"""


from __future__ import division

from Model_D_CDS import *

class Model_CD_UDS(Model_D_CDS):
    """Classe de Modelo para Equações Difusivas/Convectivas, utilizando
    o esquema UDS Upwind Discretization Scheme ."""
    
    def build_matrix(self, grid):
        "Classe que constrói a matrz que será resolvida pelo solver"
        cpoints = len(grid.cells)
        A = np.zeros((cpoints, cpoints))

        for i in range(cpoints):
            k_lh = self.center_scheme(grid.k(i-1), grid.k(i))
            k_rh = self.center_scheme(grid.k(i), grid.k(i+1))
            A_lh = self.center_scheme(grid.A(i-1), grid.A(i))
            A_rh = self.center_scheme(grid.A(i), grid.A(i+1))
            v_lh = self.center_scheme(grid.v(i-1), grid.v(i))
            v_rh = self.center_scheme(grid.v(i), grid.v(i+1))
            rho_lh = self.center_scheme(grid.rho(i-1), grid.rho(i))
            rho_rh = self.center_scheme(grid.rho(i), grid.rho(i+1))
            
            ##  Velocidades Positivas  ##
            if (grid.v(i) >= 0):
                lc_k = - A_lh*rho_lh*v_lh
                ld_k = - 2*k_lh*A_lh/(grid.dx(i) + grid.dx(i-1))
                L_k = lc_k + ld_k
                
                rc_k = 0.
                rd_k = - 2*k_rh*A_rh/(grid.dx(i+1) + grid.dx(i))
                R_k = rc_k + rd_k
                
                C_k = A_rh*rho_rh*v_rh - rd_k - ld_k
                                
            ##  Velocidades Negativas  ##    
            if (grid.v(i) < 0):
                lc_k = 0.
                ld_k = - 2*k_lh*A_lh/(grid.dx(i) + grid.dx(i-1))
                L_k = lc_k + ld_k
                
                rc_k = A_rh*rho_rh*v_rh
                rd_k = - 2*k_rh*A_rh/(grid.dx(i+1) + grid.dx(i))
                R_k = rc_k + rd_k
                
                C_k = - A_lh*rho_lh*v_lh - rd_k - ld_k
            

            if i == 0:
                A[i][i] = - L_k - R_k
                A[i][i+1] = R_k
            else:
                if i == cpoints-1:
                    A[i][i-1] = L_k
                    A[i][i] = - L_k - R_k
                else:
                    A[i][i-1] = L_k
                    A[i][i] = - L_k - R_k
                    A[i][i+1] = R_k
        return A
        
    def build_coef_vector(self, grid):
        #classe que constrói o vetor que será resolvido pelo solver
        cpoints = len(grid.cells)
        B = np.zeros(cpoints)

        for i in range(cpoints):
            #Termos Fonte
            B[i] = grid.dx(i)*grid.A(i)*grid.Source(i)

        #Condição de Contorno Esquerda
        k_lh_lbc = self.center_scheme(grid.k(-1), grid.k(0))
        A_lh_lbc = self.center_scheme(grid.A(-1), grid.A(0))
        rho_lh_lbc = self.center_scheme(grid.rho(-1), grid.rho(0))
        v_lh_lbc = self.center_scheme(grid.v(-1), grid.v(0))
        
        ##  Velocidades Positivas  ##
        if (grid.v(0) >= 0):
            lc_k = - A_lh_lbc*rho_lh_lbc*v_lh_lbc
            ld_k = - 2*k_lh_lbc*A_lh_lbc/(grid.dx(0) + grid.dx(-1))
            L_k = lc_k + ld_k
        
        ##  Velocidades Negativas  ##    
        if (grid.v(0) < 0):
            lc_k = 0.
            ld_k = - 2*k_lh_lbc*A_lh_lbc/(grid.dx(i) + grid.dx(i-1))
            L_k = lc_k + ld_k
            
        B[0] = B[0] - L_k*grid.phi(-1)

        #Condição de Contorno Direita
        k_rh_rbc = self.center_scheme(grid.k(cpoints - 1), grid.k(cpoints))
        A_rh_rbc = self.center_scheme(grid.A(cpoints - 1), grid.A(cpoints))
        rho_rh_rbc = self.center_scheme(grid.rho(cpoints - 1), grid.rho(cpoints))
        v_rh_rbc = self.center_scheme(grid.v(cpoints - 1), grid.v(cpoints))
        
        ##  Velocidades Positivas  ##
        if (grid.v(0) >= 0):
            rc_k = 0.
            rd_k = - 2*k_rh_rbc*A_rh_rbc/(grid.dx(cpoints) + grid.dx(cpoints - 1))
            R_k = rc_k + rd_k
            
        ##  Velocidades Negativas  ##    
        if (grid.v(0) < 0):
            rc_k = A_rh_rbc*rho_rh_rbc*v_rh_rbc
            rd_k = - 2*k_rh_rbc*A_rh_rbc/(grid.dx(cpoints) + grid.dx(cpoints - 1))
            R_k = rc_k + rd_k

        B[cpoints - 1] = B[cpoints - 1] - R_k*grid.phi(cpoints)

        return B
        
if __name__  == '__main__':

    from GridCD import *

    gridteste = GridCD(100, 500, 0)
    for i in range(5):
        gridteste.add_cell(10e-3, 1000, 0.1, 100, 1, 0)


    modelCD = Model_CD_UDS()

    A = modelCD.build_matrix(gridteste)
    print A
    print
    
    B = modelCD.build_coef_vector(gridteste)
    print B    
