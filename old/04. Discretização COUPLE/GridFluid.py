# -*- coding: iso-8859-1 -*-
"""
@author: Wagner Queiroz Barros
@date:   Tue Apr 15 10:36:55 2014
@email:  wagnerqb@gmail.com
@brief:  GridFluid class, used to store cell data in Fluid Flow Problems
"""
from __future__ import division
import numpy as np
from Cell import *


class GridFluid():
    "Classe grid fluido, utilizado na resolu��o das equa��es de Navier-Stokes"

    def __init__(self, lbc_t, lbc, rbc_t, rbc, msrc):
            #Atributos
            self.cells = []             # Boundary (0-pressure, 1-velocity)
            self.lbc_t = lbc_t          # Left Boundary Type
            self.lbc = lbc              # Left Boundary Condition
            self.rbc_t = rbc_t          # Right Bounday Type
            self.rbc = rbc              # Right Boundary Condition
            self.msrc = msrc            # Mass Source Term
            self.n = 0
            self.extrapolar = self.extrapolation   

    def __getitem__(self, index):
        "Sobrecarga do operador [ ]"
        return self.cells[index]

    def add_cell(self, A, dx, rho, mu, v_r, p):
        "Adiciona uma celula no grid."
        self.cells.append(CellFluid(A, dx, rho, mu, v_r, p))
        self.n += 1

    def A(self, index):
        "Area no Centro da c�lula"
        if (index < 0):
            # Sem extrapola��o
            #return self[0].A

            # Com extrapola��o
            return self.extrapolation_l(self[0].A, self[1].A)

        if (index > self.n-1):
            # Sem extrapola��o
            #return self[(len(self.cells)-1)].A

            # Com extrapola��o
            return self.extrapolation_r(self[-2].A, self[-1].A)

        return self[index].A

    def A_r(self, index):
        "Area no centro da c�lula a direita"
        return self.A(index+1)

    def A_l(self, index):
        "Area no centro da c�lula a esquerda"
        return self.A(index-1)

    def A_rh(self, index):
        "Area na face direita da c�lula"
        return (self.A(index) + self.A_r(index))/2

    def A_lh(self, index):
        "Area na face esquerda da c�lula"
        return (self.A_l(index) + self.A(index))/2

    def dx(self, index):
        "Delta_x no centro da celula index."
        if (index < 0):
            return self[0].dx

        if (index > self.n-1):
            return self[self.n-1].dx

        return self[index].dx

    def dx_r(self, index):
        "Delta_x no centro da celula direita."
        return self.dx(index+1)

    def dx_l(self, index):
        "Delta_x no centro da celula esquerda."
        return self.dx(index-1)

    def dx_rh(self, index):
        "Delta_x na face direita da c�lula"
        return (self.dx_r(index) + self.dx(index))/2

    def dx_lh(self, index):
        "Delta_x na face esquerda da c�lula."
        return (self.dx(index) + self.dx_l(index))/2

    def rho(self, index):
        "densidade no Centro da c�lula"
        if (index < 0):
            return self[0].rho

        if (index > self.n-1):
            return self[self.n-1].rho

        return self[index].rho

    def rho_r(self, index):
        "densidade no centro da celula a direita"
        return self.rho(index+1)

    def rho_l(self, index):
        "densidade no centro da celula a esquerda"
        return self.rho(index-1)

    def rho_rh(self, index):
        "densidade na face direita da celula"
        return (self.rho(index) + self.rho_r(index))/2

    def rho_lh(self, index):
        "densidade na face esquerda da celula"
        return (self.rho_l(index) + self.rho(index))/2

    def mu(self, index):
        "viscosidade da celula index."
        if (index < 0):
            return self[0].mu

        if (index > self.n-1):
            return self[self.n-1].mu

        return self[index].mu

    def mu_r(self, index):
        "viscosidade no centro da c�lula a direita"
        return self.mu(index+1)

    def mu_l(self, index):
        "viscosidade no centro da c�lula a esquerda"
        return self.mu(index-1)

    def mu_rh(self, index):
        "viscosidade na face direita da c�lula"
        return (self.mu_r(index) + self.mu(index))/2

    def mu_lh(self, index):
        "viscosidade na face esquerda da c�lula"
        return (self.mu(index) + self.mu_l(index))/2

    def v_rh(self, index):
        "velocidade na Face Direita da c�lula"
        if (index < 0):
            if (self.lbc_t == 0):
                # Sem extrapola��o
                #return self[0].v_rh

                #TIP: N�o foi testado ########################################

                # Com extrapola��o
                A_lh = self.extrapolation_lh(self[0].A, self[1].A)
                return self[0].v_r*self.A_rh(0)/A_lh
            else:
                return self.lbc
        if (index > self.n-1):
            if (self.rbc_t == 0):
                # Sem extrapola��o
                #return self[n-1].v_rh

                # Com extrapola��o
                A_rrh = self.extrapolation_rrh(self[-2].A, self[-1].A)
                ##Area interpolada dentro do tubo e considerada constante fora
                return self[-1].v_rh*self.A_r(self.n-1)/A_rrh
                
                ##Area interpolada dentro e fora do tubo
                # return self[-1].v_rh*self.A_rh(self.n-1)/A_rrh
            else:
                return self.rbc

        return self[index].v_rh

    def set_v_rh(self, index, _v_rh):
        "Atribui a velocidade na face direita da c�lula"
        self[index].v_rh = _v_rh

    def v_lh(self, index):
        "velocidade na Face Esquerda da c�lula"
        return self.v_rh(index-1)

    def v(self, index):
        "velocidade no centro da c�lula (CDS Interpolation Method)"
        return (self.v_lh(index) + self.v_rh(index))/2

    def v_r(self, index):
        "velocidade no centro da c�lula direita (CDS Method)"
        return self.v(index + 1)

    def v_l(self, index):
        "velocidade no centro da c�lula esquerda (CDS Method)"
        return self.v(index - 1)

    def p(self, index):
        "press�o no Centro da c�lula"
        if (index < 0):
            if (self.lbc_t == 0):
                return self.lbc
            else:
                # Sem extrapola��o
                #return self[0].p

                # Com extrapola��o
                return self.extrapolation_l(self[0].p, self[1].p)

        if index > self.n-1:
            if (self.rbc_t == 0):
                return self.rbc
            else:
                # Sem extrapola��o
                #return self[(len(self.cells)-1)].p

                # Com extrapola��o
                return self.extrapolation_r(self[-2].p, self[-1].p)

        return self[index].p

    def p_r(self, index):
        "press�o no centro da c�lula direita."
        return self.p(index + 1)

    def p_l(self, index):
        "press�o no centro da c�lula esquerda."
        return self.p(index - 1)

    def p_rh(self, index):
        "press�o na face direita da c�lula."
        return (self.p_r(index) + self.p(index))/2

    def p_rl(self, index):
        "press�o na face esquerda da c�lula."
        return (self.p(index) + self.p_l(index))/2

    def set_p(self, index, _p):
        "Atribui a press�o no centro da c�lula"
        self[index].p = _p

    def msrc(self, index):
        "Termo Fonte"
        return self.msrc

    def get_all_x(self):
        "Retorna um vetor com a posicao x do centro de todas as celulas,"
        "inclusive as condicoes de contorno."
        cpoints = len(self.cells)
        all_x = np.zeros(cpoints+2)

        all_x[0] = 0
        all_x[1] = self[0].dx/2
        for i in range(2, cpoints+1):
            all_x[i] = all_x[i-1] + self[i-1].dx/2 + self[i-2].dx/2

        all_x[cpoints+1] = all_x[cpoints] + self[cpoints-1].dx/2
        return all_x

    def get_all_p(self):
        """Retorna um vetor com a press�o de todas as celulas, inclusive as
        condicoes de contorno."""
        all_P = np.zeros(self.n+2)

        all_P[0] = self.lbc
        for i in range(self.n):
            all_T[i+1] = self[i].p
        all_T[self.n+1] = self.rbc

        return all_P

    def get_p(self):
        """Retorna um vetor com a press�o de todas as celulas."""
        p = []
        for i in range(self.n):
            p.append(self[i].p)
        return np.array(p)

    def get_v_rh(self):
        """Retorna um vetor com a velocidade de todas as celulas."""
        v_rh = []
        for i in range(self.n):
            v_rh.append(self[i].v_rh)
        return np.array(v_rh)

    def extrapolation(self, x1, y1, x2, y2, x):
        """Fun��o que interpola uma reta que passa pelos pontos (x1,y1) e
        (x2,y2) no ponto x."""
        return y1 + (x-x1)*(y2-y1)/(x2-x1)

    def extrapolation_r(self, y1, y2):
        """Fun��o que extrapola uma propriedade y no bloco a direita."""
        x1 = 0
        x2 = self.dx_lh(self.n-1)
        x = x2 + self.dx_rh(self.n-1)
        return self.extrapolation(x1, y1, x2, y2, x)

    def extrapolation_rrh(self, y1, y2):
        """Fun��o que extrapola uma propriedade y na face direita do vizinho a
        direita."""
        x1 = 0
        x2 = self.dx_lh(self.n-1)
        x = x2 + self.dx_rh(self.n-1) + 0.5*self.dx_r(self.n-1)
        return self.extrapolation(x1, y1, x2, y2, x)

    def extrapolation_l(self, y1, y2):
        """Fun��o que extrapola uma propriedade y no bloco a esquerda."""
        x2 = 0
        x1 = -self.dx_rh(0)
        x = x1 - self.dx_lh(0)
        return self.extrapolation(x1, y1, x2, y2, x)

    def extrapolation_lh(self, y1, y2):
        """Fun��o que extrapola uma propriedade y na face direta do vizinho a
        esquerda."""
        x2 = 0
        x1 = -self.dx_rh(0)
        x = x1 - self.dx(0)*0.5
        return self.extrapolation(x1, y1, x2, y2, x)

if __name__ == '__main__':

    grid = GridFluid(1, 200, 0, 500, 0)
    grid.add_cell(A=0.001, dx=0.1, rho=1000, mu=1, v_r=10, p=100)
    grid.add_cell(A=0.002, dx=0.2, rho=2000, mu=2, v_r=20, p=200)
    grid.add_cell(A=0.003, dx=0.3, rho=3000, mu=3, v_r=30, p=300)

    print "Teste A"
    print grid.A(-1)
    print grid.A(1)
    print grid.A(2)
    print

    print "Teste dx"
    print grid.dx(-1)
    print grid.dx(1)
    print grid.dx(2)
    print

    print "Teste rho"
    print grid.rho(-1)
    print grid.rho(1)
    print grid.rho(2)
    print

    print "Teste mu"
    print grid.mu(-1)
    print grid.mu(1)
    print grid.mu(2)
    print

    print "Teste v_hr"
    print grid.v_rh(-1)
    print grid.v_rh(0)
    print grid.v_rh(1)
    print grid.v_rh(2)
    print grid.v_rh(3)
    print

    print "Teste v_lh"
    print grid.v_lh(-1)
    print grid.v_lh(0)
    print grid.v_lh(1)
    print grid.v_lh(2)
    print grid.v_lh(3)
    print

    print "Teste p"
    print grid.p(-1)
    print grid.p(0)
    print grid.p(1)
    print grid.p(2)
    print grid.p(3)
    print
    
    print grid.extrapolar(1,1,2,3,3)
