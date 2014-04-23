# -*- coding: iso-8859-1 -*-
"""
@author: Wagner Queiroz Barros
@date:   Tue Apr 15 10:36:55 2014
@email:  wagnerqb@gmail.com
@brief:  GridFluid_SIMPLE class, used to store cell data in Fluid Flow Problems
"""
from __future__ import division
from Cell import *


class GridFluid_SIMPLE():
    "Classe grid fluido para o método SIMPLE, utilizado na resolução das equações de Navier-Stokes"

    def __init__(self, lbc_t, lbc, rbc_t, rbc, msrc):
            #Atributos
            self.cells = []             # Boundary (0-pressure, 1-velocity)
            self.lbc_t = lbc_t          # Left Boundary Type
            self.lbc = lbc              # Left Boundary Condition
            self.rbc_t = rbc_t          # Right Bounday Type
            self.rbc = rbc              # Right Boundary Condition
            self.msrc = msrc            # Mass Source Term

    def __getitem__(self, index):
        "Sobrecarga do operador [ ]"
        return self.cells[index]

    #Métodos
    def add_cell(self, A, dx, rho, mu, v_r, p):
        "Adiciona uma celula no grid."
        self.cells.append(CellFluid(A, dx, rho, mu, v_r, p))

    def A(self, index):
        "Area no Centro da célula"
        if (index < 0):
            return self[0].A

        if (index > (len(self.cells)-1)):
            return self[(len(self.cells)-1)].A

        return self[index].A

    def A_r(self, index):
        "Area no centro da célula a direita"
        return self.A(index+1)

    def A_l(self, index):
        "Area no centro da célula a esquerda"
        return self.A(index-1)

    def A_rh(self, index):
        "Area na face direita da célula"
        return (self.A(index) + self.A_r(index))/2

    def A_lh(self, index):
        "Area na face esquerda da célula"
        return (self.A_l(index) + self.A(index))/2

    def dx(self, index):
        "Delta_x no centro da celula index."
        if (index < 0):
            return self[0].dx

        if (index > (len(self.cells)-1)):
            return self[(len(self.cells)-1)].dx

        return self[index].dx

    def dx_r(self, index):
        "Delta_x no centro da celula direita."
        return self.dx(index+1)

    def dx_l(self, index):
        "Delta_x no centro da celula esquerda."
        return self.dx(index-1)

    def dx_rh(self, index):
        "Delta_x na face direita da célula"
        return (self.dx_r(index) + self.dx(index))/2

    def dx_lh(self, index):
        "Delta_x na face esquerda da célula."
        return (self.dx(index) + self.dx_l(index))/2

    def rho(self, index):
        "densidade no Centro da célula"
        if (index < 0):
            return self[0].rho

        if (index > (len(self.cells)-1)):
            return self[(len(self.cells)-1)].rho

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

        if (index > (len(self.cells)-1)):
            return self[(len(self.cells)-1)].mu

        return self[index].mu

    def mu_r(self, index):
        "viscosidade no centro da célula a direita"
        return self.mu(index+1)

    def mu_l(self, index):
        "viscosidade no centro da célula a esquerda"
        return self.mu(index-1)

    def mu_rh(self, index):
        "viscosidade na face direita da célula"
        return (self.mu_r(index) + self.mu(index))/2

    def mu_lh(self, index):
        "viscosidade na face esquerda da célula"
        return (self.mu(index) + self.mu_l(index))/2

    def v_rh(self, index):
        "velocidade na Face Direita da célula"
        if (index < 0):
            if (self.lbc_t == 0):
                return self[0].v_rh
            else:
                return self.lbc

        if (index > (len(self.cells)-1)):
            if (self.rbc_t == 0):
                return self[(len(self.cells)-1)].v_rh
            else:
                return self.rbc

        return self[index].v_rh

    def set_v_rh(self, index, _v_rh):
        "Atribui a velocidade na face direita da célula"
        self[index].v_rh = _v_rh

    def v_lh(self, index):
        "velocidade na Face Esquerda da célula"
        return self.v_rh(index-1)

    def v(self, index):
        "velocidade no centro da célula (CDS Interpolation Method)"
        return (self.v_lh(index) + self.v_rh(index))/2

    def v_r(self, index):
        "velocidade no centro da célula direita (CDS Method)"
        return self.v(index + 1)

    def v_l(self, index):
        "velocidade no centro da célula esquerda (CDS Method)"
        return self.v(index - 1)

    def p(self, index):
        "pressão no Centro da célula"
        if (index < 0):
            if (self.lbc_t == 0):
                return self.lbc
            else:
                return self[0].p

        if (index > (len(self.cells)-1)):
            if (self.rbc_t == 0):
                return self.rbc
            else:
                return self[(len(self.cells)-1)].p

        return self[index].p

    def p_r(self, index):
        "pressão no centro da célula direita."
        return self.p(index + 1)

    def p_l(self, index):
        "pressão no centro da célula esquerda."
        return self.p(index - 1)

    def p_rh(self, index):
        "pressâo na face direita da célula."
        return (self.p_r(index) + self.p(index))/2

    def p_rl(self, index):
        "pressâo na face esquerda da célula."
        return (self.p(index) + self.p_l(index))/2

    def set_p(self, index, _p):
        "Atribui a pressão no centro da célula"
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
        "Retorna um vetor com a pressão de todas as celulas, inclusive as"
        "condicoes de contorno."
        cpoints = len(self.cells)
        all_P = np.zeros(cpoints+2)

        all_P[0] = self.lbc
        for i in range(cpoints):
            all_T[i+1] = self[i].p
        all_T[cpoints+1] = self.rbc

        return all_P

if __name__ == '__main__':

    grid = GridFluid_SIMPLE(1, 200, 0, 500, 0)
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
