# -*- coding: iso-8859-1 -*-
"""
@author: Wagner Queiroz Barros
@date:   Tue Apr 15 10:36:55 2014
@email:  wagnerqb@gmail.com
@brief:  GridFluid class, used to store cell data in Fluid Flow Problems
"""
from __future__ import division
from Cell import *


class GridFluid():
    "Classe grid fluido, utilizado na resolução das equações de Navier-Stokes"

    def __init__(self, lbc, rbc, msrc):
            #Atributos
            self.cells = []
            self.lbc = lbc              # Left Boundary Condition
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


    def dx(self, index):
        "Delta_x da celula index."
        if (index < 0):
            return self[0].dx

        if (index > (len(self.cells)-1)):
            return self[(len(self.cells)-1)].dx


        return self[index].dx

    def rho(self, index):
        "densidade no Centro da célula"
        if (index < 0):
            return self[0].rho

        if (index > (len(self.cells)-1)):
            return self[(len(self.cells)-1)].rho

        return self[index].rho

    def mu(self, index):
        "Delta_x da celula index."
        if (index < 0):
            return self[0].mu

        if (index > (len(self.cells)-1)):
            return self[(len(self.cells)-1)].mu

        return self[index].mu


    def v_r(self, index):
        "velocidade na Face Direita da célula"
        if (index < 0):
            return self[0].v_r

        if (index > (len(self.cells)-1)):
            return self[(len(self.cells)-1)].v_r

        return self[index].v_r


    def v_l(self, index):
        "velocidade na Face Esquerda da célula"
        if (index <= 0):
            return self[0].v_r

        if (index > (len(self.cells)-1)):
            return self[(len(self.cells)-1)].v_r

        return self[index-1].v_r

    def p(self, index):
        "pressão no Centro da célula"
        if (index < 0):
            return self[0].p

        if (index > (len(self.cells)-1)):
            return self[(len(self.cells)-1)].p

        return self[index].p


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


    grid = GridFluid(100, 500, 0)
    grid.add_cell(0.001, 0.1, 1000, 1, 10, 100)
    grid.add_cell(0.002, 0.2, 2000, 2, 20, 200)
    grid.add_cell(0.003, 0.3, 3000, 3, 30, 300)

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
    
    print "Teste v_r"
    print grid.v_r(-1)
    print grid.v_r(0)
    print grid.v_r(1)
    print grid.v_r(2)
    print grid.v_r(3)
    print
    
    print "Teste v_l"
    print grid.v_l(-1)
    print grid.v_l(0)
    print grid.v_l(1)
    print grid.v_l(2)
    print grid.v_l(3)
    print
    
    print "Teste p"
    print grid.p(-1)
    print grid.p(1)
    print grid.p(2)
    print
    


