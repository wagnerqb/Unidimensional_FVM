# -*- coding: iso-8859-1 -*-
"""
@author: Wagner Queiroz Barros
@date:   Wed Apr 09 14:53:01 2014
@email:  wagnerqb@gmail.com
@brief:  Grid class, used to store cell data
"""
from __future__ import division
import numpy as np
from Interpolation import *
from Cell import *
import matplotlib.pyplot as plt


class Grid():
    "Classe Grid."

    def __init__(self, interpolator, lbc, rbc, Source):
        #Atributos
        self.cells = []
        self.interpolator = interpolator
        self.lbc = lbc
        self.rbc = rbc
        self.source = Source

    def __getitem__(self, index):
        "Sobrecarga do operador [ ]"
        return self.cells[index]

    #Métodos
    def add_cell(self, A, k, T, dx):
        "Adiciona uma celula no grid."
        self.cells.append(Cell(A, k, T, dx))

    def dT_dx_lh(self, index):
        "dT/dx left half"
        if index == 0:
            return 0
        cell = self[index]
        l_cell = self[index-1]
        return self.interpolator.right_derivative(l_cell.T, cell.T, l_cell.dx,
                                                  cell.dx)

    def dT_dx_rh(self, index):
        "dT/dx right half"
        if (index == len(self.cells) - 1):
            return 0
        cell = self[index]
        rightcell = self[index+1]
        return self.interpolator.right_derivative(cell.T, rightcell.T,
                                                  cell.dx, rightcell.dx)

    def k_lh(self, index):
        "Kappa na Face Esquerda."
        if (index == 0):
            return self[0].k

        cell = self[index]
        leftcell = self[index-1]
        return self.interpolator.center_scheme(leftcell.k, cell.k)

    def k_rh(self, index):
        "Kappa na Face Direita."
        if (index == len(self.cells) - 1):
            return self[index].k

        cell = self[index]
        rightcell = self[index+1]
        return self.interpolator.center_scheme(cell.k, rightcell.k)

    def A(self, index):
        "Area no Centro da célula"
        return self[index].A

    def A_lh(self, index):
        "Area na Face Esquerda"
        if (index == 0):
            return self[0].A

        cell = self[index]
        leftcell = self[index-1]
        return self.interpolator.center_scheme(leftcell.A, cell.A)

    def A_rh(self, index):
        "Area na Face Direita"
        if (index == (len(self.cells) - 1)):
            return self[index].A

        cell = self[index]
        rightcell = self[index+1]
        return self.interpolator.center_scheme(cell.A, rightcell.A)

    def dx(self, index):
        "Delta_x da celula index."
        if (index < 0) or (index > (len(self.cells)-1)):
            return 0

        return self[index].dx

    def get_T(self, index):
        "Pega a temperatura da celula index."
        if index == -1:
            return self.lbc

        if index == len(self.cells):
            return self.rbc

        return self.cells[index].T

    def get_all_T(self):
        "Retorna um vetor com a temperatura de todas as celulas, inclusive as"
        "condicoes de contorno."
        cpoints = len(self.cells)
        all_T = np.zeros(cpoints+2)

        all_T[0] = self.lbc
        for i in range(cpoints):
            all_T[i+1] = self[i].T
        all_T[cpoints+1] = self.rbc

        return all_T

    def set_T(self, index, _T):
        "Seta a temperatura da celula index."
        self[index].T = _T

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

    def print_T(self):
        "Imprime na tela uma lista com a temperatura de todas as celulas"
        cpoints = len(self.cells)
        for i in range(cpoints):
            print self[i].T

    def plot_T(self):
        "Plota um grafico com os dados da Temperatura"
        plt.figure()
        x = self.get_all_x()
        y = self.get_all_T()

        plt.plot(x, y, '-d')

        #Configuracoes de grafico
        plt.title(u"Temperatura vs. Posição")
        plt.xlabel(u"X Position [m]")
        plt.ylabel(u"Temperature [ºC]")
        plt.xlim(0, max(x))
        plt.ylim(0, max(y))
        plt.grid()

        plt.show()

    def Source(self, index):
        "Retorna o termo fonte da célula index"
        return self.source


if __name__ == '__main__':

    from Model import *

    ncells = 5

    grid = Grid(Interpolation(), 100, 500, 1000)
    model = Model()

    for i in range(ncells):
        grid.add_cell(10e-3, 1000, 100, 0.1)

    #Resolvendo o Sistema
    A = np.matrix(model.build_matrix(grid))
    b = np.matrix(model.build_coef_vector(grid))
    x = (A.I*b.T).A1

    #Atribuindo os resultados
    for i in range(ncells):
        grid.set_T(i, x[i])

    grid.print_T()

    grid.plot_T()
