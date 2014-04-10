# -*- coding: iso-8859-1 -*-
"""
@author: Wagner Queiroz Barros
@date:   Wed Apr 09 14:53:01 2014
@email:  wagnerqb@gmail.com
@brief:  Grid class, used to store cell data
"""
from __future__ import division
import numpy as np
#from Interpolation import *
from Cell import *
import matplotlib.pyplot as plt


class Grid():
    "Classe Grid."

    def __init__(self, lbc, rbc, Source):
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
    def add_cell(self, A, k, phi, dx):
        "Adiciona uma celula no grid."
        self.cells.append(Cell(A, k, phi, dx))

    def A(self, index):
        "Area no Centro da célula"
        return self[index].A

    def k(self, index):
        "Delta_x da celula index."
        if (index < 0):
            return self[0].k

        if (index > (len(self.cells)-1)):
            return self[(len(self.cells)-1)].k

        return self[index].k

    def dx(self, index):
        "Delta_x da celula index."
        if (index < 0) or (index > (len(self.cells)-1)):
            return 0

        return self[index].dx

    def get_phi(self, index):
        "Pega a temperatura da celula index."
        if index == -1:
            return self.lbc

        if index == len(self.cells):
            return self.rbc

        return self.cells[index].phi

    def get_all_phi(self):
        "Retorna um vetor com a temperatura de todas as celulas, inclusive as"
        "condicoes de contorno."
        cpoints = len(self.cells)
        all_T = np.zeros(cpoints+2)

        all_T[0] = self.lbc
        for i in range(cpoints):
            all_T[i+1] = self[i].phi
        all_T[cpoints+1] = self.rbc

        return all_T

    def set_phi(self, index, _phi):
        "Seta a temperatura da celula index."
        self[index].phi = _phi

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

    def print_phi(self):
        "Imprime na tela uma lista com a temperatura de todas as celulas"
        cpoints = len(self.cells)
        for i in range(cpoints):
            print self[i].phi

    def plot_T(self):
        "Plota um grafico com os dados da Temperatura"
        plt.figure()
        x = self.get_all_x()
        y = self.get_all_phi()

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

    grid = Grid(100, 500, 0)
    model = Model()

    for i in range(ncells):
        grid.add_cell(10e-3, 1000, 100, 0.1)

    #Resolvendo o Sistema
    A = np.matrix(model.build_matrix(grid))
    b = np.matrix(model.build_coef_vector(grid))
    x = (A.I*b.T).A1

    #Atribuindo os resultados
    for i in range(ncells):
        grid.set_phi(i, x[i])

    grid.print_phi()

    grid.plot_T()
