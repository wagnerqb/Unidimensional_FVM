# -*- coding: iso-8859-1 -*-
"""
@author: Bismarck Gomes Souza Jr e Wagner Queiroz Barros
@date:   Tue Apr 15 10:36:55 2014
@email:  bismarckgomes@gmail.com e wagnerqb@gmail.com
@brief:  Classe grid. Armazena as células.
"""
from __future__ import division
from Cell import *


class Grid():
    "Classe Grid genérica."
    #Tipos de condição de contorno
    PERIODIC_BC = -1    # Condição de contorno periódica
    VELOCITY_BC = 1     # Condição de contorno de velocidade
    PRESSURE_BC = 0     # Condição de contorno de pressão

    def __init__(self):
        "Inicializando variáveis."
        self.n = 0                  # Numero de células
        self.cells = []             # Inicializando vetor de células
        self.rbc_t = self.PERIODIC_BC    # Condição de contorno periódica
        self.lbc_t = self.PERIODIC_BC    # Condição de contorno periódica

    def __getitem__(self, index):
        "Sobrecarga do operador [ ]."
        return self.cells[index]

    def set_boundaries(self, lbc_t, lbc, rbc_t, rbc):
        """Condições de contorno.
        lbc_t: Left Boundary Type (-1:periodica, 0:pressure, 1:velocity)
        lbc: Left Boundary Value
        rbc_t: Left Boundary Type (-1:periodica, 0:pressure, 1:velocity)
        rbc: Right Boundary Value
        """
        self.lbc_t = lbc_t          # Left Boundary Type
        self.lbc = lbc              # Left Boundary Value
        self.rbc_t = rbc_t          # Right Bounday Type
        self.rbc = rbc              # Right Boundary Value

    def add_cell(self, cell):
        "Adiciona uma celula no grid."
        self.cells.append(cell)
        self.n += 1

    def msrc(self, index):
        "Fonte de massa da célula."
        return self[index].msrc

    #============================= ÁREA TRANSVERSAL ==========================#
    def A(self, index):
        "Area no centro da célula."
        if (index < 0):
            # Condição de contorno esquerda
            if self.lbc_t == self.PERIODIC_BC:
                # Condição de contorno periódica
                A = self[-1].A
            else:
                # Extrapolando a propriedade
                A = self.extrapolation_l(self[0].A, self[1].A)

        elif (index > self.n-1):
            # Condição de contorno direita
            if self.rbc_t == self.PERIODIC_BC:
                # Condição de contorno periódica
                A = self[0].A
            else:
                # Extrapolando a propriedade
                A = self.extrapolation_r(self[-2].A, self[-1].A)
        else:
            A = self[index].A

        return A

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

    #=========================== COMPRIMENTO (DX) ============================#
    def dx(self, index):
        "Delta_x no centro da celula index."
        if (index < 0):
            # Condição de contorno esquerda
            if self.lbc_t == self.PERIODIC_BC:
                # Condição de contorno periódica
                dx = self[-1].dx
            else:
                dx = self[0].dx

        elif (index > self.n-1):
            # Condição de contorno direita
            if self.rbc_t == self.PERIODIC_BC:
                # Condição de contorno periódica
                dx = self[0].dx
            else:
                dx = self[-1].dx

        else:
            dx = self[index].dx

        return dx

    def dx_r(self, index):
        "Delta_x no centro da celula direita."
        return self.dx(index+1)

    def dx_l(self, index):
        "Delta_x no centro da celula esquerda."
        return self.dx(index-1)

    def dx_rh(self, index):
        "Delta_x na face direita da célula."
        return (self.dx_r(index) + self.dx(index))/2

    def dx_lh(self, index):
        "Delta_x na face esquerda da célula."
        return (self.dx(index) + self.dx_l(index))/2

    #============================ EXTRAPOLAÇÃO ===============================#
    def extrapolation(self, x1, y1, x2, y2, x):
        """Função que interpola uma reta que passa pelos pontos (x1,y1) e
        (x2,y2) no ponto x."""
        return y1 + (x-x1)*(y2-y1)/(x2-x1)

    def extrapolation_r(self, y1, y2):
        "Função que extrapola uma propriedade y no bloco a direita."
        x1 = 0
        x2 = self.dx_lh(self.n-1)
        x = x2 + self.dx_rh(self.n-1)
        return self.extrapolation(x1, y1, x2, y2, x)

    def extrapolation_rrh(self, y1, y2):
        """Função que extrapola uma propriedade y na face direita do vizinho a
        direita."""
        x1 = 0
        x2 = self.dx_lh(self.n-1)
        x = x2 + self.dx_rh(self.n-1) + 0.5*self.dx_r(self.n-1)
        return self.extrapolation(x1, y1, x2, y2, x)

    def extrapolation_l(self, y1, y2):
        "Função que extrapola uma propriedade y no bloco a esquerda."
        x2 = 0
        x1 = -self.dx_rh(0)
        x = x1 - self.dx_lh(0)
        return self.extrapolation(x1, y1, x2, y2, x)

    def extrapolation_lh(self, y1, y2):
        """Função que extrapola uma propriedade y na face direta do vizinho a
        esquerda."""
        x2 = 0
        x1 = -self.dx_rh(0)
        x = x1 - self.dx(0)*0.5
        return self.extrapolation(x1, y1, x2, y2, x)


