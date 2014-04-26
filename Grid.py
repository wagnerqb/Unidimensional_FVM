# -*- coding: iso-8859-1 -*-
"""
@author: Bismarck Gomes Souza Jr e Wagner Queiroz Barros
@date:   Tue Apr 15 10:36:55 2014
@email:  bismarckgomes@gmail.com e wagnerqb@gmail.com
@brief:  Classe grid. Armazena as c�lulas.
"""
from __future__ import division
from Cell import *


class Grid():
    "Classe Grid gen�rica."
    #Tipos de condi��o de contorno
    PERIODIC_BC = -1    # Condi��o de contorno peri�dica
    VELOCITY_BC = 1     # Condi��o de contorno de velocidade
    PRESSURE_BC = 0     # Condi��o de contorno de press�o

    def __init__(self):
        "Inicializando vari�veis."
        self.n = 0                  # Numero de c�lulas
        self.cells = []             # Inicializando vetor de c�lulas
        self.rbc_t = self.PERIODIC_BC    # Condi��o de contorno peri�dica
        self.lbc_t = self.PERIODIC_BC    # Condi��o de contorno peri�dica

    def __getitem__(self, index):
        "Sobrecarga do operador [ ]."
        return self.cells[index]

    def set_boundaries(self, lbc_t, lbc, rbc_t, rbc):
        """Condi��es de contorno.
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
        "Fonte de massa da c�lula."
        return self[index].msrc

    #============================= �REA TRANSVERSAL ==========================#
    def A(self, index):
        "Area no centro da c�lula."
        if (index < 0):
            # Condi��o de contorno esquerda
            if self.lbc_t == self.PERIODIC_BC:
                # Condi��o de contorno peri�dica
                A = self[-1].A
            else:
                # Extrapolando a propriedade
                A = self.extrapolation_l(self[0].A, self[1].A)

        elif (index > self.n-1):
            # Condi��o de contorno direita
            if self.rbc_t == self.PERIODIC_BC:
                # Condi��o de contorno peri�dica
                A = self[0].A
            else:
                # Extrapolando a propriedade
                A = self.extrapolation_r(self[-2].A, self[-1].A)
        else:
            A = self[index].A

        return A

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

    #=========================== COMPRIMENTO (DX) ============================#
    def dx(self, index):
        "Delta_x no centro da celula index."
        if (index < 0):
            # Condi��o de contorno esquerda
            if self.lbc_t == self.PERIODIC_BC:
                # Condi��o de contorno peri�dica
                dx = self[-1].dx
            else:
                dx = self[0].dx

        elif (index > self.n-1):
            # Condi��o de contorno direita
            if self.rbc_t == self.PERIODIC_BC:
                # Condi��o de contorno peri�dica
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
        "Delta_x na face direita da c�lula."
        return (self.dx_r(index) + self.dx(index))/2

    def dx_lh(self, index):
        "Delta_x na face esquerda da c�lula."
        return (self.dx(index) + self.dx_l(index))/2

    #============================ EXTRAPOLA��O ===============================#
    def extrapolation(self, x1, y1, x2, y2, x):
        """Fun��o que interpola uma reta que passa pelos pontos (x1,y1) e
        (x2,y2) no ponto x."""
        return y1 + (x-x1)*(y2-y1)/(x2-x1)

    def extrapolation_r(self, y1, y2):
        "Fun��o que extrapola uma propriedade y no bloco a direita."
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
        "Fun��o que extrapola uma propriedade y no bloco a esquerda."
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


