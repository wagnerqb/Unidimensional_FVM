# -*- coding: iso-8859-1 -*-
"""
@author: Bismarck Gomes Souza Jr e Wagner Queiroz Barros
@date:   Tue Apr 15 10:36:55 2014
@email:  bismarckgomes@gmail.com e wagnerqb@gmail.com
@brief:  Classe grid. Armazena as células.
"""
from __future__ import division
import numpy as np
from Cell import *

#Tipos de condição de contorno
PERIODIC_BC = -1    # Condição de contorno periódica
VELOCITY_BC = 1     # Condição de contorno de velocidade
PRESSURE_BC = 0     # Condição de contorno de pressão


class Grid():
    "Classe Grid genérica."

    def __init__(self):
        "Inicializando variáveis."
        self.n = 0                  # Numero de células
        self.cells = []             # Inicializando vetor de células
        self.rbc_t = PERIODIC_BC    # Condição de contorno periódica
        self.lbc_t = PERIODIC_BC    # Condição de contorno periódica

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
            if self.lbc_t == PERIODIC_BC:
                # Condição de contorno periódica
                A = self[-1].A
            else:
                # Extrapolando a propriedade
                A = self.extrapolation_l(self[0].A, self[1].A)

        elif (index > self.n-1):
            # Condição de contorno direita
            if self.rbc_t == PERIODIC_BC:
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
            if self.lbc_t == PERIODIC_BC:
                # Condição de contorno periódica
                dx = self[-1].dx
            else:
                dx = self[0].dx

        elif (index > self.n-1):
            # Condição de contorno direita
            if self.rbc_t == PERIODIC_BC:
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


class GridWell(Grid):
    "Classe Grid que armazenará as células do poço."

    def __init__(self):
        "Inicializando variáveis."
        Grid.__init__(self)

    #========================== MASSA ESPECÍFICA =============================#
    def rho(self, index):
        "Densidade no centro da célula."
        if (index < 0):
            rho = self[0].rho

        elif (index > self.n-1):
            rho = self[self.n-1].rho

        else:
            rho = self[index].rho

        return rho

    def rho_r(self, index):
        "Densidade no centro da celula a direita."
        return self.rho(index+1)

    def rho_l(self, index):
        "Densidade no centro da celula a esquerda."
        return self.rho(index-1)

    def rho_rh(self, index):
        "Densidade na face direita da celula."
        return (self.rho(index) + self.rho_r(index))/2

    def rho_lh(self, index):
        "Densidade na face esquerda da celula."
        return (self.rho_l(index) + self.rho(index))/2

    #============================ VELOCIDADE =================================#
    def v_rh(self, index):
        "Velocidade na face direita da célula."
        if (index < 0):
            # Condição de contorno esquerda
            if (self.lbc_t == VELOCITY_BC):
                # Condição de contorno de velocidade
                v_rh = self.lbc

            elif (self.lbc_t == PERIODIC_BC):
                # Condição de contorno periódica
                v_rh = self[-1].vrh

            else:
                # Para outras condições de contorno a velocidade é extrapolada
                A_lh = self.extrapolation_lh(self[0].A, self[1].A)
                v_rh = self[0].v_rh*self.A_rh(0)/A_lh

        elif (index > self.n-1):
            # Condição de contorno direita
            if (self.rbc_t == VELOCITY_BC):
                # Condição de contorno de velocidade
                v_rh = self.rhc

            elif (self.rbc_t == PERIODIC_BC):
                # Condição de contorno periódica
                v_rh = self[0].v_rh

            else:
                # Condição de contorno de pressão
                A_rrh = self.extrapolation_rrh(self[-2].A, self[-1].A)

                ##Area interpolada dentro e fora do tubo
                #v_rh = self[-1].v_rh*self.A_rh(self.n-1)/A_rrh

                ##Area interpolada dentro do tubo e considerada constante fora
                v_rh = self[-1].v_rh*self.A_r(self.n-1)/A_rrh

        else:
            v_rh = self[index].v_rh

        return v_rh

    def v_lh(self, index):
        "Velocidade na face esquerda da célula."
        return self.v_rh(index-1)

    def v(self, index):
        "Velocidade no centro da célula (CDS Interpolation Method)."
        return (self.v_lh(index) + self.v_rh(index))/2

    def v_r(self, index):
        "Velocidade no centro da célula direita (CDS Method)."
        return self.v(index + 1)

    def v_l(self, index):
        "Velocidade no centro da célula esquerda (CDS Method)."
        return self.v(index - 1)

    def set_v_rh(self, index, v_rh):
        "Atribui a velocidade na face direita da célula."
        self[index].v_rh = v_rh

    def get_v_rh(self):
        "Retorna um vetor com a velocidade de todas as celulas."
        v_rh = []
        for i in range(self.n):
            v_rh.append(self[i].v_rh)
        return np.array(v_rh)

    #============================== PRESSÃO ==================================#
    def p(self, index):
        "Pessão no centro da célula."
        if (index < 0):
            # Condição de contorno esquerda
            if (self.lbc_t == PRESSURE_BC):
                # Condição de contorno de pressão
                p = self.lbc

            elif (self.lbc_t == PERIODIC_BC):
                # Condição de contorno periódica
                p = self[-1].p

            else:
                # Com extrapolação
                p = self.extrapolation_l(self[0].p, self[1].p)

        elif index > self.n-1:
            # Condição de contorno direita
            if (self.rbc_t == PRESSURE_BC):
                # Condição de contorno de pressão
                p = self.rbc

            elif (self.rbc_t == VELOCITY_BC):
                # Condição de contorno periódica
                p = self[0].p

            else:
                # Para outras condições de contorno a pressão é extrapolada
                p = self.extrapolation_r(self[-2].p, self[-1].p)

        else:
            p = self[index].p

        return p

    def p_r(self, index):
        "Pressão no centro da célula direita."
        return self.p(index + 1)

    def p_l(self, index):
        "Pressão no centro da célula esquerda."
        return self.p(index - 1)

    def p_rh(self, index):
        "Pressâo na face direita da célula."
        return (self.p_r(index) + self.p(index))/2

    def p_rl(self, index):
        "Pressâo na face esquerda da célula."
        return (self.p(index) + self.p_l(index))/2

    def set_p(self, index, p):
        "Atribui a pressão no centro da célula"
        self[index].p = p

    def get_p(self):
        """Retorna um vetor com a pressão de todas as celulas."""
        p = []
        for i in range(self.n):
            p.append(self[i].p)
        return np.array(p)
