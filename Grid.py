# -*- coding: iso-8859-1 -*-
"""
@author: Bismarck Gomes Souza Jr e Wagner Queiroz Barros
@date:   Tue Apr 15 10:36:55 2014
@email:  bismarckgomes@gmail.com e wagnerqb@gmail.com
@brief:  Classe grid. Armazena as células.
"""
from __future__ import division
import numpy as np
from Cell import CellWell


class GridWell():
    "Classe Grid que armazenará as células."

    def __init__(self, lbc_t, lbc, rbc_t, rbc, msrc):
        self.cells = []             # Boundary (0-pressure, 1-velocity)
        self.lbc_t = lbc_t          # Left Boundary Type
        self.lbc = lbc              # Left Boundary Condition
        self.rbc_t = rbc_t          # Right Bounday Type
        self.rbc = rbc              # Right Boundary Condition
        self.msrc = msrc            # Mass Source Term
        self.n = 0

    def __getitem__(self, index):
        "Sobrecarga do operador [ ]."
        return self.cells[index]

    def add_cell(self, A, dx, rho, mu, v_r, p):
        "Adiciona uma celula no grid."
        self.cells.append(CellWell(A, dx, rho, mu, v_r, p))
        self.n += 1

    def A(self, index):
        "Area no centro da célula."
        if (index < 0):
            # Sem extrapolação
            #A = self[0].A

            # Com extrapolação
            A = self.extrapolation_l(self[0].A, self[1].A)

        elif (index > self.n-1):
            # Sem extrapolação
            #A = self[(len(self.cells)-1)].A

            # Com extrapolação
            A = self.extrapolation_r(self[-2].A, self[-1].A)
        else:
            A = self[index].A

        return A

    def dx(self, index):
        "Delta_x no centro da celula index."
        if (index < 0):
            dx = self[0].dx

        elif (index > self.n-1):
            dx = self[self.n-1].dx

        else:
            dx = self[index].dx

        return dx

    def rho(self, index):
        "Densidade no centro da célula."
        if (index < 0):
            rho = self[0].rho

        elif (index > self.n-1):
            rho = self[self.n-1].rho

        else:
            rho = self[index].rho

        return rho

    def mu(self, index):
        "Viscosidade no centro da celula."
        if (index < 0):
            return self[0].mu

        if (index > self.n-1):
            return self[self.n-1].mu

        return self[index].mu

    def v_rh(self, index):
        "Velocidade na face direita da célula."
        if (index < 0):
            if (self.lbc_t == 0):
                # Sem extrapolação
                #return self[0].v_rh

                #TIP: Não foi testado ########################################

                # Com extrapolação
                A_lh = self.extrapolation_lh(self[0].A, self[1].A)
                return self[0].v_r*self.A_rh(0)/A_lh
            else:
                return self.lbc
        if (index > self.n-1):
            if (self.rbc_t == 0):
                # Sem extrapolação
                #return self[n-1].v_rh

                # Com extrapolação
                A_rrh = self.extrapolation_rrh(self[-2].A, self[-1].A)
                ##Area interpolada dentro do tubo e considerada constante fora
                return self[-1].v_rh*self.A_r(self.n-1)/A_rrh

                ##Area interpolada dentro e fora do tubo
                # return self[-1].v_rh*self.A_rh(self.n-1)/A_rrh
            else:
                return self.rbc

        return self[index].v_rh

    def p(self, index):
        "Pessão no centro da célula."
        if (index < 0):
            if (self.lbc_t == 0):
                return self.lbc
            else:
                # Sem extrapolação
                #return self[0].p

                # Com extrapolação
                return self.extrapolation_l(self[0].p, self[1].p)

        if index > self.n-1:
            if (self.rbc_t == 0):
                return self.rbc
            else:
                # Sem extrapolação
                #return self[(len(self.cells)-1)].p

                # Com extrapolação
                return self.extrapolation_r(self[-2].p, self[-1].p)

        return self[index].p

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

    def mu_r(self, index):
        "Viscosidade no centro da célula a direita."
        return self.mu(index+1)

    def mu_l(self, index):
        "Viscosidade no centro da célula a esquerda."
        return self.mu(index-1)

    def mu_rh(self, index):
        "Viscosidade na face direita da célula."
        return (self.mu_r(index) + self.mu(index))/2

    def mu_lh(self, index):
        "Viscosidade na face esquerda da célula."
        return (self.mu(index) + self.mu_l(index))/2

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

    def set_v_rh(self, index, v_rh):
        "Atribui a velocidade na face direita da célula."
        self[index].v_rh = v_rh

    def get_p(self):
        """Retorna um vetor com a pressão de todas as celulas."""
        p = []
        for i in range(self.n):
            p.append(self[i].p)
        return np.array(p)

    def get_v_rh(self):
        "Retorna um vetor com a velocidade de todas as celulas."
        v_rh = []
        for i in range(self.n):
            v_rh.append(self[i].v_rh)
        return np.array(v_rh)

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

#    def msrc(self, index):
#        "Termo Fonte"
#        return self.msrc
#
#    def get_all_x(self):
#        "Retorna um vetor com a posicao x do centro de todas as celulas,"
#        "inclusive as condicoes de contorno."
#        cpoints = len(self.cells)
#        all_x = np.zeros(cpoints+2)
#
#        all_x[0] = 0
#        all_x[1] = self[0].dx/2
#        for i in range(2, cpoints+1):
#            all_x[i] = all_x[i-1] + self[i-1].dx/2 + self[i-2].dx/2
#
#        all_x[cpoints+1] = all_x[cpoints] + self[cpoints-1].dx/2
#        return all_x
#
#    def get_all_p(self):
#        """Retorna um vetor com a pressão de todas as celulas, inclusive as
#        condicoes de contorno."""
#        all_P = np.zeros(self.n+2)
#
#        all_P[0] = self.lbc
#        for i in range(self.n):
#            all_T[i+1] = self[i].p
#        all_T[self.n+1] = self.rbc
#
#        return all_P
