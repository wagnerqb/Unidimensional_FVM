# -*- coding: iso-8859-1 -*-
"""
@author: Bismarck Gomes Souza Jr e Wagner Queiroz Barros
@date:   Sat Apr 26 12:22:52 2014
@email:  bismarckgomes@gmail.com e wagnerqb@gmail.com
@brief:  Classe Grid que armazenará as células do poço.
"""
from Grid import Grid
import numpy as np


class GridWell(Grid):
    "Classe Grid que armazenará as células do poço."

    def __init__(self):
        "Inicializando variáveis."
        Grid.__init__(self)

    #========================== MASSA ESPECÍFICA =============================#
    def rho(self, index, p=None, T=None):
        "Densidade no centro da célula."
        if (index < 0):
            # Condição de contorno esquerda
            if self.lbc_t == self.PERIODIC_BC:
                # Condição de contorno periódica
                rho = self[-1].rho(p, T)
            else:
                rho = self[0].rho(p, T)

        elif (index > self.n-1):
            # Condição de contorno direita
            if self.rbc_t == self.PERIODIC_BC:
                # Condição de contorno periódica
                rho = self[0].rho(p, T)
            else:
                rho = self[-1].rho(p, T)

        else:
            rho = self[index].rho(p, T)

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
            if (self.lbc_t == self.VELOCITY_BC):
                # Condição de contorno de velocidade
                v_rh = self.lbc

            elif (self.lbc_t == self.PERIODIC_BC):
                # Condição de contorno periódica
                v_rh = self[-1].v_rh

            else:
                # Para outras condições de contorno a velocidade é extrapolada
                A_lh = self.extrapolation_lh(self[0].A, self[1].A)
                v_rh = self[0].v_rh*self.A_rh(0)/A_lh

        elif (index > self.n-1):
            # Condição de contorno direita
            if (self.rbc_t == self.VELOCITY_BC):
                # Condição de contorno de velocidade
                v_rh = self.rhc

            elif (self.rbc_t == self.PERIODIC_BC):
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
            if (self.lbc_t == self.PRESSURE_BC):
                # Condição de contorno de pressão
                p = self.lbc

            elif (self.lbc_t == self.PERIODIC_BC):
                # Condição de contorno periódica
                p = self[-1].p

            else:
                # Com extrapolação
                p = self.extrapolation_l(self[0].p, self[1].p)

        elif index > self.n-1:
            # Condição de contorno direita
            if (self.rbc_t == self.PRESSURE_BC):
                # Condição de contorno de pressão
                p = self.rbc

            elif (self.rbc_t == self.VELOCITY_BC):
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
