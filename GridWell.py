# -*- coding: iso-8859-1 -*-
"""
@author: Bismarck Gomes Souza Jr e Wagner Queiroz Barros
@date:   Sat Apr 26 12:22:52 2014
@email:  bismarckgomes@gmail.com e wagnerqb@gmail.com
@brief:  Classe Grid que armazenar� as c�lulas do po�o.
"""
from Grid import Grid
import numpy as np


class GridWell(Grid):
    "Classe Grid que armazenar� as c�lulas do po�o."

    def __init__(self):
        "Inicializando vari�veis."
        Grid.__init__(self)

    #========================== MASSA ESPEC�FICA =============================#
    def rho(self, index, p=None, T=None):
        "Densidade no centro da c�lula."
        if (index < 0):
            # Condi��o de contorno esquerda
            if self.lbc_t == self.PERIODIC_BC:
                # Condi��o de contorno peri�dica
                rho = self[-1].rho(p, T)
            else:
                rho = self[0].rho(p, T)

        elif (index > self.n-1):
            # Condi��o de contorno direita
            if self.rbc_t == self.PERIODIC_BC:
                # Condi��o de contorno peri�dica
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
        "Velocidade na face direita da c�lula."
        if (index < 0):
            # Condi��o de contorno esquerda
            if (self.lbc_t == self.VELOCITY_BC):
                # Condi��o de contorno de velocidade
                v_rh = self.lbc

            elif (self.lbc_t == self.PERIODIC_BC):
                # Condi��o de contorno peri�dica
                v_rh = self[-1].v_rh

            else:
                # Para outras condi��es de contorno a velocidade � extrapolada
                A_lh = self.extrapolation_lh(self[0].A, self[1].A)
                v_rh = self[0].v_rh*self.A_rh(0)/A_lh

        elif (index > self.n-1):
            # Condi��o de contorno direita
            if (self.rbc_t == self.VELOCITY_BC):
                # Condi��o de contorno de velocidade
                v_rh = self.rhc

            elif (self.rbc_t == self.PERIODIC_BC):
                # Condi��o de contorno peri�dica
                v_rh = self[0].v_rh

            else:
                # Condi��o de contorno de press�o
                A_rrh = self.extrapolation_rrh(self[-2].A, self[-1].A)

                ##Area interpolada dentro e fora do tubo
                #v_rh = self[-1].v_rh*self.A_rh(self.n-1)/A_rrh

                ##Area interpolada dentro do tubo e considerada constante fora
                v_rh = self[-1].v_rh*self.A_r(self.n-1)/A_rrh

        else:
            v_rh = self[index].v_rh

        return v_rh

    def v_lh(self, index):
        "Velocidade na face esquerda da c�lula."
        return self.v_rh(index-1)

    def v(self, index):
        "Velocidade no centro da c�lula (CDS Interpolation Method)."
        return (self.v_lh(index) + self.v_rh(index))/2

    def v_r(self, index):
        "Velocidade no centro da c�lula direita (CDS Method)."
        return self.v(index + 1)

    def v_l(self, index):
        "Velocidade no centro da c�lula esquerda (CDS Method)."
        return self.v(index - 1)

    def set_v_rh(self, index, v_rh):
        "Atribui a velocidade na face direita da c�lula."
        self[index].v_rh = v_rh

    def get_v_rh(self):
        "Retorna um vetor com a velocidade de todas as celulas."
        v_rh = []
        for i in range(self.n):
            v_rh.append(self[i].v_rh)
        return np.array(v_rh)

    #============================== PRESS�O ==================================#
    def p(self, index):
        "Pess�o no centro da c�lula."
        if (index < 0):
            # Condi��o de contorno esquerda
            if (self.lbc_t == self.PRESSURE_BC):
                # Condi��o de contorno de press�o
                p = self.lbc

            elif (self.lbc_t == self.PERIODIC_BC):
                # Condi��o de contorno peri�dica
                p = self[-1].p

            else:
                # Com extrapola��o
                p = self.extrapolation_l(self[0].p, self[1].p)

        elif index > self.n-1:
            # Condi��o de contorno direita
            if (self.rbc_t == self.PRESSURE_BC):
                # Condi��o de contorno de press�o
                p = self.rbc

            elif (self.rbc_t == self.VELOCITY_BC):
                # Condi��o de contorno peri�dica
                p = self[0].p

            else:
                # Para outras condi��es de contorno a press�o � extrapolada
                p = self.extrapolation_r(self[-2].p, self[-1].p)

        else:
            p = self[index].p

        return p

    def p_r(self, index):
        "Press�o no centro da c�lula direita."
        return self.p(index + 1)

    def p_l(self, index):
        "Press�o no centro da c�lula esquerda."
        return self.p(index - 1)

    def p_rh(self, index):
        "Press�o na face direita da c�lula."
        return (self.p_r(index) + self.p(index))/2

    def p_rl(self, index):
        "Press�o na face esquerda da c�lula."
        return (self.p(index) + self.p_l(index))/2

    def set_p(self, index, p):
        "Atribui a press�o no centro da c�lula"
        self[index].p = p

    def get_p(self):
        """Retorna um vetor com a press�o de todas as celulas."""
        p = []
        for i in range(self.n):
            p.append(self[i].p)
        return np.array(p)
