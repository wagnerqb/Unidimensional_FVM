# -*- coding: iso-8859-1 -*-
"""
@author: Bismarck Gomes Souza Jr e Wagner Queiroz Barros
@date:   Wed Apr 09 14:53:01 2014
@email:  bismarckgomes@gmail.com e wagnerqb@gmail.com
@brief:  Classe c�lula. Armazena as propriedades dos fluidos e do tubo.
"""
from Fluid import *


class Cell():
    "Classe celula gen�rica."

    def __init__(self, A, dx, theta, msrc):
        self.A = A                  # �rea transversal
        self.dx = dx                # Comprimento
        self.msrc = msrc            # Fonte de massa
        self.theta = theta                  # Inclina��o no centro


class CellWell(Cell):
    "Classe celula que ser� utilizada nas simula��es no po�o."

    def __init__(self, A, dx, fluid, v_rh, p, T, theta, msrc):
        Cell.__init__(self, A, dx, theta, msrc)
        self.fluid = fluid                  # Fluido
        self.v_rh = v_rh                    # Velocidade na face direita
        self.p = p                          # Press�o no centro
        self.T = T                          # Temperatura no centro
        self.p_old = p              # TODO: Melhorar condi��o inicial
        self.v_rh_old = v_rh        # TODO: Melhorar condi��o inicial

    def rho(self):
        "Densidade do fluido � press�o p e temperatura T."
        return self.fluid.rho(self.p, self.T)

    def rho_old(self):
        """Densidade do fluido na press�o e temperatura
        do passo de tempo anterior"""
        return self.fluid.rho(self.p_old, self.T)  # TODO: temperatura OLD
