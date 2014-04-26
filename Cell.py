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

    def __init__(self, A, dx, msrc):
        self.A = A                  # �rea transversal
        self.dx = dx                # Comprimento
        self.msrc = msrc            # Fonte de massa


class CellWell(Cell):
    "Classe celula que ser� utilizada nas simula��es no po�o."

    def __init__(self, A, dx, fluid, v_rh, p, msrc):
        Cell.__init__(self, A, dx, msrc)    # Inicializa �rea e comprimento
        self.fluid = fluid                  # Fluido
        self.v_rh = v_rh                    # Velocidade na face direita
        self.p = p                          # Press�o no centro
        self.p_old = p              # TODO: Melhorar condi��o inicial
        self.v_rh_old = v_rh        # TODO: Melhorar condi��o inicial

    def rho(self, p, T):
        "Densidade do fluido � press�o p e temperatura T."
        return self.fluid.rho(p, T)
