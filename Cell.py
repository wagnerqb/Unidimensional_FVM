# -*- coding: iso-8859-1 -*-
"""
@author: Bismarck Gomes Souza Jr e Wagner Queiroz Barros
@date:   Wed Apr 09 14:53:01 2014
@email:  bismarckgomes@gmail.com e wagnerqb@gmail.com
@brief:  Classe c�lula. Armazena as propriedades dos fluidos e do tubo.
"""


class Cell():
    "Classe celula gen�rica."

    def __init__(self, A, dx, msrc):
        self.A = A                  # �rea transversal
        self.dx = dx                # Comprimento
        self.msrc = msrc            # Fonte de massa


class CellWell(Cell):
    "Classe celula que ser� utilizada nas simula��es no po�o."

    def __init__(self, A, dx, rho, v_rh, p, msrc):
        Cell.__init__(self, A, dx, msrc)  # Inicializa �rea e comprimento
        self.rho = rho              # Densidade
        self.v_rh = v_rh            # Velocidade na face direita
        self.p = p                  # Press�o no centro
