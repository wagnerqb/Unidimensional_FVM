# -*- coding: iso-8859-1 -*-
"""
@author: Bismarck Gomes Souza Jr e Wagner Queiroz Barros
@date:   Wed Apr 09 14:53:01 2014
@email:  bismarckgomes@gmail.com e wagnerqb@gmail.com
@brief:  Classe célula. Armazena as propriedades dos fluidos e do tubo.
"""


class Cell():
    "Classe celula genérica."

    def __init__(self, A, dx, msrc):
        self.A = A                  # Área transversal
        self.dx = dx                # Comprimento
        self.msrc = msrc            # Fonte de massa


class CellWell(Cell):
    "Classe celula que será utilizada nas simulações no poço."

    def __init__(self, A, dx, rho, v_rh, p, msrc):
        Cell.__init__(self, A, dx, msrc)  # Inicializa área e comprimento
        self.rho = rho              # Densidade
        self.v_rh = v_rh            # Velocidade na face direita
        self.p = p                  # Pressão no centro
