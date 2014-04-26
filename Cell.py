# -*- coding: iso-8859-1 -*-
"""
@author: Bismarck Gomes Souza Jr e Wagner Queiroz Barros
@date:   Wed Apr 09 14:53:01 2014
@email:  bismarckgomes@gmail.com e wagnerqb@gmail.com
@brief:  Classe célula. Armazena as propriedades dos fluidos e do tubo.
"""


class CellWell():
    "Classe celula que será utilizada nas simulações no poço."

    def __init__(self, A, dx, rho, mu, v_rh, p):
        self.A = A          # Área transversal
        self.dx = dx        # Comprimento
        self.rho = rho      # Densidade
        self.mu = mu        # Viscosidade
        self.v_rh = v_rh    # Velocidade na face direita
        self.p = p          # Pressão no centro
