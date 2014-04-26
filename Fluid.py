# -*- coding: iso-8859-1 -*-
"""
@author: Bismarck Gomes Souza Jr e Wagner Queiroz Barros
@date:   Sat Apr 26 12:08:00 2014
@email:  bismarckgomes@gmail.com e wagnerqb@gmail.com
@brief:  Classe do fluido que escoa na tubulação.
"""


class Fluid():
    "Classe Fluido."

    def __init__(self):
        pass

    def rho(self, p, T):
        pass


class FluidIncompressible():
    "Classe Fluido."

    def __init__(self, rho):
        self.rho_ = rho

    def rho(self, p=None, T=None):
        "Densidade do fluido na pressão p, e temperatura T"
        return self.rho_

    def drho_dp(self, p, T=None):
        """Derivada da Densidade do fluido em relação a pressão,
        na pressão p, e temperatura T"""
        return 0
