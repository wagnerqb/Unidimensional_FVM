# -*- coding: iso-8859-1 -*-
"""
@author: Bismarck Gomes Souza Jr e Wagner Queiroz Barros
@date:   Sat Apr 26 12:08:00 2014
@email:  bismarckgomes@gmail.com e wagnerqb@gmail.com
@brief:  Classe do fluido que escoa na tubulação.
"""

from __future__ import division

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


class FluidIdeal():
    "Classe Fluido ideal."

    R = 8.31144621  # Consante universal dos gases

    def __init__(self, MM):
        self.R_M = self.R/MM

    def rho(self, p, T):
        "Densidade do fluido na pressão p, e temperatura T."
        return p/(T*self.R_M)

    def drho_dp(self, p, T):
        """Derivada da densidade do fluido em relação a pressão,
        na pressão p, e temperatura T."""
        return 1./(T*self.R_M)

if __name__ == '__main__':

    print "TESTE"
