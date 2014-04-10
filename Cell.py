# -*- coding: iso-8859-1 -*-
"""
@author: Wagner Queiroz Barros
@date:   Wed Apr 09 14:53:01 2014
@email:  wagnerqb@gmail.com
@brief:  Cell class, used to store numerical data
"""
from __future__ import division


class Cell():
    "Classe Celula."

    def __init__(self, A, k, T, dx):
        #Atributos
        self.A = A          # Area
        self.k = k          # Kappa
        self.T = T          # Temperatura
        self.dx = dx        # Delta x


if __name__ == '__main__':
    c = Cell(1, 2, 3, 4)
    print c.A
    print c.k
    print c.T
    print c.dx

    print c.__module__
