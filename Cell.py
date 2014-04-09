# -*- coding: iso-8859-1 -*-
"""
@author: Wagner Queiroz Barros
@date:   Wed Apr 09 14:53:01 2014
@email:  wagnerqb@gmail.com
@brief:  Cell class, used to store numerical data
"""
from __future__ import division
from Model import *

class Cell():
    "Classe Celula."
    
    def __init__(self, A, kappa, T, deltax):
        
        #Atributos
        self.A = A
        self.kappa = kappa
        self.T = T
        self.deltax = deltax     
        
if __name__ == '__main__':
    c = Cell(1, 2, 3, 4)
    print c.A
    print c.kappa
    print c.T
    print c.deltax
    
    print c.__module__
