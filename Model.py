# -*- coding: iso-8859-1 -*-
"""
@author: Wagner Queiroz Barros
@date:   Wed Apr 09 15:35:41 2014
@email:  wagnerqb@gmail.com
@brief:  Class Storing all models using to discretize PDEs
"""

#from Cell import *
from __future__ import division

class Model():
    "Classe de Modelo"
    
    #Funções de Derivada
    def RightDerivative(self, leftprop, rightprop, leftx, rightx):
        return (rightprop - leftprop)/(rightx/2 + leftx/2)
        
    #Funções de Interpolação
    def CenterScheme(self, leftprop, rightprop) :
        #Classe interpola utilizando diferencas centrais
        return (rightprop - leftprop)/2
        
#    def ExplicitIteration(self, grid):
#        for g