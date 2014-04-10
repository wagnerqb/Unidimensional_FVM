# -*- coding: iso-8859-1 -*-
"""
@author: Wagner Queiroz Barros
@date:   Thu Apr 10 10:06:22 2014
@email:  wagnerqb@gmail.com
@brief:  Interpolation Classes
"""


class Interpolation():
    "Classe de Interpola��o pelo Esquema Centrado."

    #Fun��es de Interpola��o
    def center_scheme(self, leftprop, rightprop):
        "Classe interpola utilizando diferencas centrais." 
        return (rightprop + leftprop)/2

    #Fun��es para Calcular Derivadas
    def right_derivative(self, leftprop, rightprop, leftx, rightx):
        return (rightprop - leftprop)/(rightx/2 + leftx/2)
