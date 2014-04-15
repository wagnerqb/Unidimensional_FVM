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

    def __init__(self, A, k, dx):
        #Atributos
        self.A = A          # Area
        self.k = k          # Kappa        
        self.dx = dx        # Delta x


class CellD(Cell):
    "Class Célula com propriedades difusivas."
    
    def __init__(self, A, k, dx, phi ):
        #Chamando o construtor de Cell
        Cell.__init__(self, A, k, dx)

        #Atributos
        self.phi = phi      # Propriedade Transportada

        
class CellCD(CellD):
    "Classe herdeira com propriedades convectivas."
    
    def __init__(self, A, k, dx, phi, rho, v):
        #Chamando o construtor de CellD
        CellD.__init__(self, A, k, dx, phi)

        #Atributos
        self.rho = rho      # Densidade
        self.v = v          # Velocidade


class CellFluid():
    "Classe Celula Fluida, Utlizada em Navier-Stokes"

    #Construtor com Atributos
    def __init__(self, A, dx, rho, mu, v_r, p):
        self.A = A
        self.dx = dx
        self.rho = rho
        self.mu = mu
        self.v_r = v_r
        self.p = p


if __name__ == '__main__':
    
    c = Cell(1, 2, 3)
    print c.A
    print c.k
    print c.dx
    print

    cd = CellD(1, 2, 3, 4)
    print cd.A
    print cd.k
    print cd.phi
    print cd.dx
    print
    
    ccd = CellCD(1, 2, 3, 4, 5, 6)
    print ccd.A
    print ccd.k
    print ccd.phi
    print ccd.dx
    print ccd.rho
    print ccd.v
    print

    cfld = CellFluid(1, 2, 3, 4, 5, 6)
    print cfld.A
    print cfld.dx
    print cfld.rho
    print cfld.mu
    print cfld.v_r
    print cfld.p
