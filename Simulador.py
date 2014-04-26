# -*- coding: iso-8859-1 -*-
"""
@author: Bismarck Gomes Souza Jr e Wagner Queiroz Barros
@date:   Fri Apr 11 13:24:24 2014
@email:  bismarckgomes@gmail.com e wagnerqb@gmail.com
@brief:  Caso 1 
        - Escoamento monofásico
        - Fluido monocomponente
        - Sem trocas de calor com a parede
        - Escoamento isotérmico
        - Sem atrito com a parede
        - Sem gravidade
        - Fluido incomporessível
"""
from __future__ import division
import matplotlib.pyplot as plt
from Cell import CellWell
from GridWell import GridWell
from Fluid import *
from DiscretizationWell_COUPLE_DC import *


# Pipe properties
A = 0.1
dx = 0.1
ncells = 20

# Fluid properties
rho = 100                       # Fluid density
msrc = 0.                       # Mass Source term per volume unity
fsrc = 0.                       # Termo fonte da QM

# Initial properties
v_ini = 1                       # Initial Condition for v
p_ini = 0                       # Initial Condition for p

# Boundary Condition
# Left Boundary Condtion (Velocidade na entrada: v_ini)
lbc_t = 1                       # LBC Type (0 - Pressure / 1 - Velocity)
lbc = 1                         # LBC Value

# Right Boundary Condtion (Pressão na saida: p_ini)
rbc_t = 0                       # RBC Type (0 - Pressure / 1 - Velocity)
rbc = 0                         # RBC Value

#Fluido
fluid = FluidIncompressible(rho)

# Creating Grid
grid = GridWell(fluid)
grid.set_boundaries(lbc_t, lbc, rbc_t, rbc)

# Creating Model
model = DiscretizationWell_COUPLE_DC()

for i in range(ncells):
    #Criando grid
    cell = CellWell(A, dx, fluid, v_ini, p_ini, msrc)
    grid.add_cell(cell)


#Iterações
model.iterate_x(grid)
new_p = grid.get_p()
new_v = grid.get_v_rh()


# Gráfico de Velocidade
plt.subplot(211)
plt.plot(new_v, 'ko', label='Numerico')
#plt.plot(v_real, 'b', label='Real')
plt.title('Velocidade')
plt.grid()

# Gráfico de Pressao
plt.subplot(212)
plt.plot(new_p, 'ko', label='Numerico')
#plt.plot(p_real, 'b', label='Real')
plt.title(u'Pressão')
plt.grid()
plt.show()
