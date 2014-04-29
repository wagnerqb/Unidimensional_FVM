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
from DiscretizationWell_DC import *


# Pipe properties
A = 0.1
dx = 10
ncells = 5
theta = 0

# Numerical Parameters
t0 = 0
dt = 1
nt = 5

# Fluid properties
T = 10                          # Temperatura
MM = 14                         # Massa molar
msrc = 0.                       # Mass Source term per volume unity
fsrc = 0.                       # Termo fonte da QM

# Initial properties
v_ini = 0                       # Initial Condition for v
p_ini = 10                       # Initial Condition for p

# Boundary Condition
# Left Boundary Condtion (Velocidade na entrada: v_ini)
lbc_t = 1                       # LBC Type (0 - Pressure / 1 - Velocity)
lbc = 5                         # LBC Value

# Right Boundary Condtion (Pressão na saida: p_ini)
rbc_t = 0                       # RBC Type (0 - Pressure / 1 - Velocity)
rbc = 10                         # RBC Value

#Fluido
fluid = FluidIdeal(MM)

# Creating Grid
grid = GridWell(fluid)
grid.set_boundaries(lbc_t, lbc, rbc_t, rbc)

for i in range(ncells):
    #Criando grid
    cell = CellWell(A, dx, fluid, v_ini, p_ini, T, theta, msrc)
    grid.add_cell(cell)

# Creating Model
model = DiscretizationWell_DC()


# Gráfico de Velocidade
plt.subplot(211)
plt.title('Velocidade')
plt.grid()

# Gráfico de Pressao
plt.subplot(212)
plt.title(u'Pressão')
plt.grid()

print '=============='
print ' Tempo  | IT '
print '--------+----'
for i in range(nt):
    t = t0 + i*dt

    it = model.iterate_t(grid, dt)

    print '{0:07g} | {1:02d}'.format(t, it)

    p = grid.get_p()
    v_rh = grid.get_v_rh()

#    # Gráfico de Velocidade
#    plt.subplot(211)
#    plt.plot(v_rh, label='Numerico')
#
#    # Gráfico de Pressao
#    plt.subplot(212)
#    plt.plot(p, label='Numerico')

print '==============\n'

print 'p:', p
print 'v:', v_rh

# Gráfico de Velocidade
plt.subplot(211)
plt.plot(v_rh, label='Numerico')

# Gráfico de Pressao
plt.subplot(212)
plt.plot(p, label='Numerico')

plt.show()
