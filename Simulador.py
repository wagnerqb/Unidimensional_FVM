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
nt = 30

# Fluid properties
T = 10                          # Temperatura
MM = 14                         # Massa molar
msrc = 0.                       # Mass Source term per volume unity
fsrc = 0.                       # Termo fonte da QM

# Initial properties
v_ini = 8                       # Initial Condition for v
p_ini = 8                       # Initial Condition for p

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



ax_v = plt.subplot(211)
plt_v = ax_v.plot(grid.get_v_rh())[0]
plt.title('Velocidade')
plt.ylim(0, lbc*2)
plt.grid()

# Gráfico de Pressao
ax_p = plt.subplot(212)
plt_p = ax_p.plot(grid.get_p())[0]
plt.title(u'Pressão')
plt.grid()
plt.ylim(0, rbc*2)

plt.ion()
plt.show()

print '=============='
print ' Tempo  | IT '
print '--------+----'
for i in range(nt):
    t = t0 + i*dt

    it = model.iterate_t(grid, dt)

    print '{0:07g} | {1:02d}'.format(t, it)

    p = grid.get_p()
    v_rh = grid.get_v_rh()
    x = range(ncells)

    # Gráfico de Velocidade
    plt.subplot(211)
    plt_v.set_data([], [])
    plt_v = plt.plot(x, v_rh, 'b', lw = 2)[0]

    # Gráfico de Pressao
    plt.subplot(212)
    plt_p.set_data([], [])
    plt_p = plt.plot(x, p, 'r',lw = 2)[0]

    plt.pause(0.001)

    plt.draw()

print '==============\n'

print 'p:', p
print 'v:', v_rh

raw_input()


