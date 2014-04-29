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
"""
from __future__ import division

from Cell import CellWell
from GridWell import GridWell
from Fluid import *
from DiscretizationWell_DC import *
import PlotFunctions as plot


    
# Pipe properties
A = 0.1
dx = 5
ncells = 50
theta = 0

# Numerical Parameters
t0 = 0
dt = 0.7
nt = 1000000

# Fluid properties
T = 300                         # Temperatura
MM = 14                         # Massa molar
msrc = 0.                       # Mass Source term per volume unity
fsrc = 0.                       # Termo fonte da QM

# Initial properties
v_ini = 0                       # Initial Condition for v
p_ini = 10                       # Initial Condition for p

# Boundary Condition
# Left Boundary Condtion (Velocidade na entrada: v_ini)
lbc_t = 1                       # LBC Type (0 - Pressure / 1 - Velocity)
lbc = 2                         # LBC Value

# Right Boundary Condtion (Pressão na saida: p_ini)
rbc_t = 0                       # RBC Type (0 - Pressure / 1 - Velocity)
rbc = 10                        # RBC Value

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

# Iniciando o plot iterativo
plt_v, plt_p = plot.create_figure_iterative(grid, 0, 5, 0, 15)

print '==================='
print '  N  | Tempo  | IT '
print '-----+--------+----'
for i in range(nt):
    t = t0 + (i+1)*dt

    it = model.iterate_t(grid, dt)

    print '{0:04} | {1:07.3f} | {2:02d}'.format(i+1, t, it)

    p = grid.get_p()
    v_rh = grid.get_v_rh()
    x = range(ncells)

    # Atualizando o plot iterativo
    #if i % 10 == 0:
    plt_v, plt_p = plot.update_iterative_plot(plt_v, plt_p, x, v_rh, p)

print '==============\n'

print 'p:', p
print 'v:', v_rh

raw_input()
