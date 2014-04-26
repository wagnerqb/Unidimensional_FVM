# -*- coding: iso-8859-1 -*-
"""
@author: Bismarck Gomes Souza Jr e Wagner Queiroz Barros
@date:   Fri Apr 11 13:24:24 2014
@email:  bismarckgomes@gmail.com e wagnerqb@gmail.com
@brief:  Exemplo do Bocal Convergente.
"""
from __future__ import division
import matplotlib.pyplot as plt
from Grid import GridWell
from DiscretiztionWell_COUPLE_CDS import *


# Pipe properties
dx = 0.1
ncells = 20
a_ = 1                          # Area de saida
b_ = 2                          # Area de entrada

# Fluid properties
rho = 100                       # Fluid density
mu = 0                          # Fluid viscosity
msrc = 0.                       # Mass Source term per volume unity

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

# Creating Grid
grid = GridWell(lbc_t, lbc, rbc_t, rbc, msrc)

# Creating Model
model = DiscretiztionWell_COUPLE_CDS()

#Soluções reais
v_real = []
p_real = []
for i in range(ncells):
    #Criando grid
    A = b_ - (b_-a_)*(.5+i)/ncells
    grid.add_cell(A, dx, rho, mu, v_ini, p_ini)

    #Calculando as soluções analíticas
    A_rh = b_ - (b_-a_)*(1+i)/ncells
    v_real.append(b_*lbc/A_rh)
    p_real.append(rbc+0.5*rho*(lbc*b_)**2*(1/a_**2-1/A**2))


#Iterações
it, erro = model.iterate_x(grid)
new_p = grid.get_p()
new_v = grid.get_v_rh()

## Resultado
#print '\n\nErro (L2):', erro
#print 'Iterações:', it

# Gráfico de Velocidade
plt.subplot(211)
plt.plot(new_v, 'ko', label='Numerico')
plt.plot(v_real, 'b', label='Real')
plt.title('Velocidade')
plt.grid()

# Gráfico de Pressao
plt.subplot(212)
plt.plot(new_p, 'ko', label='Numerico')
plt.plot(p_real, 'b', label='Real')
plt.title(u'Pressão')
plt.grid()
plt.show()
