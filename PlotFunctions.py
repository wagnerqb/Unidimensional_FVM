# -*- coding: iso-8859-1 -*-
"""
@author: Bismarck Gomes Souza Junior
@date:   Tue Apr 29 15:54:19 2014
@email:  bismarckjunior@outlook.com
@brief:  
"""
import matplotlib.pyplot as plt


def create_figure_iterative(grid, min_v=None, max_v=None, min_p=None, max_p=None):
    '''Inicializa a criação dos eixos. Necessário para o plot iterativo.'''
    fig = plt.figure()

    # Gráfico de Velocidade
    ax_v = fig.add_subplot(211)
    plt_v = ax_v.plot(grid.get_v_rh())[0]
    ax_v.set_title('Velocidade')
    if min_v != None and max_v != None:
        ax_v.set_ylim(min_v, max_v)
    ax_v.grid()

    # Gráfico de Pressao
    ax_p = fig.add_subplot(212)
    plt_p = ax_p.plot(grid.get_p())[0]
    ax_p.set_title(u'Pressão')
    if min_p != None and max_p != None:
        ax_p.set_ylim(min_p, max_p)
    ax_p.grid()

    plt.ion()
    plt.show()

    return plt_v, plt_p


def update_iterative_plot(plt_v, plt_p, x, v_rh, p):
    '''Atualiza o plot iterativo.'''
    # Gráfico de Velocidade
    plt.subplot(211)
#    plt_p.set_data(x, v_rh)
    plt_v.set_data([], [])
    plt_v = plt.plot(x, v_rh, 'b', lw=2)[0]

    # Gráfico de Pressao
    plt.subplot(212)
#    plt_p.set_data(x, p)
    plt_p.set_data([], [])
    plt_p = plt.plot(x, p, 'r',lw=2)[0]

    plt.pause(0.0001)

    plt.draw()

    return plt_v, plt_p
