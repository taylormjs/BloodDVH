#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  1 14:57:02 2017

@author: taylorsorenson,zhaoqiwang
"""

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
from readDoses import *

def plot_2d(matrix):
    "plot a 2d matrix"
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    y_dim = len(matrix)
    x_dim = len(matrix[0])
    X = np.arange( x_dim)
    Y = np.arange(y_dim)
    X, Y = np.meshgrid(X, Y)
    plt.title('Dose plot for ' +  str(matrix) )
    plt.xlabel('z positions')
    plt.ylabel('y positions')
    plt.xlim([x_min,x_max]) 
    plt.ylim([y_min,y_max])
    Z = matrix
    surf = ax.plot_surface(X, Y, Z, cmap=cm.cool,
                       linewidth=0, antialiased=False)
    plt.show()
    
def plot_bloods_3d(bloods, c = 'r',m ='o'):
    '''plot the postion of the blood in 3d'''
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    xs = []
    ys = []
    zs = []
    for i in bloods:
        (x,y,z) = i.get_position()
        xs.append(x)
        ys.append(y)
        zs.append(z)
    ax.scatter(xs, ys, zs, c=c, marker=m)

    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')
    plt.show()
