#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 11:09:47 2017

@author: taylorsorenson
"""

#Simvascular Note - make sure to site the given papers
#Updegrove, A., Wilson, N., Merkow, J., Lan, H., Marsden, A. L. and Shadden,
#S. C., SimVascular - An open source pipeline for cardiovascular simulation, 
#Annals of Biomedical Engineering (2017) 45: 525. DOI:10.1007/s10439-016-1762-8 (2017) 

#==============================================================================
# Data from Simvascular Simulation - Aorta and arteries
#==============================================================================
import numpy as np
import matplotlib.pyplot as plt


#==============================================================================
# Create velocity fields here: Take these in as parameters for each vector_field
# object made
#==============================================================================
x = np.linspace(0, 5, 50)
y = np.linspace(0, 5, 75)
z = np.linspace(0, 5, 100)
xx, yy, zz = np.meshgrid(x, y, z, indexing = 'ij') #or is indexing = 'xy' ?
vx = xx*0
vy = yy**2
vz = zz
#h = plt.contourf(x,y,vy)





def make_into_cell(x):
    '''Useful for map function - makes all the positions integer values
    '''
    if x == 0:
        return x
    else:
        f = float(x) #make the string into a float
        return int(f*100) #multiple the float by 100 to get a cell number


def create_velocity_field(filename):
    '''returns 3D array, denoting a velocity field
    Note - each field either has vx, vy, or vz, NOT all three
    '''    
    with open(filename, 'r') as file:
        velx, vely, velz, x_pos, y_pos, z_pos = [], [], [], [], [], []
        file.readline() #to skip the first line
        for line in file:
            entry = line.split(',')
            #pull and labels just the values we want - this may need to be changed
            #if more files are added
            x,y,z = tuple(map(make_into_cell,[entry[2],entry[3],entry[4]]))
            vx,vy,vz = tuple(map(float, [entry[6], entry[7], entry[8]]))
            velx.append(vx)
            vely.append(vy)
            velz.append(vz)
            x_pos.append(x)
            y_pos.append(y)
            z_pos.append(z)
        #Find max and min values
        x_max, x_min = max(x_pos), min(x_pos)
        y_max, y_min = max(y_pos), min(y_pos)
        z_max, z_min = max(z_pos), min(z_pos)
        #Find dimensions of new matrix, create zero matrix with those dimensions
        x_dim = x_max - x_min
        y_dim = y_max - y_min
        z_dim = z_max - z_min
        vx_field = np.zeros((x_dim+1, y_dim+1, z_dim+1)) #add 1 to avoid index errors
        vy_field = np.zeros((x_dim+1, y_dim+1, z_dim+1))
        vz_field = np.zeros((x_dim+1, y_dim+1, z_dim+1))
        #modify entries of x_pos, y_pos, and z_pos to start at 0 (rather than a negative #)
        if x_min < 0:
            x_pos = np.array(x_pos) + abs(x_min)
        else:
            x_pos = np.array(x_pos) - x_min
        if y_min < 0:
            y_pos = np.array(y_pos) + abs(y_min)
        else:
            y_pos = np.array(y_pos) - y_min
        if z_min < 0:
            z_pos = np.array(z_pos) + abs(z_min)
        else:
            z_pos = np.array(z_pos) - z_min
        #Add all velocities to 3D matrices - vx, vy, and vz matrices
        for i in range(len(x_pos)):
            x, y, z = x_pos[i], y_pos[i], z_pos[i]
            vx, vy, vz = velx[i], vely[i], velz[i]
            vx_field[x][y][z] = vx
            vy_field[x][y][z] = vy            
            vz_field[x][y][z] = vz
    return (vx_field,vy_field,vz_field)


create_velocity_field('section2veldata.csv')

vx_field_artery, vy_field_artery, vz_field_artery = create_velocity_field('section2veldata.csv')
#print(np.count_nonzero(vz_field_artery))



        

