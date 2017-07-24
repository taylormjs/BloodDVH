#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 30 10:50:21 2017

@author: taylorsorenson, zhaoqiwang
"""

import numpy as np 
from matplotlib import pyplot as plt
import math
import random
import readDoses
random.seed(1)


#print('f1' in mat)

class Blood(object):
    ''' A Blood cell with an initial position (and velocity)
    '''
    def __init__(self, position, init_dose = 0):
        self.position = position
        self.dose_added = init_dose
        self.getting_dose = False;
        self.dose_recive = 0

        
    def get_dose(self):
        return self.dose_added
    
    def add_dose(self):
        self.dose_added += self.dose_recive

    def get_position(self):
        return self.position
    
    def find_new_position(self, position):
        self.position = position
        
    def current_dose_level(self,dose_matrix):
        '''return dose that each unit of blood gets in a dose_matrix
        '''
        self.dose_recive = 0
        position = self.get_position()
        x = position.x
        y = position.y
        z = position.z
        dose = dose_matrix[x][y][z]
        self.dose_recive = dose
#    dose = dose_matrix[1][2][4] 
#    def __str__(self):
#        print("Blood at position", self.get_position(), "w/ Dose of ", self.get_dose()
#       
#pos = Position(1,2,3)
#blood_1 = Blood(pos, 4)
#print(blood_1)

class Position(object):
    '''A position within the blood vessel. 
    '''

    def __init__(self, x, y, z): 
        #assumes x, y, and z are int, representing a box within the vector field
        #NOT REPRESENTING an exact location in the x, y, and z
        self.x = x
        self.y = y
        self.z = z
        self.position = (self.x, self.y, self.z)
        
    def get_position(self):
        return self.position 
    
    def get_x(self):
        return self.x
    def get_y(self):
        return self.y
    def get_z(self):
        return self.z
    def get_new_position(self, vx,vy,vz, dt):#change to units rather than actual distance
        #dt is the change in time dt
        old_x, old_y, old_z = self.get_x(), self.get_y(), self.get_z()
        dx = vx * dt
        dy = vy * dt
        dz = vz * dt
        #make the new positions units within the vector field, hence use math.floor
        new_x, new_y, new_z = math.floor(old_x + dx), math.floor(old_y + dy), \
                                        math.floor(old_z + dz)
        return Position (new_x, new_y, new_z)

        
class const_vector_field(object):
        '''Each position has an associated velocity in x,y, and z directions as 
            well as an associated dose'''
        def __init__(self, x_dim, y_dim, z_dim, vx,  vy = 0, vz = 0): #change the order of these later
            '''Assumes x_dim is a scalar of the number of units in our matrix
            In one dimension, y_dim and z_dim should just be 1
            #TODO - Adjust this later to match up with the velocity fields  in vdx files 
            '''             
            self.vx, self.vy, self.vz = vx, vy, vz
            self.x_dim, self.y_dim, self.zdim = x_dim, y_dim, z_dim
            self.velocity = (self.vx, self.vy, self.vz)
            #create three 3-D velocity matrices, each containing the velocity in one direction x,y,or z
            self.vx_field = np.zeros((z_dim,y_dim, x_dim)) + self.vx   
            self.vy_field = np.zeros((z_dim,y_dim, x_dim)) + self.vy  
            self.vz_field = np.zeros((z_dim,y_dim, x_dim)) + self.vz   
        
        def set_dose_matrix(self, dose_matrix): #maybe take out the dim parameters later
            '''assume dose matrix a numpy matrix of the same dimensions as velocity field''' 
            self.dose_matrix = dose_matrix
            
        def fit_dose_matrix(self, exp_dose_matrix): #exp stands for expanded
            '''
            #assumes dose matrix IS NOT the same dimensions as velocity field
            #simply takes in a matrix and truncates the matrix to be the same dimensions
            #as the velocity field
            '''
            pass #do later
       
        def get_vx_at_position(self,x,y,z):
            return self.vx_field[x][y][z]
        
        
        def get_vy_at_position(self,x,y,z):
            return self.vy_field[x][y][z]
        
        def get_vz_at_position(self,x,y,z):
            return self.vz_field[x][y][z]
        
        def get_velocity_fields(self):
            return [self.vx_field, self.vy_field, self.vz_field]
        
        def is_position_in_dose_field(self, position):
            '''find which x, y, and z coordinate the position is at within the field'''
            x = math.floor(position.get_x()) #math.floor makes the position an int
            y = math.floor(position.get_y()) #TODO - math.floor not really needed
            z = math.floor(position.get_z())
            return self.dose_matrix[x][y][z] != 0 #if dose is nonzero, must be in dose field
        
        def is_position_in_vector_field(self, position):
            x = position.get_x()
            y = position.get_y() 
            z = position.get_z()
            return (0 <= x <= self.x_dim and 0 <= y <= self.y_dim and 0 <= z <= self.z_dim)

#field = const_vector_field(10,1,1,5.2)
#print(field.get_field())
#print(field.get_v_at_position((0,0,1)))



def make_blood(num_blood_cells,x_max=100,y_max=100,z_max=100):
    '''makes num_blood_cells objects and gives them all an initial position
    '''
    bloods =[]
    x = np.random.randint(0,x_max,num_blood_cells)
    y = np.random.randint(0,y_max,num_blood_cells)
    z = np.random.randint(0,z_max,num_blood_cells)
    for i in range(num_blood_cells):
        position = Position(x[i],y[i],z[i])
        blood = Blood(position)
        bloods += [blood]
        
    return bloods
        
  
def add_dose_for_allblood(all_bloods,dose_matrix):
    '''
    add a constant dose of radiation to all blood within the radiation beam
    blood - a list of blood objects that are
    '''

    for i in all_bloods:
        i.current_dose_level(dose_matrix)
        i.add_dose()
    
def bloods_flow(all_bloods, vector_field,dt):
    for i in all_bloods:
        position = i.get_position()
        x = position.x
        y = position.y
        z = position.z
        vx = vector_field.get_vx_at_position(x,y,z)
        vy = vector_field.get_vy_at_position(x,y,z)
        vz = vector_field.get_vz_at_position(x,y,z)
        new_position = position.get_new_position(vx,vy,vz, dt)
        i.find_new_position(new_position)
    

def simulate_blood_flow(all_bloods,vector_field,dose_matrix,total_time, dt):
    "simulate the blood flow in some given dose matrix over some time"
    t = 0;
    while t <= total_time:
        add_dose_for_allblood(all_bloods,dose_matrix) 
        bloods_flow(all_bloods,vector_field,dt)
        t = t + dt
 
    return all_bloods 


def test_blood_flow(total_t,dt):
    '''generate blood object and vector fields and run the simulations'''
    all_bloods = make_blood(10000,x_max=100,y_max=100,z_max=100)
    vector_field = const_vector_field(120,120,120,1)
    #make a matrix of random values between 0 and 1
    dose_matrix =  np.random.random((120,120, 120))
    
    new_all_bloods = simulate_blood_flow(all_bloods,vector_field,dose_matrix,total_t, dt)
    return new_all_bloods


        
def make_pdf(blood_cells):
    '''make a probability density function which will graph blood dose vs. 
    volume of blood (which will be fraction of total blood for 1D)
    Assumes blood_cells is a list of blood objects, all of which have varying
    doses
    '''
    #find the doses of all the cells, append to doses list
    total = len(blood_cells)
    doses = []
    for cell in blood_cells:
        doses.append(cell.get_dose())
    hist, bins = np.histogram(doses, bins='auto', density=False) #normed = True instead?
    bin_centers = (bins[1:]+bins[:-1])*0.5
    return (bin_centers,hist) #returned this way to make easier to plot later
    
def plot_pdf(blood_cells):
    '''plots the data from make_pdf'''
    #plot these doses on a histogram
    bin_centers, hist = make_pdf(blood_cells)
#    bin_centers, hist = make_pdf(blood_cells)
    plt.figure()
    plt.title("Probabilty Density Function")
    plt.xlabel("Dose (Gray)")
    plt.ylabel("Frequency")
    plt.plot(bin_centers, hist)
    plt.grid(True)
    plt.show()
    


def test_pdf(total_t, dt):
    '''tests makepdf and plots pdf
    '''
    blood_cells = test_blood_flow(total_t,dt)
    plot_pdf(blood_cells)
         


def make_cdf(blood_cells):
    '''make the cumulative density function'''
    #create new numpy array to plot
    bin_centers, hist = make_pdf(blood_cells)
    cumsum = np.cumsum(hist)/len(blood_cells)
    plt.figure()
    plt.plot(bin_centers, cumsum)
    plt.title("Cumulative Density Function")
    plt.xlabel("Dose (Gray)")
    plt.ylabel("% of Blood Cells")
    plt.grid(True)
    
 
def test_cdf(total_t, dt):
    '''tests make_cdf and plots cdf
    '''
    blood_cells = test_blood_flow(total_t, dt)
    make_cdf(blood_cells)


def make_dvh(blood_cells): 
    '''
    Note - this is independent from the pdf and cdf functions now to avoid
    going through all the blood cells multiple times
    '''
    total = len(blood_cells)
    doses = []
    for cell in blood_cells:
        doses.append(cell.get_dose())
    hist, bins = np.histogram(doses, bins='auto', density=False) #normed = True instead?
    bin_centers = (bins[1:]+bins[:-1])*0.5
    dvh = np.cumsum(hist)/len(blood_cells)
    for i in range(len(dvh)):
        dvh[i] = 1 - dvh[i]
    return (bin_centers,dvh)
    
    
        
def test_dvh(total_t, dt):      
    blood_cells = test_blood_flow(total_t, dt)
    bin_centers, dvh = make_dvh(blood_cells)
    plt.figure()
    plt.title("Dose-Volume Histogram")
    plt.xlabel("Dose (Gray)")
    plt.ylabel("Volume (%)")
    plt.plot(bin_centers, dvh)
    plt.grid(True)
    

if __name__ == '__main__': 
#    num_blood_cells = 100
#    test_pdf(10,.1)
#    test_cdf(10,.1)
    test_dvh(10,.1) 
    

        






