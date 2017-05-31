#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 30 10:50:21 2017

@author: taylorsorenson, zhaoqiwang
"""

import numpy as np 
import math
import random
import readDoses

#Attempting to push to github



class Blood(object):
    ''' A Blood cell with an initial position (and velocity)
    '''
    def __init__(self, position, init_dose = 0):
        self.position = position
        self.dose_added = init_dose
        self.getting_dose = False;

        
    def get_dose(self):
        return self.dose_added
    
    def add_dose(self,dose):
        self.dose_added += dose

    def get_position(self):
        return self.position
    
    def find_new_position(self, position):
        self.position = position
        
    
#    dose = dose_matrix[1][2][4] 
        
    

class Position(object):
    '''A position within the blood vessel. 
    '''
    def __init__(self, x, y, z):
        #assumes velocity and position are tuples in all three dimensions
        #in 1D case, vy and vz are 0 --> velocity = (vx, 0, 0)

    def __init__(self, x, y, z): #assumes x, y, and z are floats
        self.x = x
        self.y = y
        self.z = z
        self.position = (self.x, self.y, self.z)
        
    def get_position(self):
        return self.position 
    
    def get_x():
        return self.x
    
    def get_y():
        return self.y
    
    def get_z():
        return self.z
    

    def get_new_position(self, vx,vy,vz, dt):#change to units rather than actual distance
        #dt is the change in time dt
        old_x, old_y, old_z = self.get_x(), self.get_y(), self.get_z()
        dx = vx * dt
        dy = vy * dt
        dz = vz * dt
        
        new_x, new_y, new_z = (old_x + dx), (old_y + dy), (old_z + dz)
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
            y = math.floor(position.get_y())
            z = math.floor(position.get_z())
            return self.dose_matrix[x][y][z] != 0

#field = const_vector_field(10,1,1,5.2)
#print(field.get_field())
#print(field.get_v_at_position((0,0,1)))



def make_blood(num_blood_cells):
    '''makes num_blood_cells objects and gives them all an initial position
    '''
    
        
def find_blood():
    '''finds which blood cells are in the path of the beam
    ''' 
    if position.dose(blood.position) !=  0:
        blood.getting_dose = True
        blood.dose = position.dose(blood.position)
        
        
               
        
def add_constant_dose(all_bloods):
    '''
    add a constant dose of radiation to all blood within the radiation beam
    blood - a list of blood objects that are
    '''

    for i in all_bloods:
        if i.getting_dose == True;
            i.setdose(i.dose)
    
def bloodflow(blood,position,dt):
    vx = positon.vx(blood.position)
    vy = positon.vy(blood.position)
    vz = positon.vz(blood.position)
    blood.x = blood.x + dt * vx
    blood.y = blood.y + dt * vy
    blood.z = blood.z + dt * vz

def simulate_blood_flow(all_bloods,Positions,total_time, dt):
    t = 0;
    while t <= total_time:
        add_constant_dose(all_bloods,dose) 
        blood_flow(all_bloods,)
        t = t + dt
    #time?
    pass 
        
def make_pdf(blood):
    pass


def make_cdf():
    pass

def make_dvh():
    pass
        
        
    


#field = const_vector_field(10,1,1,5.2)
#print(field.get_field())
#print(field.get_v_at_position((0,0,1)))




