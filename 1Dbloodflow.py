#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 30 10:50:21 2017

@author: taylorsorenson, zhaoqiwang
"""

import numpy as np 
import math
import random

class Blood(object):
    ''' A Blood cell with an initial position (and velocity)
    '''
    def __init__(self, position, init_dose = 0):
        self.position = position
        self.dose_added = init_dose
        
    def get_dose(self):
        return self.dose_added
     
    def set_dose(self,dose):
        self.dose_added += dose
        
    def get_position(self):
        return self.position
    

class Position(object):
    '''A position within the blood vessel
    '''
    def __init__(self, x, y, z, dose):
        #assumes velocity and position are tuples in all three dimensions
        #in 1D case, vy and vz are 0 --> velocity = (vx, 0, 0)
        self.x = x
        self.y = y
        self.z = z
        self.position = (self.x, self.y, self.z)
        self.dose = dose #does this need to be here?
        
    def get_position(self):
        return self.position    
    
    def update_position(self, velocity):
        pass
        
class const_vector_field(object):
        #assumes the velocity everywhere is constant in one direction
        def __init__(self, x_dim, y_dim, z_dim, vx,  vy = 0, vz = 0): #change the order of these later
            '''Assumes x_dim is a scalar of the number of units in our matrix
            In one dimension, y_dim and z_dim should just be 1
            #TODO - Adjust this later to match up with the velocity fields  in vdx files 
            '''             
            self.vx = vx
            self.vy = vy
            self.vz = vz
            self.velocity = (self.vx, self.vy, self.vz)
            #create vector field, all of which have the same (and constant) velocity 
            self.field = np.zeros((z_dim,y_dim, x_dim)) + self.vx            
        
        def get_velocity(self):
            return self.velocity
        
        def get_v_at_position(self,pos):
            x = pos.get_position()[0]
            y = pos.get_position()[1]
            z = pos.get_position()[2]
            return self.field[x][y][z]
        
        def get_field(self):
            return self.field


field = const_vector_field(10,1,1,5.2)
print(field.get_field())
print(field.get_v_at_position((0,0,1)))



def make_blood(num_blood_cells):
    '''makes num_blood_cells objects and gives them all an initial position
    '''
    
        
def find_blood():
    '''finds which blood cells are in the path of the beam
    '''       
        
def add_constant_dose(blood, dose):
    '''
    add a constant dose of radiation to all blood within the radiation beam
    blood - a list of blood objects that are
    '''
    #time?
    pass 

def simulate_blood_flow(total_time, dt):
    pass
        
def make_pdf(blood):
    pass


def make_cdf():
    pass

def make_dvh():
    pass
        
        
    
        


