#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 30 10:50:21 2017

@author: taylorsorenson, zhaoqiwang
"""

#Attempting to push to github


class Blood(object):
    ''' '''
    def __init__(self):
        self.dose_added = 0
        
    def set_dose(self,dose):
        self.dose_added += dose
        
    def get_dose(self):
        return self.dose_added


class Position(object):
    def __init__(self, dose):
        self.vx = 0
        self.vy = 0
        self.vz = 0
        self.velocity = (self.vx,self.vy,self.vz)
        self.x = 0
        self.y = 0
        self.z = 0
        self.position = (self.x, self.y, self.z)

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
    
def simulate_blood_flow(total_time, dt):
    pass
        
def make_pdf(blood):
    pass


def make_cdf():
    pass

def make_dvh():
    
        
        
    
        


