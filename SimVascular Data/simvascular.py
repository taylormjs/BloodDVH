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


vx = open('vxdata.csv', 'r')
vy = open('vydata.csv', 'r')
vz = open('vzdata.csv', 'r')

#for line in vx:
#    print(line)

def create_velocity_field(filename):
    pass