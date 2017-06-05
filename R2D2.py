"""
Created on Thu Jun  1 16:38:04 2017

@author: taylorsorenson, zhaoqiwang, 
"""

import numpy as np 
from matplotlib import pyplot as plt
import math
import random
from readDoses import *
import time


class Blood(object):
    ''' A Blood voxel with an position and dose
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
        x = int(position.x) #TODO - these may need to be getter functions later
        y = int(position.y) #ie y = position.get_y()
        z = int(position.z)
        dose = dose_matrix[x][y][z]
        self.dose_recive = dose
        
    def is_in_field(self):
        pass
    
    def __str__(self):
        return "Blood at position" + str(self.get_position()) + "w/ Dose of " + \
                str(self.get_dose())
                
        

#
#pos = Position(1,2,3)
#print(pos)
#blood_1 = Blood(pos, 4)
##print(blood_1)

class Position(object):
    '''A position within the blood vessel. 
    '''

    def __init__(self, x, y, z): 
        #assumes x, y, and z are int, representing a voxel within the vector field
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
    
    def __str__(self):
        pass

class Vector_field(object):
    '''TODO  - make this the superclass of const_vector_field. Make another
    subclass of vector_field called varying_field(?)
    '''
    '''Each position has an associated velocity in x,y, and z directions as 
    well as an associated dose'''
    def __init__(self, x_dim, y_dim, z_dim): #change the order of these later
        '''Assumes x_dim is a scalar of the number of units in our matrix
        In one dimension, y_dim and z_dim should just be 1
        #TODO - Adjust this later to match up with the velocity fields  in vdx files 
        '''             
        self.x_dim, self.y_dim, self.z_dim = x_dim, y_dim, z_dim
        #create three 3-D velocity matrices, each containing the velocity in one direction x,y,or z  
        
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
       
    
    def get_dimensions(self):
        return (self.x_dim, self.y_dim, self.z_dim)
    
    def get_x_dim(self):
        return self.x_dim
    
    def get_y_dim(self):
        return self.y_dim
    
    def get_z_dim(self):
        return self.z_dim
        
        
    def is_position_in_dose_field(self, position):
        '''Returns true if the input position is within the dose_field
        '''
        x = math.floor(position.get_x()) #math.floor makes the position an int
        y = math.floor(position.get_y()) #TODO - math.floor not really needed
        z = math.floor(position.get_z())
        return self.dose_matrix[x][y][z] != 0 #if dose is nonzero, must be in dose field
        
    def is_position_in_vector_field(self, position):
        '''Returns true if position is within the blood vessel
        '''
        x = position.get_x()
        y = position.get_y() 
        z = position.get_z()
        return (0 <= x < self.x_dim and 0 <= y < self.y_dim and 0 <= z < self.z_dim) 
        #less than or equal to? <=     

    def __str__(self):
            return "Vector_field w/ dim: " + str((self.x_dim, self.y_dim, \
                                                         self.z_dim))

class Const_vector_field(Vector_field):
        '''Each position has an associated velocity in x,y, and z directions as 
            well as an associated dose'''
            #NOTE - vy does not have a default value of 0 in 2d version
        def __init__(self, x_dim, y_dim, z_dim, vx, vy, vz): #change the order of these later
            '''Assumes x_dim is a scalar of the number of units in our matrix
            In one dimension, y_dim and z_dim should just be 1
            #TODO - Adjust this later to match up with the velocity fields  in vdx files 
            '''     
            Vector_field.__init__(self,x_dim, y_dim, z_dim)
            self.vx, self.vy, self.vz = vx, vy, vz
            self.velocity = (self.vx, self.vy, self.vz)
            #create three 3-D velocity matrices, each containing the velocity in one direction x,y,or z
            self.vx_field = np.zeros((x_dim,y_dim, z_dim)) + self.vx   
            self.vy_field = np.zeros((x_dim,y_dim, z_dim)) + self.vy  
            self.vz_field = np.zeros((x_dim,y_dim, z_dim)) + self.vz  
        
        def get_vx_at_position(self,x,y,z):
            return self.vx_field[int(x)][int(y)][int(z)]       
        
        def get_vy_at_position(self,x,y,z):
            return self.vy_field[int(x)][int(y)][int(z)] 
        
        def get_vz_at_position(self,x,y,z):
            return self.vz_field[int(x)][int(y)][int(z)]   

        def get_velocity_fields(self):
            return [self.vx_field, self.vy_field, self.vz_field]                         
        

def make_blood(num_blood_cells,x_min = 0, y_min = 0, z_min = 0, \
               x_max=1,y_max=50,z_max=80):
    '''makes num_blood_cells objects and gives them all an initial position
    '''
    bloods =[]
    x = np.random.uniform(x_min,x_max,num_blood_cells)#changed from randit to uniform
    y = np.random.uniform(y_min,y_max,num_blood_cells)
    z = np.random.uniform(z_min,z_max,num_blood_cells)
    for i in range(num_blood_cells):
        position = Position(x[i],y[i],z[i])
        blood = Blood(position)
        bloods += [blood]
        
    return bloods

def plot_positions(blood, t):
    '''Blood is a list of all the blood voxels objects, t is a scalar representing
    at what time the blood is plotted
    '''
    pass   
  
def add_dose_for_allblood(all_bloods,dose_matrix, vector_field):
    '''
    add a constant dose of radiation to all blood within the radiation beam
    blood - a list of blood voxel objects
    '''

    for i in all_bloods:
        if vector_field.is_position_in_dose_field(i.get_position()):
            i.current_dose_level(dose_matrix)
            i.add_dose()

def generate_new_blood(lefted_bloods,velocity_field):
    '''generate new blood as to simulate the blood flow into the region from y = 0, based 
        on the velocity and the cross section area in yz-plane
    '''
    num_blood = len(lefted_bloods)
    vy_field = velocity_field.get_velocity_fields()[1]
    mean_vy = vy_field.
    
    bloods = make_blood(num_blood,y_max = 0)
    
    return bloods
    
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
        if vector_field.is_position_in_vector_field(new_position):
            i.find_new_position(new_position)
        else:
            #create a new blood object in a random y position at z = 0 (and x = 0 for 2d)
            random_y = random.random()*vector_field.get_y()
            new_blood = Blood(Position(0, random_y, 0))
            all_bloods.append(new_blood)
        

 

def blood_flow_with_beam(all_bloods,vector_field,dose_matrix,total_time, dt):
    """simulate the flow of blood while the radiation beam IS ON, assumes
    a dose_matrix is a two dimension matrix of the dose at each point in a 
    2d space
    """
    t = 0;
    while t <= total_time:
        add_dose_for_allblood(all_bloods,dose_matrix, vector_field) 
        bloods_flow(all_bloods,vector_field,dt)
        t = t + dt
 
    return all_bloods 

def blood_flow_no_beam(all_bloods, vector_field, time_gap, dt):
    '''simulate blood flow while radiation beam IS OFF. 
    time_gap is the amount of time the beam is off (seconds), dt is length
    of each time step.
    '''
    t = 0
    while t<= time_gap:
        #only update the position of each blood voxel
        bloods_flow(all_bloods, vector_field, dt)
        t += dt
    return all_bloods

def simulate_blood_flow(dose_fields, times, time_gaps, dt, blood_density= 1):
    '''generate blood objects and velocity vector fields, then runs the simulation
    NOTE - assumes the starting point is when the first field turns on
    dose_fields: a list of dose_field matrices
    times: a list of the total time each dose matrix is 'on'
    time_gaps: the time between doses, while blood is still moving
    '''
    vector_field = Const_vector_field(1,120,120,0,1,3) #TODO = change this vy
    #make a matrix of random values between 0 and 1  
    num_bloods = blood_density * vector_field.size()#total number of the blood
    all_bloods = make_blood(num_bloods)
    for i in range(len(times)):
        vector_field.set_dose_matrix(dose_fields[i])
        blood_after_dose = blood_flow_with_beam(all_bloods, vector_field, dose_fields[i], times[i], dt)
        try:
            all_bloods = blood_flow_no_beam(blood_after_dose, vector_field, time_gaps[i], dt)
        except IndexError: #because len(time_gaps) < len(times); 2 compared to 3, for ex.
            pass
    return blood_after_dose #Note - don't return all_bloods


        
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
    


def test_pdf(dose_fields, times, time_gaps, dt, n_bloods = 100):
    '''tests makepdf and plots pdf
    '''
    blood_cells = simulate_blood_flow(dose_fields, times, time_gaps,dt, num_bloods = n_bloods)
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
    
 
def test_cdf(dose_fields, times, time_gaps, dt, n_bloods = 100):
    '''tests make_cdf and plots cdf
    '''
    blood_cells = simulate_blood_flow(dose_fields, times, time_gaps,dt, num_bloods = n_bloods)
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
    
    
        
def test_dvh(dose_fields, times, time_gaps, dt, n_bloods = 100):      
    blood_cells = simulate_blood_flow(dose_fields, times, time_gaps,dt, num_bloods = n_bloods)
    bin_centers, dvh = make_dvh(blood_cells)
    plt.figure()
    plt.title("Dose-Volume Histogram\n # of Blood Voxels: " + str(n_bloods))
    plt.xlabel("Dose (Gray)")
    plt.ylabel("Volume (%)")
    plt.plot(bin_centers, dvh)
    plt.grid(True)
    

if __name__ == '__main__': 
    pass
#    num_blood_cells = 100
#    test_pdf(doses, time_on, time_off, .1, n_bloods = 1000)
#    test_cdf(doses, time_on, time_off, .1, n_bloods = 1000)
    for n in [10,100,1000, 10000]:
        test_dvh(doses, time_on, time_off, .1, n_bloods = n) 
    

        
