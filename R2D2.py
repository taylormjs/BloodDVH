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
from dose_field_plot import *


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
        
    def current_dose_level(self,dose_matrix,dt):
        '''return dose that each unit of blood gets in a dose_matrix
        '''
        self.dose_recive = 0
        position = self.get_position()
        x = int(position.x) #TODO - these may need to be getter functions later
        y = int(position.y) #ie y = position.get_y()
        z = int(position.z)
        dose = dose_matrix[x][y][z] * dt
        self.dose_recive = dose
        
    def is_in_field(self):
        pass
    
    def __str__(self):
        return "Blood at " + str(self.get_position()) + " w/ Dose of " + \
                str(self.get_dose())
                

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
        new_x, new_y, new_z = (old_x + dx), (old_y + dy), (old_z + dz)
        return Position (new_x, new_y, new_z)
    
    def __str__(self):
        return "Position: " + str((self.x, self.y, self.z))
        #'(' + str(self.x) + str(self.y) + str(self.z) + ')'

class Vector_field(object):
    '''TODO  - make this the superclass of const_vector_field. Make another
    subclass of vector_field called varying_field(?)
    '''
    '''Each position has an associated velocity in x,y, and z directions as 
    well as an associated dose'''
    def __init__(self, vx_field, vy_field, vz_field): #inputs vx,vy,vz from vtk or csv files
        '''vx_field, vy_field, vz_field are ndarrays
        #TODO - Adjust this later to match up with the velocity fields  in vdx files 
        ''' 
        self.vx_field = vx_field
        self.vy_field = vy_field
        self.vz_field = vz_field
        self.v_square = self.get_v_square()        
        (self.x_dim, self.y_dim, self.z_dim) = vx_field.shape
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
        pass #TODO - later
       
    def get_vx_at_position(self,x,y,z):
        return self.vx_field[int(x)][int(y)][int(z)]       
        
    def get_vy_at_position(self,x,y,z):
        return self.vy_field[int(x)][int(y)][int(z)] 
        
    def get_vz_at_position(self,x,y,z):
        return self.vz_field[int(x)][int(y)][int(z)]   
        
    def get_vx_field(self):
        return self.vx_field                       

    def get_vy_field(self):
        return self.vy_field  
        
    def get_vz_field(self):
        return self.vz_field  
        
    def get_v_square(self):
        return self.vx_field**2 + self.vy_field**2 + self.vz_field**2
        
    def get_velocity_fields(self):
        return [self.vx_field, self.vy_field, self.vz_field]             
    
    def get_dimensions(self):
        return (self.x_dim, self.y_dim, self.z_dim)
    
    def get_x_dim(self):
        return self.x_dim
    
    def get_y_dim(self):
        return self.y_dim
    
    def get_z_dim(self):
        return self.z_dim
    
    def get_size(self):
        return self.vx_field.size
               
    def is_position_in_dose_field(self, position):
        '''Returns true if the input position is within the dose_field
        '''
        x = int(position.get_x()) 
        y = int(position.get_y()) #Changed math.floor to int to save time
        z = int(position.get_z())
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
            #NOTE - vy dose not have a default value of 0 in 2d version
        def __init__(self, x_dim, y_dim, z_dim, vx, vy, vz): #change the order of these later
            '''Assumes x_dim is a scalar of the number of units in our matrix
            In one dimension, y_dim and z_dim should just be 1
            #TODO - Adjust this later to match up with the velocity fields  in vdx files 
            ''' 
            self.vx, self.vy, self.vz = vx, vy, vz
            velocity = (self.vx, self.vy, self.vz)
            vx_field = np.zeros((x_dim,y_dim, z_dim)) + self.vx   
            vy_field = np.zeros((x_dim,y_dim, z_dim)) + self.vy  
            vz_field = np.zeros((x_dim,y_dim, z_dim)) + self.vz  
            Vector_field.__init__(self,vx_field, vy_field, vz_field)
  
            
        

def make_blood(num_blood_cells,x_min = 0, y_min = 0, z_min = 0, \
               x_max=120,y_max=120,z_max=120):
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

  
  
def add_dose_for_allblood(all_bloods,dose_matrix, vector_field,dt):
    '''
    add a constant dose of radiation to all blood within the radiation beam
    blood - a list of blood voxel objects
    '''
    vector_field.set_dose_matrix(dose_matrix)
    for i in all_bloods:
        try:
            if vector_field.is_position_in_dose_field(i.get_position()):
                i.current_dose_level(dose_matrix,dt)
                i.add_dose()
        except IndexError:
            pass

def gen_new_blood(out_blood_per_t,vector_field, in_boundary, axis):
    '''generate new blood as to simulate the blood flow into the region from y = 0, based 
    on the velocity and the cross section area in yz-plane
    out_bloods (int)  - number of bloods leaving each time step, will vary
    total_velocity (float) - the sum of the magnitudes of all velocities in a specified
    area
    in_boundary (tuple of ranges y_lo, y_hi,z_lo, z_hi))- the boundary 
    at which the new blood cells can enter
    axis (str) - can be either 'y' or 'z', defining the axis along which the
    blood cells will be added (ie - if the blood flows in, crossing the z axis,
    axis = 'z')
    '''    
    #TODO - if this is taking too long later, separate into two functions
    
    #Initiliaze num_bloods, boundary limits, and fields
    #num_bloods = len(out_bloods)
    lo = in_boundary[0]
    hi = in_boundary[1]
    v_square = vector_field.get_v_square()
    #Find the velocities at the units adjacent to the boundary
    if axis == 'y':
        mini_field = v_square[0,lo:hi+1, 0]        
    elif axis == 'z':
        mini_field = v_square[0, 0, lo:hi+1]
    elif axis == 'x': # TODO  - note that x won't work for 2D case as x = 0 
        mini_field = v_square[lo:hi+1, 0, 0] #TODO - 
    else:
        raise ValueError("must input 'x', 'y', or 'z' as the axis. 'y' is default")
    #find probability matrix to find indices at which the blood will enter
    total = sum(mini_field)
    prob_field =  mini_field / total
    #generate a 1D list of number_bloods leaving using np.random.choice()
    #with probabilities defined, this will be the positions of the new blood
    #units
    choose_ind = np.random.choice(len(prob_field),out_blood_per_t, p = prob_field) #Note - a list of indices
    ran = np.random.random(out_blood_per_t)                          
    gen_bloods = []
    for i in range(len(choose_ind)):
        if axis == 'y':
            y = choose_ind[i]+ ran[i] + lo           
            gen_blood = Blood(Position(0,y,0))
        elif axis == 'z':
            z = choose_ind[i]+ ran[i] + lo  
            gen_blood = Blood(Position(0,0,z))
        elif axis == 'x':
            x = choose_ind[i]+ ran[i] + lo 
            gen_blood = Blood(Position(x,0,0))
        gen_bloods.append(gen_blood)
    
    return gen_bloods

def gen_new_blood_3d(out_blood_per_t,vector_field, in_boundary,direction='z'):
    '''generate new blood as to simulate the blood flow into the region from y = 0, z = 0 based 
    on the velocity and the cross section area in yz-plane
    out_bloods (int)  - number of bloods leaving each time step, will vary
    vector_field(vector_field object) the vector field object contains the velocity of f
    in_boundary (a list of tuple of ranges[(x_lo,x_hi),(y_lo, y_hi),(z_lo, z_hi)])- the boundary 
    at which the new blood cells can enter
    axis (str) - can be either 'y' or 'z', defining the axis along which the
    blood cells will be added (ie - if the blood flows in, crossing the z axis,
    axis = 'z')
    '''    
    #TODO - if this is taking too long later, separate into two functions
    
    #Initiliaze num_bloods, boundary limits, and fields
    #num_bloods = len(out_bloods)
    x_boundary = in_boundary[0]
    y_boundary = in_boundary[1]
    z_boundary = in_boundary[2]
    x_lo = x_boundary[0]
    x_hi = x_boundary[1]
    x_len = x_hi -x_lo
    y_lo = y_boundary[0]
    y_hi = y_boundary[1]
    y_len = y_hi -y_lo
    z_lo = z_boundary[0]
    z_hi = z_boundary[1]
    z_len = z_hi -z_lo
    v_square = vector_field.v_square
    gen_bloods = []
    if direction == 'z':
        mini_field = v_square[x_lo:x_hi+1,y_lo:y_hi+1, 0]   
        total = mini_field.sum()
        prob_field =  mini_field / total
    #generate a 1D list of number_bloods leaving using np.random.choice()
    #with probabilities defined, this will be the positions of the new blood
    #units
        choose_ind = np.random.choice(prob_field.size,out_blood_per_t, \
                                      p = prob_field.flatten()) #flatten the matrix to 1-d                           

        for i in range(out_blood_per_t):
            x = choose_ind[i]//y_len #reconsturct the matrix
            y = choose_ind[i] % y_len
            x = x + np.random.random() + x_lo 
            y = y + np.random.random() + y_lo                        
            gen_blood = Blood(Position(x,y,0))                  
            gen_bloods.append(gen_blood)
    
    elif direction == 'y':
        mini_field = v_square[x_lo:x_hi+1 ,0, z_lo:z_hi+1]   
        total = mini_field.sum()
        prob_field =  mini_field / total
        choose_ind = np.random.choice(prob_field.size,out_blood_per_t, \
                                      p = prob_field.flatten())                         

        for i in range(out_blood_per_t):
            x = choose_ind[i]//z_len #reconsturct the matrix
            z = choose_ind[i] % z_len
            x = x + np.random.random() + x_lo 
            z = z + np.random.random() + z_lo                        
            gen_blood = Blood(Position(x,0,z))                  
            gen_bloods.append(gen_blood)
    elif direction == 'x':
        mini_field = v_square[0, y_lo:y_hi+1, z_lo:z_hi+1]   
        total = mini_field.sum()
        prob_field =  mini_field / total
        choose_ind = np.random.choice(prob_field.size,out_blood_per_t, \
                                      p = prob_field.flatten()) 
        #flatten the matrix to 1d
        for i in range(out_blood_per_t):
            y = choose_ind[i]//z_len #reconsturct the matrix
            z = choose_ind[i] % z_len
            y = y + np.random.random() + y_lo 
            z = z + np.random.random() + z_lo                         
            gen_blood = Blood(Position(0,y,z))                 
            gen_bloods.append(gen_blood)      
    
    return gen_bloods

    
    
    
def bloods_flow(in_bloods, out_bloods, vector_field,dt, in_boundary = [(0,20),(10,60),(0,20)], direction = 'z'):
    '''Updates the positions of all the bloods still within the vector field
    in one time step with length dt (a float)
    in_bloods - those still within defined vector field
    out_bloods - those that have left vector field
    vector_field - vector_field object
    in_boundary - the boundary the blood is flowing into, a tuple (ylo, yhi) if axis == 'y'
    '''

    #Move all bloods within the field
    out_blood_count = 0
    for i in in_bloods:
        position = i.get_position()
        x = position.x
        y = position.y
        z = position.z
        vx = vector_field.get_vx_at_position(x,y,z)
        vy = vector_field.get_vy_at_position(x,y,z)
        vz = vector_field.get_vz_at_position(x,y,z)
        new_position = position.get_new_position(vx,vy,vz, dt)
        #if new position still in field, update position
        if vector_field.is_position_in_vector_field(new_position):
            i.find_new_position(new_position)
        else:
            #move that blood from one list to another (all_bloods -> leaving_blods)
            out_bloods.append(i)
            in_bloods.remove(i)
            out_blood_count += 1
    #Generate new bloods, one blood unit for each blood unit out
    
    in_bloods += gen_new_blood_3d(out_blood_count,vector_field, in_boundary, direction)
#    print('the # bloods flows out' , out_blood_count)
# 

def blood_flow_with_beam(in_bloods, out_bloods, vector_field, \
                         dose_matrix,total_time, dt, in_boundary ):
    """simulate the flow of blood while the radiation beam IS ON, assumes
    a dose_matrix is a two dimension matrix of the dose at each point in a 
    2d space
    """
    t = 0;
    while t <= total_time:
        dose_per_time = dose_matrix / total_time #dose_per_time a dose matrix per unit time
        add_dose_for_allblood(in_bloods,dose_per_time, vector_field,dt) 
        bloods_flow(in_bloods, out_bloods, vector_field,dt,in_boundary) #TODO - add inboundary as argument later?
        t = t + dt
        
    return in_bloods

def blood_flow_no_beam(in_bloods, out_bloods, vector_field, \
                       time_gap, dt, in_boundary):
    '''simulate blood flow while radiation beam IS OFF. 
    time_gap is the amount of time the beam is off (seconds), dt is length
    of each time step.
    '''
    t = 0
    while t<= time_gap:
        #only update the position of each blood voxel
        bloods_flow(in_bloods, out_bloods, vector_field, dt,in_boundary)
        t += dt
    return in_bloods


def simulate_blood_flow(dose, vector_field, dt, blood_density= 1, \
                        in_boundary = [(0,20),(10,60),(0,20)]):
    '''generate blood objects and velocity vector fields, then runs the simulation
    NOTE - assumes the starting point is when the first field turns on
    dose_fields: a list of dose_field matrices
    times: a list of the total time each dose matrix is 'on'
    time_gaps: the time between doses, while blood is still moving
    '''
    #initialize vector field, make # of blood voxels based on blood density
    num_bloods = int(blood_density * vector_field.get_size()) #total number of the blood
    (x_dim , y_dim ,z_dim) = vector_field.get_dimensions()
    dose_fields = dose.get_dose_field()
    times = dose.get_dose_time()
    time_gaps = dose.get_dose_time_gap()
    in_bloods = make_blood(num_bloods,x_max =x_dim,y_max=y_dim,z_max=z_dim)
    out_bloods = []
    plot_bloods_3d(in_bloods, c = 'r', m = '^')# initial position of the bloods 
    for i in range(len(times)):
        #Add dose
        in_bloods_after_dose = blood_flow_with_beam(in_bloods, out_bloods, vector_field, dose_fields[i], \
                                                    times[i], dt, in_boundary)
        try:
            #then circulate blood without dose
            in_bloods = blood_flow_no_beam(in_bloods_after_dose, out_bloods, vector_field, time_gaps[i],\
                                           dt,in_boundary)
        except IndexError: #because len(time_gaps) < len(times); 2 compared to 3, for ex.
            pass
    #bloods still in vector fields and those that have left should be returned in the 
    #same list
    total_time = sum(times) + sum(time_gaps)
    num_out_per_time = len(out_bloods)/total_time
    print("# of bloods OUT of velocity field: " ,len(out_bloods)) 
    print("# of bloods OUT per time", num_out_per_time)
    print("# of bloods IN velocity field: " ,len(in_bloods)) 
    print("Original # of bloods", num_bloods)
    print("# of bloods generated", (len(out_bloods) + len(in_bloods) - num_bloods))
    plot_bloods_3d(in_bloods, c = 'b',m = 'o')# fianl postion of the bloods 

    return in_bloods + out_bloods

        

def make_pdf(bloods):
    '''make a probability density function which will graph blood dose vs. 
    volume of blood (which will be fraction of total blood for 1D)
    Assumes blood_cells is a list of blood objects, all of which have varying
    doses
    '''
    #find the doses of all the cells, append to doses list
    total = len(bloods)
    doses = []
    for cell in bloods:
        doses.append(cell.get_dose())
    hist, bins = np.histogram(doses, bins='fd') #normed = True instead?
    bin_centers = (bins[1:]+bins[:-1])*0.5
    return (bin_centers,hist) #returned this way to make easier to plot later
    
def plot_pdf(bloods):
    '''plots the data from make_pdf'''
    #plot these doses on a histogram
    bin_centers, hist = make_pdf(bloods)
#    bin_centers, hist = make_pdf(blood_cells)
    plt.figure()
    plt.title("Probabilty Density Function")
    plt.xlabel("Dose (Gray)")
    plt.ylabel("Frequency")
    plt.plot(bin_centers, hist)
    plt.grid(True)
    plt.show()


def make_dvh(bloods): 
    '''
    Note - this is independent from the pdf and cdf functions now to avoid
    going through all the blood cells multiple times
    '''
    total = len(bloods)
    print("plotting ", total, " cells")
    doses = []
    for cell in bloods:
        doses.append(cell.get_dose())
    hist, bins = np.histogram(doses, bins='fd') #normed or density = True instead?
    bin_centers = (bins[1:]+bins[:-1])*0.5
    dvh = np.cumsum(hist)
    return (bin_centers,dvh)
    
           
#def plot_dvh(dose_fields, times, time_gaps, dt, blood_d = 1):      
#    blood_cells = simulate_blood_flow(dose_fields, times, time_gaps, \
#                                      dt, blood_density = blood_d)
#    bin_centers, dvh = make_dvh(blood_cells)
#    plt.figure()
#    num_bloods = len(blood_cells)
#    plt.title("Dose-Volume Histogram\n # of Blood Voxels: " + str(num_bloods) + \
#              "\nBlood Density: " + str(blood_d))            
#    plt.xlabel("Dose (Gray)")
#    plt.ylabel("Volume (%)")
#    plt.ylim(0,100)
#    print(num_bloods)
#    plt.plot(bin_centers, (num_bloods - dvh)/num_bloods*100, c='green')
#    plt.grid(True)

def plot_dvh(bloods,blood_density):
    bin_centers, dvh = make_dvh(bloods)
    plt.figure()
    num_bloods = len(bloods)
    plt.title("Dose-Volume Histogram\n # of Blood Voxels: " + str(num_bloods) + \
              "\nBlood Density: " + str(blood_density))            
    plt.xlabel("Dose (Gray)")
    plt.ylabel("Volume (%)")
    plt.ylim(0,100)
    print(num_bloods)
    plt.plot(bin_centers, (num_bloods - dvh)/num_bloods*100, c='green')
    plt.grid(True)
    
def plot_positions(blood, color, x_min = 0, x_max = 120, y_min = 0, y_max = 120):
    #NOTE - default arguments were added her to define the dimensions of the graph
    #these dimensions should be the same as the velocity vector field
    '''Plots the initial position of the blood voxels on a 2d grid
    Blood is a list of all the blood voxels objects 
    '''
    x_pos = []
    y_pos = []
    z_pos = []
    for b in blood:
        x = b.get_position().get_x()
        y = b.get_position().get_y()
        z = b.get_position().get_z()
        x_pos.append(x)
        y_pos.append(y)
        z_pos.append(z)
    plt.figure()
    plt.title('Positions of ' +  str(len(blood)) +  ' Blood Voxels')
    plt.xlabel('z positions')
    plt.ylabel('y positions')
    plt.xlim([x_min,x_max]) 
    plt.ylim([y_min,y_max])
    plt.plot(z_pos, y_pos, color)
 
    
def test_plot_positions(num_bloods, vector_field, time, dt,\
                        in_boundary = [(0,10),(30,60),(70,100)],direction ='x'):
    in_bloods = make_blood(num_bloods)
    out_bloods = []
    #plot initial positions
    plot_bloods_3d(in_bloods, c = 'r', m = '^')
    #flow blood and plot final positions
    t = 0
    while t <time:
        bloods_flow(in_bloods, out_bloods, vector_field, dt, in_boundary ,direction)
        t += dt
    plot_bloods_3d(in_bloods, c = 'b',m = 'o')
            
def test_blood():
    '''run the blood simulation and plot various figures'''
    times = [1]
    dose_field = [np.random.rand(50,50,50)]
    time_gaps = [1]
    dt = 0.1
    vector_field = Const_vector_field(50,50,50,1,2,3)
    blood_density = 1
    dose = Dose(dose_field,times,time_gaps)
    bloods = simulate_blood_flow(dose, vector_field, dt, blood_density, \
                        in_boundary = [(0,5),(10,30),(20,35)])
   
    plot_dvh(bloods,blood_density)
def animate_blood():
    '''TODO - animates the flow of blood through a space
    See matplotlab.animate module
    '''
    pass
    

if __name__ == '__main__': 
    start = time.time()
    test_blood()
    end = time.time()
    print("Time to run: ", (end-start), "seconds")

    
    
    
    

##   start time
#    start = time.time()
##    for n in [1]:
##        plot_dvh(doses, time_on, time_off, .1, blood_d = n) #time_on and time_off found in readDoses.py
#    vector_field = Const_vector_field(120,120,120,1,2,3)
#    test_plot_positions(100,vector_field, 20, .1,in_boundary = [(0,10),(30,60),(70,100)],direction ='z')    
#    #stop time
#    end = time.time()
#    print("Time to run: ", (end-start), "seconds")


    

        
