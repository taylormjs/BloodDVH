from VectorFields import VectorFields
from BloodPlot import *
from Position import Position
from Dose import Dose
from BloodSimulation import  BloodSimulation
import scipy.io
import numpy as np



def main():

    simulation = BloodSimulation()
    mat_v = scipy.io.loadmat('VELOCITYVECSbloodflow.mat')  # read velocity field from mat file
    vx = np.array(mat_v['u'])
    vy = np.array(mat_v['v'])
    vz = np.array(mat_v['w'])
    mat_d = scipy.io.loadmat('DOSESbloodflow.mat')  # read dose field
    velocity_field = VectorFields(vx,vy,vz)
    mat_d = scipy.io.loadmat('DOSESbloodflow.mat')  # read dose field
    dr = np.array(mat_d['dr'])
    doses_object = Dose([dr],[1],[1])
    blood_density = 1
    dt = 0.1
    simulator = BloodSimulation( velocity_field, doses_object, blood_density, dt)