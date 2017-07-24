from VectorFields import *
from BloodPlot import *
from Position import Position
from Doses import Doses
from BloodSimulation import BloodSimulation
import scipy.io
import numpy as np
import time
import random


def main():
    print('=============================================\n============Blood Simulator===================\n \
    =============================================')
    blood_density = int(input('Choose an integer value for the blood density (1-10'))
    dt = float(input('Choose a value for the time step'))
    decision = input('Would you like to vary the velocity field? (y/n)')
    if decision == 'y':
        multipliers = tuple(input('Choose multipliers for the velocity field in the form: x_mult, y_mult, z_mult'))
    else:
        multipliers = (1, 1, 1)
    x_mult, y_mult, z_mult = multipliers
    print('-----Loading Velocity Fields-------')
    v_begin_time = time.time()
    mat_v = scipy.io.loadmat('VELOCITYVECSbloodflow.mat')  # read velocity field from mat file
    vx = np.array(mat_v['u'])
    vy = np.array(mat_v['v'])
    vz = np.array(mat_v['w'])
    velocity_field = VectorFields(vx,vy,vz)
    velocity_field.multiplyVelocityField(x_mult, y_mult, z_mult)
    # velocity_field = Const_vector_field(4, 4, 10, .01, .01, 20)
    v_end_time = time.time()
    print('   Time to load: ', v_end_time - v_begin_time)

    print('-----Loading Dose Fields-----')
    d_begin_time = time.time()
    mat_d = scipy.io.loadmat('DOSESbloodflow.mat')  # read dose fields
    dr = np.array(mat_d['dr'])
    # dr = np.random.rand(4, 4, 10)
    doses_object = Doses([dr], [1] ,[1])
    blood_density = 2
    dt = 0.1
    d_end_time = time.time()
    print('   Time to load: ', d_end_time - d_begin_time)

    print('-----Running Simulation-----')
    simulator = BloodSimulation(velocity_field, doses_object, blood_density, dt)
    s_begin_time = time.time()
    blood_voxels = simulator.simulate_blood_flow(blood_density, plot_positions=False)
    s_end_time = time.time()
    print('   Time to complete simulation: ', s_end_time - s_begin_time)

    print('-----Creating Blood DVH Plot-----')
    blood_density = simulator.blood_density
    plot_dvh(blood_voxels, dt, blood_density=blood_density,save_plot=True)
    print('   Finished plotting and saving plot in DVHGraphs/   ')



if __name__ == '__main__':
    main()