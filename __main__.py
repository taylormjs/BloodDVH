from VectorFields import *
from BloodPlot import *
from Position import Position
from Dose import Dose
from BloodSimulation import BloodSimulation
import scipy.io
import numpy as np
import time
import random


def main():
    print('=============================================\n============Blood Simulator==================\n \
    =============================================')
    blood_density = int(input('Choose an integer value for the blood density (1-10): '))
    dt = float(input('Choose a value for the time step: '))
    decision = input('Would you like to vary the velocity field? (y/n): ')
    vx_multiplier, vy_multiplier, vz_multiplier = 1,1,1
    if decision == 'y':
        vx_multiplier = float(input('Choose multiplier for the velocity field in X Direction (Default is 1): '))
        vy_multiplier = float(input('Choose multiplier for the velocity field in Y Direction (Default is 1: '))
        vz_multiplier = float(input('Choose multiplier for the velocity field in Z Direction (Default is 1): '))
    print('-----Loading Velocity Fields-------')
    v_begin_time = time.time()
    mat_v = scipy.io.loadmat('Blood3dLiver/VELOCITYVECSbloodflow.mat')  # read velocity field from mat file
    vx = np.array(mat_v['u'])
    vy = np.array(mat_v['v'])
    vz = np.array(mat_v['w'])
    velocity_field = VectorFields(vx,vy,vz)
    velocity_field.multiplyVelocityField(vx_multiplier, vy_multiplier, vz_multiplier)
    # velocity_field = Const_vector_field(4, 4, 10, .01, .01, 20)
    v_end_time = time.time()
    print('   Time to load: ', v_end_time - v_begin_time)

    print('-----Loading Dose Fields-----')
    d_begin_time = time.time()
    mat_d = scipy.io.loadmat('Blood3dLiver/DOSESbloodflow.mat')  # read dose fields
    dr = np.array(mat_d['dr'])
    # dr = np.random.rand(4, 4, 10)
    doses = Dose([dr], [1] ,[1])
    #TODO - allow user to vary dose
    d_end_time = time.time()
    print('   Time to load: ', d_end_time - d_begin_time)

    print('-----Running Simulation-----')
    simulator = BloodSimulation(velocity_field, doses, blood_density, dt)
    s_begin_time = time.time()
    blood_voxels = simulator.simulate_blood_flow(blood_density, plot_positions=False)
    num_bloods = len(blood_voxels)
    s_end_time = time.time()
    print('   Time to complete simulation: ', s_end_time - s_begin_time)

    print('-----Saving Blood DVH Data-----')
    return (make_dvh(blood_voxels), num_bloods)


    # print('-----Creating Blood DVH Plot-----')
    # # blood_density = simulator.blood_density
    # plot_dvh(blood_voxels, dt, blood_density=blood_density, save_plot=True)
    # print('   Finished plotting and saving plot in DVHGraphs/   ')



if __name__ == '__main__':
    data_sets = [i for i in range(3)]
    num_bloods = [i for i in range(3)]
    for i in range(3):
        data_sets[i], num_bloods[i] = main()
    graphAndSaveDVHPlots(data_sets, .1, num_bloods, ['r','k','b'], ['v=1','v=2','v=5'])
