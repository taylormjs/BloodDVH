from VectorFields import *
from BloodPlot import *
from Position import Position
from Dose import Dose
from BloodSimulation import BloodSimulation
import scipy.io
import numpy as np
import time
import random


def runSimulation(blood_density, dt, multiplier):
    print('=============================================\n============Blood Simulator==================\n=============================================')
    # blood_density = int(input('Choose an integer value for the blood density (1-10): '))
    # dt = float(input('Choose a value for the time step: '))
    # decision = input('Would you like to vary the velocity field? (y/n): ')
    vx_multiplier, vy_multiplier, vz_multiplier = multiplier
    # if decision == 'y':
    #     vx_multiplier = float(input('Choose multiplier for the velocity field in X Direction (Default is 1): '))
    #     vy_multiplier = float(input('Choose multiplier for the velocity field in Y Direction (Default is 1: '))
    #     vz_multiplier = float(input('Choose multiplier for the velocity field in Z Direction (Default is 1): '))
    print('-----Loading Velocity Fields-------')
    v_begin_time = time.time()
    mat_v = scipy.io.loadmat('Blood3dLiver/VELOCITYVECSbloodflow.mat')  # read velocity field from mat file
    vx = np.array(mat_v['u'])
    vy = np.array(mat_v['v'])
    vz = np.array(mat_v['w'])
    # vx = np.random.randint(0,2,(512,512,168))
    # vy = np.random.randint(0, 1, (512, 512, 168))
    # vz = np.random.randint(0, 1, (512, 512, 168))
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
    blood_voxels = simulator.simulate_blood_flow(blood_density, plot_positions=True)
    num_bloods = len(blood_voxels)
    s_end_time = time.time()
    print('   Time to complete simulation: ', s_end_time - s_begin_time)

    print('-----Saving Blood DVH Data-----')
    bins, dvh = make_dvh(blood_voxels)
    return (bins, dvh, num_bloods)  #returns bincenters,dvh (x and y points)


    # print('-----Creating Blood DVH Plot-----')
    # # blood_density = simulator.blood_density
    # plot_dvh(blood_voxels, dt, blood_density=blood_density, save_plot=True)
    # print('   Finished plotting and saving plot in DVHGraphs/   ')

def main():
    blood_density = 1
    dt = .01
    multiplier = [(1,1,1),(5,5,5)] #(30,30,30),(1,1,1)]
    bin_list = []
    dvh_list = []
    data_lengths = []
    for m in range(len(multiplier)):
        bins, dvh, num_bloods = runSimulation(blood_density,dt,multiplier[m])
        bin_list.append(bins)
        dvh_list.append(dvh)
        data_lengths.append(num_bloods)
    return (bin_list, dvh_list, data_lengths)  #a list of tuples, each containing (bins,dvh)





if __name__ == '__main__':
    bin_list, dvh_list, data_lengths = main()
    #((bins,centers),num_bloods)
    graphAndSaveDVHPlots(bin_list, dvh_list, data_lengths, ['r','k'], ['v = 1','v = 5'])
    print('dt = .01')
