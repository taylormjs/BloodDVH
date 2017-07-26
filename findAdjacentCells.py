from Position import Position
from VectorFields import *
from BloodPlot import *
from Dose import Dose
from BloodSimulation import BloodSimulation
import scipy.io
import numpy as np
import time
import random


matrix = np.array([[[1,2,3],[4,5,6],[7,8,9]],[[10,11,12],[13,14,15],[16,17,18]],[[19,20,21],[22,23,24],[25,26,27]]])
# print(matrix)
#
# print('center = ',matrix[1] [1][1])


def find_value_at_position(position,matrix):
    x,y,z = position
    return matrix[int(x)][int(y)][int(z)]

def get_adjacent_cells(position, s):
    x_coord,y_coord,z_coord = position
    result = []
    if s ==1:
        for x,y,z in [(x_coord+i,y_coord+j,z_coord+k) for i in (-1,0,1) for j in (-1,0,1) \
                  for k in (-1,0,1) if i != 0 or j != 0 or k != 0]: #this is a list of indices surrounding a point
            result.append((x,y,z))
    elif s==2:
        for x,y,z in [(x_coord + i, y_coord + j, z_coord + k) for i in (-2, -1, 0, 1, 2) for j in (-2, -1, 0, 1, 2) \
                        for k in (-2, -1, 0, 1, 2) if i != 0 or j != 0 or k != 0]:
            result.append((x,y,z))
    else:
        result = None
    return result

# l = get_adjacent_cells((5,5,5),2)
# print(l)
# print(len(l))


vx = np.random.randint(0,1,(7,7,7))
vy = np.random.randint(0,1,(7,7,7))
vz = np.random.randint(0,1,(7,7,7))
vx[3][3][2] = 0
vy[1][1][1] = 1

vfield = VectorFields(vx,vy,vz)
prev_pos = Position(2,2,2)
pos = Position(3,3,2)
# print(vfield.get_v_square())
print(vfield.searchShellAndGetNewPosition(pos,prev_pos))


mat_v = scipy.io.loadmat('Blood3dLiver/VELOCITYVECSbloodflow.mat')  # read velocity field from mat file
vx = np.array(mat_v['u'])
vy = np.array(mat_v['v'])
vz = np.array(mat_v['w'])


print(vx)