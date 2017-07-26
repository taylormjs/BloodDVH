import numpy as np
from Position import Position
import random

class VectorFields(object):
    '''TODO  - make this the superclass of const_vector_field. Make another
    subclass of vector_field called varying_field(?)
    '''
    '''Each position has an associated velocity in x,y, and z directions as 
    well as an associated dose'''

    def __init__(self, vx_field, vy_field, vz_field):  # inputs vx,vy,vz from vtk or csv files
        '''vx_field, vy_field, vz_field are ndarrays
        #TODO - Adjust this later to match up with the velocity fields  in vdx files
        '''
        self.vx_field = vx_field
        self.vy_field = vy_field
        self.vz_field = vz_field
        self.v_square = self.get_v_square()
        self.shape = vx_field.shape
        if len(self.shape) == 2:
            self.x_dim, self.y_dim = self.shape
        else:
            self.x_dim, self.y_dim, self.z_dim = self.shape
            # create three 3-D velocity matrices, each containing the velocity in one direction x,y,or z

    def multiplyVelocityField(self, x_mult, y_mult, z_mult):
        self.vx_field *= x_mult
        self.vy_field *= y_mult
        self.vz_field *= z_mult

    def set_dose_matrix(self, dose_matrix):  # maybe take out the dim parameters later
        '''assume dose matrix a numpy matrix of the same dimensions as velocity field'''
        self.dose_matrix = dose_matrix

    def fit_dose_matrix(self, exp_dose_matrix):  # exp stands for expanded
        '''
        #assumes dose matrix IS NOT the same dimensions as velocity field
        #simply takes in a matrix and truncates the matrix to be the same dimensions
        #as the velocity field
        '''
        pass  # TODO - later

    def get_vx_at_position(self, x, y, z):
        return self.vx_field[int(x)][int(y)][int(z)]

    def get_vy_at_position(self, x, y, z):
        return self.vy_field[int(x)][int(y)][int(z)]

    def get_vz_at_position(self, x, y, z):
        return self.vz_field[int(x)][int(y)][int(z)]


    def get_adjacent_cells(self,x_coord,y_coord,z_coord,shell_num): #shell_num should usually be 1 (positive int)
        surrounding_shell = []
        if shell_num == 1:
            for x, y, z in [(x_coord + i, y_coord + j, z_coord + k) for i in (-shell_num, 0, shell_num) for j in (-shell_num, 0, shell_num) \
                        for k in (-shell_num, 0, shell_num) if i != 0 or j != 0 or k != 0]:  # this is a list of indices surrouding a point
                surrounding_shell.append((x,y,z))
        elif shell_num == 2:
            for x, y, z in [(x_coord + i, y_coord + j, z_coord + k) for i in (-2, -1, 0, 1, 2) for j in (-2, -1, 0, 1, 2) \
                            for k in (-2, -1, 0, 1, 2) if i != 0 or j != 0 or k != 0]:
                surrounding_shell.append((x,y,z))

        return surrounding_shell #or should this be a yield?

    def searchShellAndGetNewPosition(self,position):
        '''takes in the position assumed to be outside of the velocity field, finds the nearest point with a
        nonzero velocity, and returns a position at that point (assumed to be in the velocity field)
        Returns None if no position can be found within the two shells immediately around the cell.'''
        i,j,k = position.get_index_of_position() #TODO - make sure this is the actual cell the blood is in, NOT one off
        v_magnitude_new = self.v_square[i][j][k]
        if v_magnitude_new != 0:
            return Position(i + random.random(), j + random.random(), k + random.random())
        shell_num = 1
        nonzero_velocity_list = []
        while nonzero_velocity_list == []:
            surrounding_shell = self.get_adjacent_cells(i, j, k, shell_num)
            if surrounding_shell != []:
                for i, j, k in surrounding_shell:
                    if self.v_square[i][j][k] != 0:
                        nonzero_velocity_list.append((i,j,k)) #add all viable options to a list
            else:
                print("No near position found")
                return None
            shell_num += 1
        new_position = self.findBestNewPosition(nonzero_velocity_list, position)
        return new_position


    def findBestNewPosition(self, index_list, position):
        best_difference = 100
        best_index = index_list[0]
        for index in index_list:
            x1,y1,z1 = position.get_position() #a tuple of exact position
            x, y, z = index
            difference = (x1 - x)**2 + (y1 - y)**2 + (z1 - z)**2
            if difference < best_difference:
                best_difference = difference
                best_index = index
            else:
                continue
        return Position(best_index)




    def get_v_around(self, x, y, z):
        '''find the velocity of its closet neighbors and return the average'''
        vx_around = []
        vy_around = []
        vz_around = []
        for i in [-1, 1]:
            for j in [-1, 1]:
                for k in [-1, 1]:
                    x_new = int(x + i)
                    y_new = int(y + j)
                    z_new = int(z + k)
                    if self.is_position_in_vector_field(Position(x_new, y_new, z_new)):
                        vx_around.append(self.vx_field[x_new, y_new, z_new])
                        vy_around.append(self.vy_field[x_new, y_new, z_new])
                        vz_around.append(self.vz_field[x_new, y_new, z_new])
                    else:
                        vx_around.append(0)
                        vy_around.append(0)
                        vz_around.append(0)
        v_around = (np.array(vx_around).mean(), np.array(vy_around).mean(), \
                    np.array(vz_around).mean())

        return v_around

    def get_vx_field(self):
        return self.vx_field

    def get_vy_field(self):
        return self.vy_field

    def get_vz_field(self):
        return self.vz_field

    def get_v_square(self):
        return self.vx_field ** 2 + self.vy_field ** 2 + self.vz_field ** 2

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

    def find_gen_field(self):
        '''this method is used to decided where to generate new blood as one leave the field
        '''
        v_square = self.get_v_square()
        surface1 = self.get_vx_field()[0, :, :]  # the surface where x = 0
        surface2 = self.get_vx_field()[-1, :, :]  # x = -1
        surface3 = self.get_vy_field()[:, 0, :]  # y = 0
        surface4 = self.get_vy_field()[:, -1, :]  # y = -1
        surface5 = self.get_vz_field()[:, :, 65]  # z = 0
        surface6 = self.get_vz_field()[:, :, 114]  # z = -1
        list_of_index = []
        mini_field =[]
        v_iterator1 = np.nditer(surface1, flags=['multi_index'])
        v_iterator2 = np.nditer(surface2, flags=['multi_index'])
        v_iterator3 = np.nditer(surface3, flags=['multi_index'])
        v_iterator4 = np.nditer(surface4, flags=['multi_index'])
        v_iterator5 = np.nditer(surface5, flags=['multi_index'])
        v_iterator6 = np.nditer(surface6, flags=['multi_index'])
        for voxel in v_iterator1:
            if voxel > 0:
                x = 0
                (y, z) = v_iterator1.multi_index
                list_of_index.append((x, y, z))
                mini_field.append(v_square[x, y, z])

        for voxel in v_iterator2:
            if voxel < 0:
                x = -1
                (y, z) = v_iterator2.multi_index
                list_of_index.append((x, y, z))
                mini_field.append(v_square[x, y, z])

        for voxel in v_iterator3:
            if voxel > 0:
                y = 0
                (x, z) = v_iterator3.multi_index
                list_of_index.append((x, y, z))
                mini_field.append(v_square[x, y, z])

        for voxel in v_iterator4:
            if voxel < 0:
                y = -1
                (x, z) = v_iterator4.multi_index
                list_of_index.append((x, y, z))
                mini_field.append(v_square[x, y, z])

        for voxel in v_iterator5:
            if voxel > 0:
                z = 65
                (x, y) = v_iterator5.multi_index
                list_of_index.append((x, y, z))
                mini_field.append(v_square[x, y, z])

        for voxel in v_iterator6:
            if voxel < 0:
                z = 114
                (x, y) = v_iterator6.multi_index
                list_of_index.append((x, y, z))
                mini_field.append(v_square[x, y, z])
        total = sum(mini_field)
        prob_field = [i/total for i in mini_field]  # need to check if the number are too small

        return (list_of_index, prob_field)

    def is_position_in_dose_field(self, position):
        '''Returns true if the input position is within the dose_field
        '''
        x = int(position.get_x())
        y = int(position.get_y())  # Changed math.floor to int to save time
        z = int(position.get_z())
        return self.dose_matrix[x][y][z] != 0  # if dose is nonzero, must be in dose field

    def is_position_in_vector_field(self, position):
        '''Returns true if position is within the blood vessel
        '''
        x = position.get_x()
        y = position.get_y()
        z = position.get_z()
        return (0 <= x < self.x_dim and 0 <= y < self.y_dim and 65 <= z < 115 )
        #TODO z has been hard coded to work with the given matrix
        # less than or equal to? <=

    def __str__(self):
        return "Vector_field w/ dim: " + str((self.x_dim, self.y_dim, \
                                              self.z_dim))

class Const_vector_field(VectorFields):
    '''Each position has an associated velocity in x,y, and z directions as
        well as an associated dose'''
    # NOTE - vy dose not have a default value of 0 in 2d version
    def __init__(self, x_dim, y_dim, z_dim, vx, vy, vz):  # change the order of these later
        '''Assumes x_dim is a scalar of the number of units in our matrix
        In one dimension, y_dim and z_dim should just be 1
        #TODO - Adjust this later to match up with the velocity fields  in vdx files
        '''
        self.vx, self.vy, self.vz = vx, vy, vz
        self.velocity = (self.vx, self.vy, self.vz)
        velocity = (self.vx, self.vy, self.vz)
        vx_field = np.zeros((x_dim, y_dim, z_dim)) + self.vx
        vy_field = np.zeros((x_dim, y_dim, z_dim)) + self.vy
        vz_field = np.zeros((x_dim, y_dim, z_dim)) + self.vz
        VectorFields.__init__(self, vx_field, vy_field, vz_field)