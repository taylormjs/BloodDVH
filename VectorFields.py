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
        v_square = vector_field.v_square
        surface1 = self.get_vx_field()[0, :, :]  # the surface where x = 0
        surface2 = self.get_vx_field()[-1, :, :]  # x = -1
        surface3 = self.get_vy_field()[:, 0, :]  # y = 0
        surface4 = self.get_vy_field()[:, -1, :]  # y = -1
        surface5 = self.get_vz_field()[:, :, 0]  # z = 0
        surface6 = self.get_vz_field()[:, :-1]  # z = -1
        list_of_index = []
        mini_field = []
        v_iterator1 = np.nditer(surface1, flags=['multi_index'])
        v_iterator2 = np.nditer(surface2, flags=['multi_index'])
        v_iterator3 = np.nditer(surface3, flags=['multi_index'])
        v_iterator4 = np.nditer(surface4, flags=['multi_index'])
        v_iterator5 = np.nditer(surface5, flags=['multi_index'])
        v_iterator6 = np.nditer(surface6, flags=['multi_index'])
        for voxel in v_iterator1:
            if voxel > 0:
                (x, y, z) = v_iterator1.multi_index
                list_of_index.append((x, y, z))
                mini_field.append(v_square[x, y, z])

        for voxel in v_iterator2:
            if voxel < 0:
                (x, y, z) = v_iterator2.multi_index
                list_of_index.append((x, y, z))
                mini_field.append([v_squarex, y, z])

        for voxel in v_iterator3:
            if voxel > 0:
                (x, y, z) = v_iterator3.multi_index
                list_of_index.append((x, y, z))
                mini_field.append(v_square[x, y, z])

        for voxel in v_iterator4:
            if voxel < 0:
                (x, y, z) = v_iterator4.multi_index
                list_of_index.append((x, y, z))
                mini_field.append(v_square[x, y, z])

        for voxel in v_iterator5:
            if voxel > 0:
                (x, y, z) = v_iterator5.multi_index
                list_of_index.append((x, y, z))
                mini_field.append(v_square[x, y, z])

        for voxel in v_iterator6:
            if voxel < 0:
                (x, y, z) = v_iterator6.multi_index
                list_of_index.append((x, y, z))
                mini_field.append([x, y, z])

        total = mini_field.sum()
        prob_field = mini_field / total  # need to check if the number are too small

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
        return (0 <= x < self.x_dim and 0 <= y < self.y_dim and 0 <= z < self.z_dim)
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
        Vector_field.__init__(self, vx_field, vy_field, vz_field)