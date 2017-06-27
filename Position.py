class Position(object):
    '''A position within the blood vessel.
    '''

    def __init__(self, x, y, z):
        # assumes x, y, and z are float, representing a voxel within the vector field
        # NOT REPRESENTING an exact location in the x, y, and z
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

    def get_new_position(self, vx, vy, vz, dt):  # change to units rather than actual distance
        # dt is the change in time dt
        old_x, old_y, old_z = self.get_x(), self.get_y(), self.get_z()
        dx = vx * dt
        dy = vy * dt
        dz = vz * dt
        # make the new positions units within the vector field, hence use math.floor
        new_x, new_y, new_z = (old_x + dx), (old_y + dy), (old_z + dz)
        return Position(new_x, new_y, new_z)

    def set_velocity(self, velocity):
        '''Useful for reading in velocity files
        '''
        self.velocity = velocity

    def get_velocity(self):
        return self.velocity

    def __str__(self):
        return "Position: " + str((self.x, self.y, self.z))
        # '(' + str(self.x) + str(self.y) + str(self.z) + ')'