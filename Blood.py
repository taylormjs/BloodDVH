from Position import Position

class Blood(object):
    ''' A Blood voxel with an position and dose'''
    def __init__(self, position, init_dose=0):
        self.position = position
        self.dose_added = init_dose
        # self.getting_dose = False;
        self.dose_reciving = 0

    def get_dose(self):
        return self.dose_added

    def add_dose(self):
        self.dose_added += self.dose_reciving

    def get_position(self):
        return self.position

    def find_new_position(self, position):
        self.position = position

    def current_dose_level(self, dose_matrix):
        '''return dose that each unit of blood gets in a dose_matrix
        '''
        self.dose_reciving = 0
        position = self.get_position()
        x = int(position.x)  # TODO - these may need to be getter functions later
        y = int(position.y)  # ie y = position.get_y()
        z = int(position.z)
        dose = dose_matrix[x][y][z]
        self.dose_reciving = dose

    def __str__(self):
        return "Blood at " + str(self.get_position()) + " w/ Dose of " + \
               str(self.get_dose())

