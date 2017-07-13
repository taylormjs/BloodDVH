from Blood import Blood
from VectorFields import VectorFields
from Doses import Doses

class BloodSimulation(object):
    def __init__(self, velocity_field, doses_object, blood_density, dt):
        #TODO - should dose_matrix be a list of dose_matrices? Like a Doses object?
        #TODO - include time gaps and times beams is on in initialization or just in simulate_blood_flow method?
        self.velocity_field = velocity_field
        self.doses = doses_object
        self.blood_density = blood_density
        self.dt = dt
        self.in_bloods = []
        self.out_bloods = []

    def create_blood(self):
        '''Takes in a velocity_field object and the blood density
        Creates blood objects within each of the cells that has a nonzero velocity
        '''
        #TODO - include blood_density in this calculation?
        bloods =[]
        #v_sum is a composition of all vector fields to make sure blood is added to
        #every voxel with nonzero elements
        v_sum = self.velocity_field.get_vx_field() + self.velocity_field.get_vy_field() \
                + self.velocity_field.get_vz_field()
        #iterate through every voxel within matrix
        v_iterator = np.nditer(v_sum, flags=['multi_index'])
        for voxel in v_iterator:
            if voxel != 0:
                x,y,z = v_iterator.multi_index
                x = x + np.random.random()
                y = y + np.random.random()
                z = z + np.random.random()
                position = Position(x,y,z)
                blood = Blood(position)
                self.in_bloods.append(blood)
        return self.in_bloods #TODO - may not need to be returned

    def add_dose_for_allblood(self, dose_matrix):
        '''
        add a constant dose of radiation to all blood within the radiation beam
        blood - a list of blood voxel objects
        '''
        self.vector_field.set_dose_matrix(dose_matrix)
        for i in self.in_bloods:
            try:
                if self.vector_field.is_position_in_dose_field(i.get_position()):
                    i.current_dose_level(dose_matrix, self.dt)
                    i.add_dose()
            except IndexError:
                pass


    def generate_new_blood(self,out_blood_per_t,list_of_index,prob_field): #TODO - modify this function
        generated_bloods = []
        selected_ind = np.random.choice(list_of_index,out_blood_per_t, \
                                          p = prob_field)
        for i in selected_ind:
            x,y,z = i
            x = x + np.random.random()
            y = y + np.random.random()
            z = z + np.random.random()
            blood = Blood(Position(x, y, z))
            generated_bloods.append(blood)

        return generated_bloods


    def bloods_flow(self):
        '''Updates the positions of all the bloods still within the vector field
        in one time step with length dt (a float)
        in_bloods - those still within defined vector field
        out_bloods - those that have left vector field
        vector_field - vector_field object
        in_boundary - the boundary the blood is flowing into, a list cotains x_low,x_high, y_low, y_high,.etc
        '''
        #TODO - should in_bloods and out_bloods be attributes of simulation class?
        # Move all bloods within the field
        out_blood_count = 0
        for blood in self.in_bloods:
            position = blood.get_position()
            x = position.x
            y = position.y
            z = position.z
            vx = self.vector_field.get_vx_at_position(x, y, z)
            vy = self.vector_field.get_vy_at_position(x, y, z)
            vz = self.vector_field.get_vz_at_position(x, y, z)
            if (vx != 0) | ((vy != 0)) | ((vz != 0)): #TODO - is this proper python syntax?
                new_position = position.get_new_position(vx, vy, vz,self.dt)
                # if new position still in field, update position
                if self.vector_field.is_position_in_vector_field(new_position):
                    blood.find_new_position(new_position)
                else:
                    # move that blood from one list to another (all_bloods -> leaving_blods)
                    self.out_bloods.append(blood) #TODO - this will present a problem later if blood leaves out the walls
                    self.in_bloods.remove(blood)
                    out_blood_count += 1
                    # Generate new bloods, one blood unit for each blood unit out
            else:
                vx, vy, vz = vector_field.get_v_around(x, y, z)
                # try to get back to the vector field by going into the opposite directions
                new_position = position.get_new_position(-vx, -vy, -vz, dt)

                if self.vector_field.is_position_in_vector_field(new_position):
                    blood.find_new_position(new_position)
                else:
                    # move that blood from one list to another (all_bloods -> leaving_bloods)
                    self.out_bloods.append(blood)
                    self.in_bloods.remove(blood)
                    out_blood_count += 1
        self.in_bloods += self.generate_new_blood(out_blood_count, list_of_index, prob_field)
        #TODO - where are list_of_index and prob_field defined?

    #    print('the # bloods flows out' , out_blood_count)
    #
    def blood_flow_with_beam(self, dose_matrix, total_time):
        """simulate the flow of blood while the radiation beam IS ON, assumes
        a dose_matrix is a two dimension matrix of the dose at each point in a
        2d space
        """
        t = 0;
        while t <= total_time:
            dose_per_time = dose_matrix / total_time  # dose_per_time a dose matrix per unit time
            self.add_dose_for_allblood(dose_matrix) #TODO - should this have dose per time argument? was deleted from add_dose_for_allbloods method
            self.bloods_flow()  # TODO - add inboundary as argument later?
            t = t + self.dt

        return self.in_bloods #TODO - may not need to be returned if self.in_bloods is an attribute (same with no beam method)

    def blood_flow_no_beam(self, time_gap, in_boundary):
        '''simulate blood flow while radiation beam IS OFF.
        time_gap is the amount of time the beam is off (seconds), dt is length
        of each time step.
        '''
        t = 0
        while t<= time_gap:
            #only update the position of each blood voxel
            self.bloods_flow() #TODO - add inboundary argument later?
            t += self.dt
        return self.in_bloods


    def simulate_blood_flow(self, in_boundary = [(0,20),(10,60),(0,20)]):
        '''generate blood objects and velocity vector fields, then runs the simulation
        NOTE - assumes the starting point is when the first field turns on
        dose_fields: a list of dose_field matrices
        times: a list of the total time each dose matrix is 'on'
        time_gaps: the time between doses, while blood is still moving
        '''
        #initialize vector field, make # of blood voxels based on blood density
        num_bloods = len(self.in_bloods)#total number of the blood
        (x_dim , y_dim ,z_dim) = self.vector_field.get_dimensions()
        dose_fields = self.doses.get_dose_fields()
        times = self.doses.get_dose_time()
        time_gaps = self.doses.get_dose_time_gap()
    #    in_bloods = make_blood(num_bloods,x_max =x_dim,y_max=y_dim,z_max=z_dim)
        out_bloods = []
        plot_bloods_3d(in_bloods, c = 'r', m = '^')# initial position of the bloods
        for i in range(len(times)):
                #Add dose
            if (dose_fields[i].shape) != (x_dim , y_dim ,z_dim):
                print('the shape of dose_field does not match with the velocity fields')
                #reshape the field if the dose field and the velocity field do not match
                dose_fields[i] = match_field(dose_fields[i],(x_dim , y_dim ,z_dim))
            self.blood_flow_with_beam(dose_fields[i], times[i]) #TODO - should this have an inboundary arg?
            try:
                #then circulate blood without dose
                self.blood_flow_no_beam(time_gaps[i], in_boundary) #TODO - note an in_bloods_after_dose arg was here before
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

        return self.in_bloods + self.out_bloods