

from Blood import Blood
from VectorFields import VectorFields

def create_blood(velocity_field, blood_density = 1):
    '''Takes in a velocity_field object and the blood density
    Creates blood objects within each of the cells that has a nonzero velocity
    '''
    bloods =[]
    #v_sum is a composition of all vector fields to make sure blood is added to
    #every voxel with nonzero elements
    v_sum = velocity_field.get_vx_field() + velocity_field.get_vy_field() + velocity_field.get_vz_field()
    #iterate through every voxel within matrix
    v_iterator = np.nditer(v_sum,flags=['multi_index'])
    for voxel in v_iterator:
        if voxel != 0:
            x,y,z = v_iterator.multi_index
            x = x + np.random.random()
            y = y + np.random.random()
            z = z + np.random.random()
            position = Position(x,y,z)
            blood = Blood(position)
            bloods.append(blood)
    return bloods

def add_dose_for_allblood(all_bloods,dose_matrix, vector_field,dt):
    '''
    add a constant dose of radiation to all blood within the radiation beam
    blood - a list of blood voxel objects
    '''
    vector_field.set_dose_matrix(dose_matrix)
    for i in all_bloods:
        try:
            if vector_field.is_position_in_dose_field(i.get_position()):
                i.current_dose_level(dose_matrix,dt)
                i.add_dose()
        except IndexError:
            pass


def gen_new_blood(out_blood_per_t,list_of_index,prob_field):
    gen_bloods = []
    selected_ind = np.random.choice(list_of_index,out_blood_per_t, \
                                      p = prob_field)
    for i in selected_ind:
        x,y,z = i

        x = x + np.random.random()
        y = y + np.random.random()
        z = z + np.random.random()
        blood = Blood(Position(x, y, z))
        gen_bloods.append(blood)

    return gen_bloods


def bloods_flow(in_bloods, out_bloods, vector_field, dt):
    '''Updates the positions of all the bloods still within the vector field
    in one time step with length dt (a float)
    in_bloods - those still within defined vector field
    out_bloods - those that have left vector field
    vector_field - vector_field object
    in_boundary - the boundary the blood is flowing into, a list cotains x_low,x_high, y_low, y_high,.etc
    '''

    # Move all bloods within the field
    out_blood_count = 0
    for i in in_bloods:
        position = i.get_position()
        x = position.x
        y = position.y
        z = position.z
        vx = vector_field.get_vx_at_position(x, y, z)
        vy = vector_field.get_vy_at_position(x, y, z)
        vz = vector_field.get_vz_at_position(x, y, z)
        if (vx != 0) | ((vy != 0)) | ((vz != 0)):
            new_position = position.get_new_position(vx, vy, vz, dt)
            # if new position still in field, update position
            if vector_field.is_position_in_vector_field(new_position):
                i.find_new_position(new_position)
            else:
                # move that blood from one list to another (all_bloods -> leaving_blods)
                out_bloods.append(i)
                in_bloods.remove(i)
                out_blood_count += 1
                # Generate new bloods, one blood unit for each blood unit out
        else:
            vx, vy, vz = vector_field.get_v_around(x, y, z)
            # try to get back to the vector field by going into the opposite directions
            new_position = position.get_new_position(-vx, -vy, -vz, dt)

            if vector_field.is_position_in_vector_field(new_position):
                i.find_new_position(new_position)
            else:
                # move that blood from one list to another (all_bloods -> leaving_blods)
                out_bloods.append(i)
                in_bloods.remove(i)
                out_blood_count += 1
    in_bloods += gen_new_blood(out_blood_count, list_of_index, prob_field)

#    print('the # bloods flows out' , out_blood_count)
#
def blood_flow_with_beam(in_bloods, out_bloods, vector_field, \
                         dose_matrix, total_time, dt, ):
    """simulate the flow of blood while the radiation beam IS ON, assumes
    a dose_matrix is a two dimension matrix of the dose at each point in a
    2d space
    """
    t = 0;
    while t <= total_time:
        dose_per_time = dose_matrix / total_time  # dose_per_time a dose matrix per unit time
        add_dose_for_allblood(in_bloods, dose_per_time, vector_field, dt)
        bloods_flow(in_bloods, out_bloods, vector_field, dt, )  # TODO - add inboundary as argument later?
        t = t + dt

    return in_bloods

def blood_flow_no_beam(in_bloods, out_bloods, vector_field, \
                       time_gap, dt, in_boundary):
    '''simulate blood flow while radiation beam IS OFF.
    time_gap is the amount of time the beam is off (seconds), dt is length
    of each time step.
    '''
    t = 0
    while t<= time_gap:
        #only update the position of each blood voxel
        bloods_flow(in_bloods, out_bloods, vector_field, dt,in_boundary)
        t += dt
    return in_bloods


def simulate_blood_flow(in_bloods,dose, vector_field, dt, blood_density= 1, \
                        in_boundary = [(0,20),(10,60),(0,20)]):
    '''generate blood objects and velocity vector fields, then runs the simulation
    NOTE - assumes the starting point is when the first field turns on
    dose_fields: a list of dose_field matrices
    times: a list of the total time each dose matrix is 'on'
    time_gaps: the time between doses, while blood is still moving
    '''
    #initialize vector field, make # of blood voxels based on blood density
    num_bloods = len(in_bloods)#total number of the blood
    (x_dim , y_dim ,z_dim) = vector_field.get_dimensions()
    dose_fields = dose.get_dose_field()
    times = dose.get_dose_time()
    time_gaps = dose.get_dose_time_gap()
#    in_bloods = make_blood(num_bloods,x_max =x_dim,y_max=y_dim,z_max=z_dim)
    out_bloods = []
    plot_bloods_3d(in_bloods, c = 'r', m = '^')# initial position of the bloods
    for i in range(len(times)):
            #Add dose
        if (dose_fields[i].shape) != (x_dim , y_dim ,z_dim):
            print('the shape of dose_field does not match with the velocity fields')
            #reshape the field if the dose field and the velocity field do not match
            dose_fields[i] = match_field(dose_fields[i],(x_dim , y_dim ,z_dim))
        in_bloods_after_dose = blood_flow_with_beam(in_bloods, out_bloods, vector_field, dose_fields[i], \
                                                        times[i], dt, in_boundary)
        try:
                #then circulate blood without dose
            in_bloods = blood_flow_no_beam(in_bloods_after_dose, out_bloods, vector_field, time_gaps[i],\
                                               dt,in_boundary)
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

    return in_bloods + out_bloods