
def generate_vector_field(dim, initial_p=(10, 10, 0)):
    '''generate a vector field to simulate a real blood flow'''
    x_dim, y_dim, z_dim = dim
    vx_field = np.zeros(dim)
    vy_field = np.zeros(dim)
    vz_field = np.zeros(dim)
    (x, y, z) = initial_p
    (vx, vy, vz) = (0.1, 0.5, 1.0)  # maximum equal to 1
    bloods = []
    while (x < x_dim) & (y < y_dim) & (z < z_dim):
        vx_field[int(x) - 2:int(x) + 2, int(y) - 2:int(y) + 2, int(z)].fill(vx)
        vy_field[int(x) - 2:int(x) + 2, int(y) - 2:int(y) + 2, int(z)].fill(vy)
        vz_field[int(x) - 2:int(x) + 2, int(y) - 2:int(y) + 2, int(z)].fill(vz)
        new_bloods = make_blood(36, int(x) - 2, int(y) - 2, int(z), \
                                int(x) + 2, int(y) + 2, int(z))
        bloods += new_bloods
        x += vx
        y += vy
        z += vz

        # add more condition for xz,yz plane

    vector_field = Vector_field(vx_field, vy_field, vz_field)
    return (vector_field, bloods)