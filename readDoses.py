import scipy.io
import numpy as np 
mat = scipy.io.loadmat('DOSEFIELDS.mat')


f1 = np.array([mat['f1']]) #has dimenstion 1*75*100
f2 = np.array([mat['f2']]) #has dimenstion 1*75*100
f3 = np.array([mat['f3']]) #has dimenstion 1*75*100
def stack_field(field,times=1,axis=0):
    '''stack a field many times to make it thicker'''
    f1 = np.copy(field) 
    new_field = np.copy(field)
    for i in range(times-1):
       new_field = np.concatenate((new_field,f1),axis)
    return new_field
        

        
t1 = mat['t1'][0][0]
t2 = mat['t2'][0][0]
t3 = mat['t3'][0][0]
gap12 = mat['gap12'][0][0]
gap23 = mat['gap23'][0][0]
total_time = t1 + t2 +t3 + gap12 + gap23

time_on = [t1,t2,t3]
time_off = [gap12,gap23]

doses = [f1,f2,f3]

#creat a class for Dose

class Dose(object):
    ''' A Blood voxel with an position and dose
    '''
    def __init__(self, dose_field , time_on, time_gap = [0]):
        self.dose_field = dose_field # a list of dose_field
        self.time_on = time_on
        self.time_gap = time_gap
        (self.x_dim, self.y_dim, self.z_dim) = self.dose_field[0].shape
    
        
    def get_dose_field(self):
        return self.dose_field  
    
    def get_shape(self):
        return (self.x_dim, self.y_dim, self.z_dim)
    
    def get_dose_time(self):
        return self.time_on
    
    def get_dose_time_gap(self):
        return self.time_gap
    
    def add_new_dose(self, dose):
        new_dose = dose.get_dose_field()
        new_time_on = dose.get_dose_time()
        new_time_gap = dose.get_dose_time_gap()
        self.dose_field += new_dose 
        self.time_on += new_time_on
        self.time_gap += new_time_gap
        
    
    def rotate(self, axis = 0):
        '''rotate the dose field along x,y or z axis'''
        pass
    
    def match_field(self, dimensions):
        '''resize the dose_field to have required dimensions''' 
        if (self.x_dim, self.y_dim, self.z_dim) != dimensions:
            new_dose_field = []
            for i in self.dose_field: 
                new_field = np.zeros(dimensions)
                new_field[0:self.x_dim, 0:self.y_dim, 0:self.z_dim] = i
                new_dose_field.append(new_field)                
            self.dose_field = new_dose_field
        
        
        