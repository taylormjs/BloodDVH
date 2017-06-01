import scipy.io
import numpy as np 
mat = scipy.io.loadmat('DOSEFIELDS.mat')


f1 = [mat['f1']] #has dimenstion 1*75*100
f2 = [mat['f2']] #has dimenstion 1*75*100
f2 = [mat['f3']] #has dimenstion 1*75*100
t1 = mat['t1'][0][0]
t2 = mat['t2'][0][0]
t3 = mat['t3'][0][0]
gap12 = mat['gap12'][0][0]
gap23 = mat['gap23'][0][0]
total_time = t1 + t2 +t3 + gap12 + gap23

time_on = [t1,t2,t3]
time_off = [gap12,gap23]
