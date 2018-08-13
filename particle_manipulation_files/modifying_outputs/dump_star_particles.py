import sys, h5py
import numpy as np

file_name = sys.argv[1]

f = h5py.File('./'+file_name+'.hdf5')
f_write = h5py.File(file_name+'.stars.hdf5','w')
f_wg = f_write.create_group('PartType4')
print f_write.keys()

for data_key in f['PartType4']:    
    f.copy('PartType4/'+data_key,f_wg)

print f_write['PartType4'].keys()
