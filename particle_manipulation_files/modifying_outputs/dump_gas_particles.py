import sys, h5py
import numpy as np

file_name = sys.argv[1]

f = h5py.File('./'+file_name+'.hdf5')
f_write = h5py.File(file_name+'.gas.hdf5','w')
f_wg = f_write.create_group('PartType1')
#print f_write.keys()

for data_key in f['PartType1']:    
    f.copy('PartType1/'+data_key,f_wg)

#print f_write['PartType1'].keys()
