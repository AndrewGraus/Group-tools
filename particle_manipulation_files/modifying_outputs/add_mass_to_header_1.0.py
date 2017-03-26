#This is just a simple script to add the mass 
#of the hi res dark matter particles if you don't
#do this then rockstar won't run on it.

import h5py, sys
import numpy as np

if len(sys.argv) == 1:
    print 'usage: python add_mass_to_header_1.0.py <any number of hdf5 files to modify separated by spaces>'

if len(sys.argv) == 2:
    #its only one file
    f = h5py.File(str(sys.argv[1]))
    m_hi = f['PartType1']['Masses'][:][0]
    m_list = np.asarray([0.0,m_hi,0.0,0.0,0.0,0.0])
    f['Header'].attrs['MassTable'] = m_list

elif len(sys.argv) > 2:
    for part_file in sys.argv[2:len(sys.argv)]:
        print part_file
        f =h5py.File(part_file)
        m_hi = f['PartType1']['Masses'][:][0]
        m_list = np.asarray([0.0,m_hi,0.0,0.0,0.0,0.0])
        f['Header'].attrs['MassTable'] = m_list
