#NOTE FROM AG: This is a script that takes a gizmo hdf5 file
#and converts it to a gadget binary. This is necessary for
#running AHF because it is not currently configured to run
#on AHF files
#
#Written by Shea (I think?)

#!/usr/bin/python

import h5py
from struct import pack,unpack
import sys,os
from numpy import *


if len(sys.argv) != 3:
    print "Usage:  python gh52gb.py <input file> <output name>"
    sys.exit(1337)

fname = sys.argv[1]
outname = sys.argv[2]


f = h5py.File(fname)
out = open(outname,'w')

ptypes = f.keys()
ptypes.remove('Header')

fhead = f['Header'].attrs
nfile = fhead['NumPart_ThisFile'][:]
ntot = sum(nfile)
masstable = fhead['MassTable'][:]

poses = []
vels = []
ids = []
masses = []     #this one might not get filled in at all
if 'PartType0' in ptypes: 
    gas = True
    if 'RateOfChangeOfEntropy' in f['PartType0'].keys():
        dsdt = True
    else:
        dsdt = False
else:   
    gas = False
    dsdt = False
if 'Potential' in f[ptypes[0]].keys():    pot = True
else:   pot = False
if 'Acceleration' in f[ptypes[0]].keys():    accel = True
else:    accel = False
if 'Timestep' in f[ptypes[0]].keys():    time = True
else:   time = False

if 'PartType1' in ptypes:    halo = True
else:   halo = False
if 'PartType2' in ptypes:    disk = True
else:   disk = False
if 'PartType3' in ptypes:    bulge = True
else:   bulge = False
if 'PartType4' in ptypes:    star = True
else:   star = False
if 'PartType5' in ptypes:    bndry = True
else:   bndry = False

#start with the header:
packedheader = ""
ar = array(fhead['NumPart_ThisFile'][:],dtype='I')
packedheader = packedheader + ar.tostring()

mtable = array(fhead['MassTable'][:],dtype='d')
nmass = 0
if gas:
    if mtable[0] == 0:
        if (f['PartType0']['Masses'][:] == f['PartType0']['Masses'][0]).all():
            mtable[0] = f['PartType0']['Masses'][0]
        else:
            nmass = nmass + f['PartType0']['Masses'].shape[0]
if halo:
    if mtable[1] == 0:
        if (f['PartType1']['Masses'][:] == f['PartType1']['Masses'][0]).all():
            mtable[1] = f['PartType1']['Masses'][0]
        else:
            nmass = nmass + f['PartType1']['Masses'].shape[0]
if disk:
    if mtable[2] == 0:
        if (f['PartType2']['Masses'][:] == f['PartType2']['Masses'][0]).all():
            mtable[2] = f['PartType2']['Masses'][0]
        else:
            nmass = nmass + f['PartType2']['Masses'].shape[0]
if bulge:
    if mtable[3] == 0:
        if (f['PartType3']['Masses'][:] == f['PartType3']['Masses'][0]).all():
            mtable[3] = f['PartType3']['Masses'][0]
        else:
            nmass = nmass + f['PartType3']['Masses'].shape[0]
if star:
    if mtable[4] == 0:
        if (f['PartType4']['Masses'][:] == f['PartType4']['Masses'][0]).all():
            mtable[4] = f['PartType4']['Masses'][0]
        else:
            nmass = nmass + f['PartType4']['Masses'].shape[0]
if bndry:
    if mtable[5] == 0:
        if (f['PartType5']['Masses'][:] == f['PartType5']['Masses'][0]).all():
            mtable[5] = f['PartType5']['Masses'][0]
        else:
            nmass = nmass + f['PartType5']['Masses'].shape[0]
        
packedheader = packedheader + mtable.tostring()

ar = array(fhead['Time'],dtype='d')
packedheader = packedheader + ar.tostring()

ar = array(fhead['Redshift'],dtype='d')
packedheader = packedheader + ar.tostring()

ar = array(fhead['Flag_Sfr'],dtype='i')
packedheader = packedheader + ar.tostring()

ar = array(fhead['Flag_Feedback'],dtype='i')
packedheader = packedheader + ar.tostring()

ar = array(fhead['NumPart_Total'][:],dtype='i')
packedheader = packedheader + ar.tostring()


ar = array(fhead['Flag_Cooling'],dtype='i')
packedheader = packedheader + ar.tostring()

ar = array(fhead['NumFilesPerSnapshot'],dtype='i')
packedheader = packedheader + ar.tostring()

ar = array(fhead['BoxSize'],dtype='d')
packedheader = packedheader + ar.tostring()

ar = array(fhead['Omega0'],dtype='d')
packedheader = packedheader + ar.tostring()

ar = array(fhead['OmegaLambda'],dtype='d')
packedheader = packedheader + ar.tostring()

ar = array(fhead['HubbleParam'],dtype='d')
packedheader = packedheader + ar.tostring()

ar = array(fhead['Flag_StellarAge'],dtype='i')
packedheader = packedheader + ar.tostring()

ar = array(fhead['Flag_Metals'],dtype='i')
packedheader = packedheader + ar.tostring()

try:
    ar = array(fhead['NumPart_Total_HW'][:],dtype='i')
except KeyError:
    ar = array(fhead['NumPart_Total_HighWord'][:],dtype='i')
packedheader = packedheader + ar.tostring()


if 'Flag_Entropy_ICs' in fhead.keys():
    try:        #This is an array in at least one file that I have, so attempt to do it that way, and if it fails, do it as a float
        ar = array(fhead['Flag_Entropy_ICs'][:],dtype='i')
        packedheader = packedheader + ar.tostring()
    except TypeError:
        ar = array(fhead['Flag_Entropy_ICs'],dtype='i')
        packedheader = packedheader + ar.tostring()
else:
    print "Using Flag_IC_Info instead of Flag_Entropy_ICs."
    ar = array(fhead['Flag_IC_Info'],dtype='i')
    packedheader = packedheader + ar.tostring()


header_bytes_left = 256 - len(packedheader)
for i in range(header_bytes_left):
    packedheader = packedheader + pack('<x')

#now to write it out into a binary file
out.write(pack('<I',256))
out.write(packedheader)
out.write(pack('<I',256))

#Now to do coordinates, order of gas halo disk bulge star boundary
vec_size = array([12*ntot],dtype='I')
out.write(vec_size.tostring())

if gas:
    ar = array(f['PartType0']['Coordinates'][:],dtype='f') #The dtype is so that tostring doesn't give me float64
    out.write(ar.tostring())
if halo:
    ar = array(f['PartType1']['Coordinates'][:],dtype='f')
    out.write(ar.tostring())
if disk:
    ar = array(f['PartType2']['Coordinates'][:],dtype='f')
    out.write(ar.tostring())
if bulge:
    ar = array(f['PartType3']['Coordinates'][:],dtype='f')
    out.write(ar.tostring())
if star:
    ar = array(f['PartType4']['Coordinates'][:],dtype='f')
    out.write(ar.tostring())
if bndry:
    ar = array(f['PartType5']['Coordinates'][:],dtype='f')
    out.write(ar.tostring())
out.write(vec_size.tostring())

print "Finished with coordinates"

#Now for velocities
out.write(vec_size.tostring())
if gas:
    ar = array(f['PartType0']['Velocities'][:],dtype='f') #The dtype is so that tostring doesn't give me float64
    out.write(ar.tostring())
if halo:
    ar = array(f['PartType1']['Velocities'][:],dtype='f')
    out.write(ar.tostring())
if disk:
    ar = array(f['PartType2']['Velocities'][:],dtype='f')
    out.write(ar.tostring())
if bulge:
    ar = array(f['PartType3']['Velocities'][:],dtype='f')
    out.write(ar.tostring())
if star:
    ar = array(f['PartType4']['Velocities'][:],dtype='f')
    out.write(ar.tostring())
if bndry:
    ar = array(f['PartType5']['Velocities'][:],dtype='f')
    out.write(ar.tostring())
out.write(vec_size.tostring())

print "Finished with velocities"

#Now for particle IDs:
float_size = array([4*ntot],dtype='I')
out.write(float_size.tostring())
if gas:
    ar = array(f['PartType0']['ParticleIDs'][:],dtype='I') #The dtype is so that tostring doesn't give me float64
    out.write(ar.tostring())
if halo:
    ar = array(f['PartType1']['ParticleIDs'][:],dtype='I')
    out.write(ar.tostring())
if disk:
    ar = array(f['PartType2']['ParticleIDs'][:],dtype='I')
    out.write(ar.tostring())
if bulge:
    ar = array(f['PartType3']['ParticleIDs'][:],dtype='I')
    out.write(ar.tostring())
if star:
    ar = array(f['PartType4']['ParticleIDs'][:],dtype='I')
    out.write(ar.tostring())
if bndry:
    ar = array(f['PartType5']['ParticleIDs'][:],dtype='I')
    out.write(ar.tostring())
out.write(float_size.tostring())

print "Finished with particle IDs"

#Now I have to check if there are variable particle masses anywhere in the file
#nmass = 0
#for p in ptypes:
#    if "Masses" in f[p].keys():
#        nmass = nmass + len(f[p]['Masses'][:])

if nmass > 0:
    nmass_size = array([4*nmass],dtype='I')
    out.write(nmass_size.tostring())
    if gas:
        #if 'Masses' in f['PartType0'].keys():
        if mtable[0] == 0:
            ar = array(f['PartType0']['Masses'][:],dtype='f')
            out.write(ar.tostring())
    if halo:
        #if 'Masses' in f['PartType1'].keys():
        if mtable[1] == 0:
            ar = array(f['PartType1']['Masses'][:],dtype='f')
            out.write(ar.tostring())
    if disk:
        #if 'Masses' in f['PartType2'].keys():
        if mtable[2] == 0:
            ar = array(f['PartType2']['Masses'][:],dtype='f')
            out.write(ar.tostring())
    if bulge:
        #if 'Masses' in f['PartType3'].keys():
        if mtable[3] == 0:
            ar = array(f['PartType3']['Masses'][:],dtype='f')
            out.write(ar.tostring())
    if star:
        #if 'Masses' in f['PartType4'].keys():
        if mtable[4] == 0:
            ar = array(f['PartType4']['Masses'][:],dtype='f')
            out.write(ar.tostring())
    if bndry:
        #if 'Masses' in f['PartType5'].keys():
        if mtable[5] == 0:
            ar = array(f['PartType5']['Masses'][:],dtype='f')
            out.write(ar.tostring())
    out.write(pack('<I',4*nmass))
    print "Done writing masses for {0} particles".format(nmass)

#Now for gas specific stuff
if gas:
    ngas = nfile[0]
    path = f['PartType0']
    gas_size = array([4*ngas],dtype='I')
    
    ar = array(path['InternalEnergy'][:],dtype='f')
    out.write(gas_size.tostring())
    out.write(ar.tostring())
    out.write(gas_size.tostring())
    
    ar = array(path['Density'][:],dtype='f')
    out.write(gas_size.tostring())
    out.write(ar.tostring())
    out.write(gas_size.tostring())
    
    ar = array(path['SmoothingLength'][:],dtype='f')
    out.write(gas_size.tostring())
    out.write(ar.tostring())
    out.write(gas_size.tostring())

#Now for Makefile enabled things:
if pot:
    out.write(float_size.tostring())
    if gas:
        ar = array(f['PartType0']['Potential'][:],dtype='f') #The dtype is so that tostring doesn't give me float64
        out.write(ar.tostring())
    if halo:
        ar = array(f['PartType1']['Potential'][:],dtype='f')
        out.write(ar.tostring())
    if disk:
        ar = array(f['PartType2']['Potential'][:],dtype='f')
        out.write(ar.tostring())
    if bulge:
        ar = array(f['PartType3']['Potential'][:],dtype='f')
        out.write(ar.tostring())
    if star:
        ar = array(f['PartType4']['Potential'][:],dtype='f')
        out.write(ar.tostring())
    if bndry:
        ar = array(f['PartType5']['Potential'][:],dtype='f')
        out.write(ar.tostring())
    out.write(float_size.tostring())
    
if accel:
    out.write(vec_size.tostring())
    if gas:
        ar = array(f['PartType0']['Acceleration'][:],dtype='f') #The dtype is so that tostring doesn't give me float64
        out.write(ar.tostring())
    if halo:
        ar = array(f['PartType1']['Acceleration'][:],dtype='f')
        out.write(ar.tostring())
    if disk:
        ar = array(f['PartType2']['Acceleration'][:],dtype='f')
        out.write(ar.tostring())
    if bulge:
        ar = array(f['PartType3']['Acceleration'][:],dtype='f')
        out.write(ar.tostring())
    if star:
        ar = array(f['PartType4']['Acceleration'][:],dtype='f')
        out.write(ar.tostring())
    if bndry:
        ar = array(f['PartType5']['Acceleration'][:],dtype='f')
        out.write(ar.tostring())
    out.write(vec_size.tostring())
    
if dsdt and gas:
    out.write(gas_size.tostring())
    ar = array(f['PartType0']['RateOfChangeOfEntropy'][:],dtype='f')
    out.write(ar.tostring())
    out.write(gas_size.tostring())

if time:
    out.write(float_size.tostring())
    if gas:
        ar = array(f['PartType0']['TimeStep'][:],dtype='f') #The dtype is so that tostring doesn't give me float64
        out.write(ar.tostring())
    if halo:
        ar = array(f['PartType1']['TimeStep'][:],dtype='f')
        out.write(ar.tostring())
    if disk:
        ar = array(f['PartType2']['TimeStep'][:],dtype='f')
        out.write(ar.tostring())
    if bulge:
        ar = array(f['PartType3']['TimeStep'][:],dtype='f')
        out.write(ar.tostring())
    if star:
        ar = array(f['PartType4']['TimeStep'][:],dtype='f')
        out.write(ar.tostring())
    if bndry:
        ar = array(f['PartType5']['TimeStep'][:],dtype='f')
        out.write(ar.tostring())
    out.write(float_size.tostring())

f.close()
out.close()
