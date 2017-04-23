import numpy as np
import h5py, sys

#From rockstar
#
#This is a simple program that takes in the rockstar halo file (ascii format), and 
#grabs a halo closest in mass to a target mass, and then takes the snapshot (hdf5)
#and calculates how far away low res particles are
if len(sys.argv) < 4:
    print 'usage: python rock_verification.py <halo mass> <location of rock halo file> <location of hdf5 file>'

f_rock = np.loadtxt(sys.argv[2])

h=0.702

#f_hires = f_rock[:,37]
M_rock = f_rock[:,2]/h

Rvir_rock = f_rock[:,4]/h
Rmax_rock = f_rock[:,6]/h
Vmax_rock = f_rock[:,5]
X_rock, Y_rock, Z_rock = f_rock[:,8]*1000.0/h,f_rock[:,9]*1000.0/h,f_rock[:,10]*1000.0/h
Vx_rock, Vy_rock, Vz_rock = f_rock[:,11],f_rock[:,12],f_rock[:,13]

print Rvir_rock[0]
print X_rock[0], Y_rock[0], Z_rock[0]

mass_diff = abs(M_rock-float(sys.argv[1]))

mass_select = (mass_diff==min(mass_diff))

mass_rock_select = M_rock[mass_select]
#f_hires_select = f_hires[mass_select]
print len(mass_rock_select)
print 'Vmax = '+str(Vmax_rock[mass_select][0])+' km/s'
print 'mass = '+'{0:.2e}'.format(mass_rock_select[0])+' Msun'
print 'Rvir = '+str(Rvir_rock[mass_select][0])+' kpc'
print 'Rmax = '+str(Rmax_rock[mass_select][0])+' kpc'
print 'position = '+str(X_rock[mass_select][0])+', '+str(Y_rock[mass_select][0])+', '+str(Z_rock[mass_select][0])+' kpc'
print 'velocity = '+str(Vx_rock[mass_select][0])+', '+str(Vy_rock[mass_select][0])+', '+str(Vz_rock[mass_select][0])+' km/s'
#print 'high res fraction = '+str(f_hires_select)+'\n'

print mass_rock_select

center = X_rock[mass_select][0], Y_rock[mass_select][0], Z_rock[mass_select][0]
R_gal = Rvir_rock[mass_select]

f = h5py.File(sys.argv[3])
#f['PartType4'].keys()
#star_coords = f['PartType4']['Coordinates'][:]
#star_mass = f['PartType4']['Masses'][:]
#star_age = f['PartType4']['StellarFormationTime'][:] #units are a
#star_dist = np.sqrt((center[0]-star_coords[:,0])**2.0+(center[1]-star_coords[:,1])**2.0+(center[2]-star_coords[:,2])**2.0)
#star_dist_select = (star_dist<R_gal)
#star_dist_gal = star_dist[star_dist_select]

halo_coords = f['PartType1']['Coordinates'][:]/h
halo_masses = f['PartType1']['Masses'][:]
halo_dist = np.sqrt((center[0]-halo_coords[:,0])**2.0+(center[1]-halo_coords[:,1])**2.0+(center[2]-halo_coords[:,2])**2.0)
halo_dist_select = (halo_dist<R_gal)
halo_dist_gal = halo_dist[halo_dist_select]
print min(halo_dist), max(halo_dist)

disk_coords = f['PartType2']['Coordinates'][:]/h
disk_masses = f['PartType2']['Masses'][:]
disk_dist = np.sqrt((center[0]-disk_coords[:,0])**2.0+(center[1]-disk_coords[:,1])**2.0+(center[2]-disk_coords[:,2])**2.0)
disk_dist_select = (disk_dist<R_gal)
disk_dist_gal = disk_dist[disk_dist_select]
print min(disk_dist), max(disk_dist)

print 'distance to closest low res particle n*rvir '+str(min(disk_dist)/R_gal)
print 'number of disk particles in Rvir '+str(len(disk_dist_gal))
print 'hi-res particles in Rvir '+str(len(halo_dist_gal))
#print 'star particles in Rvir '+str(len(star_dist_gal))
