import numpy as np
import h5py, sys

#This program just loops over all of the massive halos in a sim and 
#calculates the radius out to which they are all uncontaminated
#We do this because:
#
#1) do we get any extra halos for free
#2) are we putting computing time into something that's low res

if len(sys.argv) < 4:
    print 'usage: python rock_check_hi_res_nearby.py <host mass> <rockstar halo file> <snapshot>'

f_rock = np.loadtxt(sys.argv[2])

h=0.71

id_rock = f_rock[:,0]
M_rock = f_rock[:,2]/h
Rvir_rock = f_rock[:,4]/h
Rmax_rock = f_rock[:,6]/h
Vmax_rock = f_rock[:,5]
centers_rock = np.zeros((len(f_rock[:,8]),3))
centers_rock[:,0], centers_rock[:,1], centers_rock[:,2] = f_rock[:,8]*1000.0/h, f_rock[:,9]*1000.0/h, f_rock[:,10]*1000.0/h
velocity_rock = np.zeros((len(f_rock[:,8]),3))
velocity_rock[:,0], velocity_rock[:,1], velocity_rock[:,2] = f_rock[:,11]/h, f_rock[:,12]/h, f_rock[:,13]/h

host_mass_diff = abs(float(sys.argv[1])-M_rock)
host_select = (host_mass_diff==min(host_mass_diff))
host_mass = M_rock[host_select][0]
host_center = centers_rock[host_select][0]
host_rvir = Rvir_rock[host_select][0]

mass_select = (M_rock>1.0e10)

id_rock_select = id_rock[mass_select]
mass_rock_select = M_rock[mass_select]
vmax_rock_select = Vmax_rock[mass_select]
rvir_rock_select = Rvir_rock[mass_select]
rmax_rock_select = Rmax_rock[mass_select]
center_rock_select = centers_rock[mass_select]

one_Rvir_list = []
two_Rvir_list = []
three_Rvir_list = []
contaminated_list = []
three_Rvir_mass_list = []

print 'loading particle data'

f_part = h5py.File(sys.argv[3])
halo_particle_pos = f_part['PartType1']['Coordinates'][:]/h
halo_masses = f_part['PartType1']['Masses'][:][0]

disk_particle_pos = f_part['PartType2']['Coordinates'][:]/h
#disk_masses = f_part['PartType2']['Masses'][:]

print 'running through halos\n'

print 'for host of mass: '+str(host_mass)+'\n'

for j in range(len(mass_rock_select)):
    center = center_rock_select[j]
    R_gal = rvir_rock_select[j]
    id_gal = id_rock_select[j]
    mass_gal = mass_rock_select[j]

    dist_from_host = np.sqrt((center[0]-host_center[0])**2.0+(center[1]-host_center[1])**2.0+(center[2]-host_center[2])**2.0)

    disk_dist = np.sqrt((center[0]-disk_particle_pos[:,0])**2.0+(center[1]-disk_particle_pos[:,1])**2.0+(center[2]-disk_particle_pos[:,2])**2.0)
    disk_dist_select = (disk_dist<R_gal)
    disk_dist_gal = disk_dist[disk_dist_select]
    halo_dist = np.sqrt((center[0]-halo_particle_pos[:,0])**2.0+(center[1]-halo_particle_pos[:,1])**2.0+(center[2]-halo_particle_pos[:,2])**2.0)
    halo_dist_select = (halo_dist<R_gal)
    halo_dist_gal = halo_dist[halo_dist_select]
    
    print 'id_gal: '+str(id_gal)
    print 'mass: '+'{0:.2e}'.format(mass_gal)
    print 'distance from host: '+str(dist_from_host)+' kpc, or '+str(dist_from_host/host_rvir)+' Rvir'
    print 'Uncontaminated out to '+str(min(disk_dist)/R_gal)+' Rvir'
    print 'fraction high res: '+str(float(len(halo_dist_gal))*halo_masses*10**10.0/(h*mass_gal))+'\n'
