import numpy as np
import h5py
import astropy.constants as const

def Calculate_temperature(snapshot,data_slice=1):
    #Okay I want to give a snapshot and then it autocalculates the temperatures
    
    proton_mass = const.m_p  #in kg by default
    gamma = 5.0/3.0
    k_B = const.k_B #in J/K by default

    pt0 = snapshot['PartType0']
    EA = pt0['ElectronAbundance'][::data_slice]
    IE = pt0['InternalEnergy'][::data_slice]
    rho = pt0['Density'][::data_slice]
    helium_mass_fracs = pt0['Metallicity'][:,1][::data_slice]
    ys_helium = helium_mass_fracs / (4.0 * (1.0 - helium_mass_fracs))
    mus = (1.0+4.0*ys_helium) / (1.0+ys_helium+EA)

    molecular_weights = proton_mass*mus
    
    Temp = (1.0e3)**2.0 * molecular_weights * (gamma-1.0) * IE / k_B #need to convert IE from (km/s)^2 to (m/s)^2

    return Temp

def return_gas_temp_and_coords(snap_dir,snap_start,starting_snap_num='052',file_name='gas_properties.hdf5',data_slice=1):
    '''
    The purpose of this program is to test the gas properties of the gas that forms stars
    in the reionization runs. The plan is to take all of the stars at z = 6 (snapshot 
    number 052 in the standard 600 FIRE snaps framework) and go back in time and find the 
    gas particles that those stars formed out of. Then track the properties of those gas 
    particle back to the beginning of the simulation.
    '''

    output_hdf5 = h5py.File(file_name)
    snap_loc = snap_dir+str(snap_start)+str(starting_snap_num)+'.hdf5'

    #Now I want all the numbers from your snap to 000
    snapshot_numbers = ["%03d" % x for x in np.arange(1,int(starting_snap_num)+1,1)]
    #This takes the number of snaps and assures they are in the 000 format of snaps
    
    for snap_id in snapshot_numbers:
        print snap_id
        snap_loc_z = snap_loc.replace(str(starting_snap_num),str(snap_id))

        print 'loading snap'
        f_z = h5py.File(snap_loc_z)

        'at z= '+str(f_z['Header'].attrs['Redshift'])

        redshift = f_z['Header'].attrs['Redshift']

        gas_ids_at_z = f_z['PartType0']['ParticleIDs'][::data_slice]
        
        print 'total number of particles: '+str(len(f_z['PartType0']['ParticleIDs'][:]))
        print 'particles written out: '+str(len(gas_ids_at_z))

        rho_at_z = f_z['PartType0']['Density'][::data_slice]
        coords_at_z = f_z['PartType0']['Coordinates'][::data_slice]
        NH_at_z = f_z['PartType0']['NeutralHydrogenAbundance'][::data_slice]

        print 'calculating temp'
        T_at_z = Calculate_temperature(f_z,data_slice=data_slice)

        print 'adding to output file'
        group = output_hdf5.create_group(str(round(redshift,3)))
        dset_i = group.create_dataset('ids',data=gas_ids_at_z)
        dset_rho = group.create_dataset('rho',data=rho_at_z)
        dset_NH = group.create_dataset('NeutralHydrogenAbundance',data=NH_at_z)
        dset_T = group.create_dataset('T',data=T_at_z)
        dset_coords = group.create_dataset('coords',data=coords_at_z)

    output_hdf5.close()

snap_dir = './m10b/'
snap_start = 'snapshot_m10b_fg_z_6_'
data_slice=100

return_gas_temp_and_coords(snap_dir,snap_start,file_name='gas_properties_100.hdf5',data_slice=data_slice)
