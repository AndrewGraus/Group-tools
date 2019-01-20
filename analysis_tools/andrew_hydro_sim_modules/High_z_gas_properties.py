import h5py, os
import numpy as np
import numpy.ma as ma
import astropy.constants as const
import pickle

def High_z_gas_properties(snap_dir,snap_start,starting_snap_num='052'):
    '''
    The purpose of this program is to test the gas properties of the gas that forms stars
    in the reionization runs. The plan is to take all of the stars at z = 6 (snapshot 
    number 052 in the standard 600 FIRE snaps framework) and go back in time and find the 
    gas particles that those stars formed out of. Then track the properties of those gas 
    particle back to the beginning of the simulation.
    '''

    snap_loc = snap_dir+str(snap_start)+str(starting_snap_num)+'.hdf5'

    f_initial = h5py.File(snap_loc)
    star_ids = f_initial['PartType4']['ParticleIDs'][:]

    star_ids_sort = np.sort(star_ids)

    print 'stars at z = 6: '+str(len(star_ids))

    output_hdf5 = h5py.File('gas_particles.hdf5')

    #Now I want all the numbers from your snap to 000
    snapshot_numbers = ["%03d" % x for x in np.arange(1,int(starting_snap_num),1)]
    #This takes the number of snaps and assures they are in the 000 format of snaps

    for snap_id in snapshot_numbers:
        snap_loc_z = snap_loc.replace(str(starting_snap_num),str(snap_id))
        f_z = h5py.File(snap_loc_z)

        'at z= '+str(f_z['Header'].attrs['Redshift'])

        redshift = f_z['Header'].attrs['Redshift']

        gas_ids_at_z = f_z['PartType0']['ParticleIDs'][:]
        rho_at_z = f_z['PartType0']['Density'][:]
        EA_at_z = f_z['PartType0']['ElectronAbundance'][:]
        IE_at_z = f_z['PartType0']['InternalEnergy'][:]
        Z_at_z = f_z['PartType0']['Metallicity'][:,1] #Note that metals are an array so it doesn't really work with how the arrays are constructed
        NH_at_z = f_z['PartType0']['NeutralHydrogenAbundance'][:]

        #Now I want to sort these array BY the gas ids array
        
        sort_mask = np.argsort(gas_ids_at_z)
        gas_ids_sort = gas_ids_at_z[sort_mask]
        rho_sort = rho_at_z[sort_mask]
        EA_sort = EA_at_z[sort_mask]
        IE_sort = IE_at_z[sort_mask]
        Z_sort = Z_at_z[sort_mask]
        NH_sort = NH_at_z[sort_mask]
    
        star_ids_masked = ma.masked_where(np.in1d(star_ids_sort,gas_ids_sort,invert=True),star_ids_sort)
    
        ma.set_fill_value(star_ids_masked,0.0)
    
        gas_final = star_ids_masked.filled(0.0).astype('float64')
        rho_final = star_ids_masked.filled(0.0).astype('float64')
        EA_final = star_ids_masked.filled(0.0).astype('float64')
        IE_final = star_ids_masked.filled(0.0).astype('float64')
        Z_final = star_ids_masked.filled(0.0).astype('float64')
        NH_final = star_ids_masked.filled(0.0).astype('float64')

        gas_mask = np.in1d(gas_ids_sort,star_ids_sort) #the GAS ids which are in the star sample at this z
        rho_from_gas_mask = rho_sort[gas_mask].astype('float64') #The densities of gas particles which turn into stars by z =6
        EA_from_gas_mask = EA_sort[gas_mask].astype('float64')
        IE_from_gas_mask = IE_sort[gas_mask].astype('float64')
        Z_from_gas_mask = Z_sort[gas_mask].astype('float64')
        NH_from_gas_mask = NH_sort[gas_mask].astype('float64')

        print EA_from_gas_mask

        np.place(rho_final, rho_final!=0.0,rho_from_gas_mask)
        np.place(EA_final, EA_final!=0.0,EA_from_gas_mask)
        np.place(IE_final, IE_final!=0.0,IE_from_gas_mask)
        np.place(Z_final, Z_final!=0.0,Z_from_gas_mask)
        np.place(NH_final, NH_final!=0.0,NH_from_gas_mask)

        print EA_final
        print np.sum(EA_final!=0.0)
        print np.sum(IE_final!=0.0)
        print np.sum(Z_final!=0.0)
        print np.sum(NH_final!=0.0)

        #okay I think we just for every redshift make a group,
        #then make a subgroup at each redshift that is
        #the star ids then all of the gas properties

        group = output_hdf5.create_group(str(round(redshift,3)))
        dset_i = group.create_dataset('ids',data=star_ids_sort)
        dset_rho = group.create_dataset('rho',data=rho_final)
        dset_EA = group.create_dataset('ElectronAbundance',data=EA_final)
        dset_IE = group.create_dataset('InternalEnergy',data=IE_final)
        dset_Z = group.create_dataset('heliummassfraction',data=Z_final)
        dset_NH = group.create_dataset('NeutralHydrogenAbundance',data=NH_final)

        assert np.sum(rho_final!=0.0)==np.sum(np.in1d(star_ids_sort,gas_ids_sort))
        #Dont use Electron abundance, because some values are supposed to be zero
        #assert np.sum(EA_final!=0.0)==np.sum(np.in1d(star_ids_sort,gas_ids_sort))
        assert np.sum(IE_final!=0.0)==np.sum(np.in1d(star_ids_sort,gas_ids_sort))
        assert np.sum(Z_final!=0.0)==np.sum(np.in1d(star_ids_sort,gas_ids_sort))
        assert np.sum(NH_final!=0.0)==np.sum(np.in1d(star_ids_sort,gas_ids_sort))

    output_hdf5.close()

def Calculate_temperature(snapshot):
    #Okay I want to give a snapshot and then it autocalculates the temperatures
    
    proton_mass = const.m_p  #in kg by default
    gamma = 5.0/3.0
    k_B = const.k_B #in J/K by default

    pt0 = snapshot['PartType0']
    EA = pt0['ElectronAbundance'][:]
    IE = pt0['InternalEnergy'][:]
    rho = pt0['Density'][:]
    helium_mass_fracs = pt0['Metallicity'][:,1]
    ys_helium = helium_mass_fracs / (4.0 * (1.0 - helium_mass_fracs))
    mus = (1.0+4.0*ys_helium) / (1.0+ys_helium+EA)

    molecular_weights = proton_mass*mus
    
    Temp = (1.0e3)**2.0 * molecular_weights * (gamma-1.0) * IE / k_B #need to convert IE from (km/s)^2 to (m/s)^2

    return Temp

def return_gas_temp_and_coords(snap_dir,snap_start,starting_snap_num='001',ending_snap_num='052',file_name='gas_properties.hdf5'):
    '''
    This program prints the gas temperatures and coordinates 
    
    Inputs:
    snap_dir - the directory containing the snapshots
    snap_start - how the snapshot file names start (before the number, including the _)
    snap_num - the snapshot where it will end (52 is z = 6 for the standard FIRE naming scheme)
    file_name - what to name the output file
    '''

    output_hdf5 = h5py.File(file_name)
    snap_loc = snap_dir+str(snap_start)+str(starting_snap_num)+'.hdf5'

    #Now I want all the numbers from your snap to 000
    snapshot_numbers = ["%03d" % x for x in np.arange(int(starting_snap_num),int(ending_snap_num)+1,1)]
    #This takes the number of snaps and assures they are in the 000 format of snaps
    
    for snap_id in snapshot_numbers:
        print snap_id
        snap_loc_z = snap_loc.replace(str(starting_snap_num),str(snap_id)) #replace number with new number

        print 'loading snap'
        f_z = h5py.File(snap_loc_z)

        'at z= '+str(f_z['Header'].attrs['Redshift'])

        redshift = f_z['Header'].attrs['Redshift']

        gas_ids_at_z = f_z['PartType0']['ParticleIDs'][::10]
        rho_at_z = f_z['PartType0']['Density'][::10]
        coords_at_z = f_z['PartType0']['Coordinates'][::10]
        NH_at_z = f_z['PartType0']['NeutralHydrogenAbundance'][::10]

        print 'calculating temp'
        T_at_z = Calculate_temperature(f_z) #refer to temperature calculation function above

        print 'adding to output file'
        group = output_hdf5.create_group(str(round(redshift,3)))
        dset_i = group.create_dataset('ids',data=gas_ids_at_z)
        dset_rho = group.create_dataset('rho',data=rho_at_z)
        dset_NH = group.create_dataset('NeutralHydrogenAbundance',data=NH_at_z)
        dset_T = group.create_dataset('T',data=T_at_z)

    output_hdf5.close()

def return_gas_temp_and_coords_struct(snap_dir,snap_start,starting_snap_num='001',ending_snap_num='052',file_name='gas_properties'):
    '''
    This program prints the gas temperatures and coordinates 

    This module version changes the output file to a structured array (much smaller than hdf5)
    
    Inputs:
    snap_dir - the directory containing the snapshots
    snap_start - how the snapshot file names start (before the number, including the _)
    snap_num - the snapshot where it will end (52 is z = 6 for the standard FIRE naming scheme)
    file_name - what to name the output file (will append .npy)
    '''

    snap_dict = {}

    snap_loc = snap_dir+str(snap_start)+str(starting_snap_num)+'.hdf5'

    #Now I want all the numbers from your snap to 000
    snapshot_numbers = ["%03d" % x for x in np.arange(int(starting_snap_num),int(ending_snap_num)+1,1)]
    #This takes the number of snaps and assures they are in the 000 format of snaps
    
    for snap_id in snapshot_numbers:
        print snap_id
        snap_loc_z = snap_loc.replace(str(starting_snap_num),str(snap_id)) #replace number with new number

        print 'loading snap'
        f_z = h5py.File(snap_loc_z)

        'at z= '+str(f_z['Header'].attrs['Redshift'])

        redshift = f_z['Header'].attrs['Redshift']

        gas_ids_at_z = f_z['PartType0']['ParticleIDs'][::10]
        rho_at_z = f_z['PartType0']['Density'][::10]
        coords_at_z = f_z['PartType0']['Coordinates'][::10]
        NH_at_z = f_z['PartType0']['NeutralHydrogenAbundance'][::10]

        print 'calculating temp'
        T_at_z = Calculate_temperature(f_z) #refer to temperature calculation function above

        print 'adding to output file'
        #group = output_hdf5.create_group(str(round(redshift,3)))
        
        snap_dict[str(round(redshift,3))] = {'ids':None,'rho':None,'NeutralHydrogenAbundance':None,'T':None}

        snap_dict[str(round(redshift,3))]['ids'] = gas_ids_at_z
        snap_dict[str(round(redshift,3))]['rho'] = rho_at_z
        snap_dict[str(round(redshift,3))]['NeutralHydrogenAbundance'] = NH_at_z
        snap_dict[str(round(redshift,3))]['T'] = T_at_z

        #dset_i = group.create_dataset('ids',data=gas_ids_at_z)
        #dset_rho = group.create_dataset('rho',data=rho_at_z)
        #dset_NH = group.create_dataset('NeutralHydrogenAbundance',data=NH_at_z)
        #dset_T = group.create_dataset('T',data=T_at_z)

    #save with pickle
    with open('temp_test.pkl', 'wb') as f:
        pickle.dump(snap_dict,f,pickle.HIGHEST_PROTOCOL) #HIGHEST_PROTOCOL saves it as a binary

    #save with numpy
    np.save('temp_test',snap_dict)
        
    

