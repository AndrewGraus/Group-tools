#The idea behind this is to make some simple functions that I can use to easily
#analyze a hydro sim.  given the Z = 0 snap and halo finder file for the dark
#matter only and full hydro sims

#Problem, going to need to figure out how to handle identifying the host halo

def Load_Particle_Data(giz_hdf5,add_dm=True,add_gas=False,add_stars=False,add_low_res=False):
    #This is just going to load particles from a snapshot into a custom dictionary
    #This works almost exactly like the gizmo hdf5s but I do this to save memory
    #if I only need certain particle species, and it also converts units
    #and calculates [Fe/H]

    #eventually I'll modify this to load in different redshifts as well

    import numpy as np
    import yt, h5py, re, os
    from math import log10
    import astropy.constants as const
    from astropy.cosmology import FlatLambdaCDM

    f = h5py.File(giz_hdf5)
    h = f['Header'].attrs['HubbleParam'] #grab h
    z = f['Header'].attrs['Redshift'] #grab z

    PD_dict = {}
    if add_dm==True:
        #add high res dm particles
        halo_coords = f['PartType1']['Coordinates'][:]/h
        halo_masses = f['PartType1']['Masses'][:]*10**10.0/h
        halo_vels = f['PartType1']['Velocities'][:]
    
        PD_dict['halo'] = {'coords':halo_coords,'masses':halo_masses,'velocities':halo_vels}

    if add_gas==True:
        #add gas particles
        gas_part_coords = f['PartType0']['Coordinates'][:]/h
        gas_part_masses = f['PartType0']['Masses'][:]*10**10.0/h
        gas_vels = f['PartType0']['Velocities'][:]
        gas_metal = f['PartType0']['Metallicity'][:]
        gas_IE = f['PartType0']['InternalEnergy'][:]
        gas_EA = f['PartType0']['ElectronAbundance'][:]
        gas_He= gas_metal[:,1] #number fraction
        gas_y_He = gas_He / (4.0 * (1.0 - gas_He))
        Molecular_weights = (1.0 + 4.0*gas_y_He)/(1.0+gas_y_He+gas_EA)*const.m_p.value
        T_gas = (1000.0)**2.0*(5.0/3.0-1.0)*Molecular_weights*gas_IE / const.k_B.value #This should give temp in K Unless I screwed something up (I probably screwed something up)

        PD_dict['gas'] = {'coords':gas_part_coords,'masses':gas_part_masses,
                          'velocities':gas_vels,'temperature':T_gas}

    if add_low_res==True:
        #add disk particles
        disk_coords = f['PartType2']['Coordinates'][:]/h
        disk_masses = f['PartType2']['Masses'][:]*10**10.0/h
        disk_vels = f['PartType2']['Velocities'][:]

        PD_dict['disk'] = {'coords':disk_coords,'masses':disk_masses,'velocities':disk_vels}

    if add_stars==True:
        #add stars and calculate Fe/H
        star_part_coord = f['PartType4']['Coordinates'][:]/h
        star_part_age = f['PartType4']['StellarFormationTime'][:]
        star_part_metal = f['PartType4']['Metallicity'][:,0] #Z
        star_part_Fe = f['PartType4']['Metallicity'][:,10]
        star_part_He = f['PartType4']['Metallicity'][:,1] #Y
        star_part_masses = f['PartType4']['Masses'][:]*10**10.0/h
        star_part_H = [1.0-star_part_metal[xx]-star_part_He[xx] for xx in range(len(star_part_metal))]
        Fe_H_sun = log10(4.0e-5)
        Fe_H = np.asarray([log10(star_part_Fe[xx]*1.00794/(star_part_H[xx]*55.845)) - Fe_H_sun for xx in range(len(star_part_H))])
        star_vels = f['PartType4']['Velocities'][:]

        PD_dict['star'] = {'coords':star_part_coord,'masses':star_part_masses,'ages':star_part_age,'metallicity':star_part_metal,'Fe_H':Fe_H,'velocities':star_vels}

    return PD_dict


def Load_Halo_Data(giz_halo_file,h=0.702):
    import numpy as np
    import yt, h5py, re, os
    from astropy.cosmology import FlatLambdaCDM

    #This loads the halo data from a rockstar catalog into a dictionary with properties that I want
    #Ill probably need to modify this to get more things or to use AHF instead eventually

    #Note h is NOT defined I could either find some way to parse it from the rockstar header
    #or require the snapshot as well (not ideal)

    
    f_rock = np.loadtxt(giz_halo_file)

    assert len(f_rock[0]) == 53

    id_rock = f_rock[:,0]
    M_rock = f_rock[:,2]/h
    Rvir_rock = f_rock[:,4]/h
    Rmax_rock = f_rock[:,6]/h
    Vmax_rock = f_rock[:,5]
    centers_rock = np.zeros((len(f_rock[:,8]),3))
    centers_rock[:,0], centers_rock[:,1], centers_rock[:,2] = f_rock[:,8]*1000.0/h, f_rock[:,9]*1000.0/h, f_rock[:,10]*1000.0/h #convert to kpc
    velocity_rock = np.zeros((len(f_rock[:,8]),3))
    velocity_rock[:,0], velocity_rock[:,1], velocity_rock[:,2] = f_rock[:,11], f_rock[:,12], f_rock[:,13]

    return {'ids':id_rock, 'masses':M_rock, 'rvir':Rvir_rock, 'rmax':Rmax_rock, 'vmax':Vmax_rock, 'centers':centers_rock,'velocities':velocity_rock}

def Load_Halo_Data_AHF(giz_halo_file,h=0.702):
    import numpy as np
    import yt, h5py, re, os
    from astropy.cosmology import FlatLambdaCDM

    #This loads the halo data from a rockstar catalog into a dictionary with properties that I want
    #Ill probably need to modify this to get more things or to use AHF instead eventually

    #Note h is NOT defined I could either find some way to parse it from the rockstar header
    #or require the snapshot as well (not ideal)
    
    f_ahf = np.loadtxt(giz_halo_file)

    assert len(f_ahf[0]) == 83

    id_ahf = f_ahf[:,0]
    M_ahf = f_ahf[:,3]/h
    Rvir_ahf = f_ahf[:,11]/h
    Rmax_ahf = f_ahf[:,12]/h
    Vmax_ahf = f_ahf[:,16]
    centers_ahf = np.zeros((len(f_ahf[:,8]),3))
    centers_ahf[:,0], centers_ahf[:,1], centers_ahf[:,2] = f_ahf[:,5]/h, f_ahf[:,6]/h, f_ahf[:,7]/h #convert to kpc (for AHF centers are already kpc)
    velocity_ahf = np.zeros((len(f_ahf[:,8]),3))
    velocity_ahf[:,0], velocity_ahf[:,1], velocity_ahf[:,2] = f_ahf[:,8], f_ahf[:,9], f_ahf[:,10]
    f_hires = f_ahf[:,37]

    return {'ids':id_ahf, 'masses':M_ahf, 'rvir':Rvir_ahf, 'rmax':Rmax_ahf, 'vmax':Vmax_ahf, 'centers':centers_ahf,'velocities':velocity_ahf,'f_hires':f_hires}

def Identify_Host(giz_hdf5,halo_file,add_velocity=False,print_values=False,print_hi_res_halos=False,subhalo_limit=10**10.0):
    import numpy as np
    import yt, h5py, re, os
    from astropy.cosmology import FlatLambdaCDM

    #This script identifies all the halos above 1.0e10 (that's hard coded) and calculates all of the things
    #that are hi-res and then returns the center of the highest mass halo that doesn't have a low res particle
    #within Rvir

    #NOTE: the plan is to make this a legacy function because I need it to do more stuff
    #but I'm leaving it in here because I might have used it in some previous program
    #and if I delete this / modify it it may break those
    #
    #There's nothing wrong with it, but I want to have it return more values without
    #having to add in a bunch of additional options

    print 'loading data'

    halo_dict = Load_Halo_Data(halo_file) #load in the rockstar file
    PD_dict = Load_Particle_Data(giz_hdf5,add_low_res=True) #load in the particle data'

    id_rock = halo_dict['ids'] 
    M_rock = halo_dict['masses']
    Vmax_rock = halo_dict['vmax']
    Rvir_rock = halo_dict['rvir']
    Rmax_rock = halo_dict['rmax']
    centers_rock = halo_dict['centers']
    velocities_rock = halo_dict['velocities']

    mass_select = (M_rock>subhalo_limit) #This is hard coded in

    id_rock_select = id_rock[mass_select]
    mass_rock_select = M_rock[mass_select]
    vmax_rock_select = Vmax_rock[mass_select]
    rvir_rock_select = Rvir_rock[mass_select]
    rmax_rock_select = Rmax_rock[mass_select]
    center_rock_select = centers_rock[mass_select]
    velocities_rock_select = velocities_rock[mass_select]

    particle_coords = PD_dict['halo']['coords'] #grab hi and low res particles
    low_res_coords = PD_dict['disk']['coords']

    M_hi_res = 0.0
    total_hi_res = []
    closest_low_res_list = []

    print 'looping over '+str(len(mass_rock_select))+' halos.'

    for j in range(len(mass_rock_select)):
        center = center_rock_select[j]
        R_gal = rvir_rock_select[j]
        id_gal = id_rock_select[j]
        mass_gal = mass_rock_select[j]
        vel_gal = velocities_rock_select[j]
        vmax_gal = vmax_rock_select[j]

        low_res_dist = np.sqrt((center[0]-low_res_coords[:,0])**2.0+(center[1]-low_res_coords[:,1])**2.0+(center[2]-low_res_coords[:,2])**2.0)
        if min(low_res_dist) > R_gal:
            #count the total number of high res halos above 1e10
            total_hi_res.append(id_gal)
            closest_low_res_list.append(min(low_res_dist)/R_gal)
            if mass_gal > M_hi_res:
                #if this is the most massive hi res halo record some stats
                M_hi_res = mass_gal
                id_hi_res = id_gal
                rvir_hi_res = R_gal
                center_hi_res = center
                closest_low_res = min(low_res_dist)/R_gal
                vel_hi_res = vel_gal
                vmax_gal_hi_res = vmax_gal
    if print_values==True:
        print '\n '
        print 'total number of large galaxies: {0:d}'.format(len(mass_rock_select))
        print 'total number of high res galaxies: {0:d}'.format(len(total_hi_res))
        print '\n'
        print '{:#^30} \n'.format('HOST PROPERTIES')
        print 'id_gal: '+str(id_hi_res)
        print 'mass: '+'{0:.2e} Msun'.format(M_hi_res)
        print 'Vmax: '+'{0:.2e} kms-1'.format(vmax_gal_hi_res)
        print 'radius: '+'{0:.2e} kpc'.format(rvir_hi_res)
        print 'closest low res particle: {0:.2e} Rvir '.format(closest_low_res)

    #return 'hosts' center and Rvir

    if print_hi_res_halos==True:
        return total_hi_res, closest_low_res_list
    if add_velocity == False:
        return center_hi_res, rvir_hi_res
    else:
        return center_hi_res, rvir_hi_res, vel_hi_res

def Identify_Host_and_Subhalos(giz_hdf5,halo_file,print_values=True,subhalo_limit=10*10.0):
    import numpy as np
    import yt, h5py, re, os
    from astropy.cosmology import FlatLambdaCDM

    #This script identifies all the halos above 1.0e10 (that's hard coded) and calculates all of the things
    #that are hi-res and then returns the center of the highest mass halo that doesn't have a low res particle
    #within Rvir

    halo_dict = Load_Halo_Data(halo_file) #load in the rockstar file
    PD_dict = Load_Particle_Data(giz_hdf5,add_low_res=True) #load in the particle data

    id_rock = halo_dict['ids'] 
    M_rock = halo_dict['masses']
    Vmax_rock = halo_dict['vmax']
    Rvir_rock = halo_dict['rvir']
    Rmax_rock = halo_dict['rmax']
    centers_rock = halo_dict['centers']
    velocities_rock = halo_dict['velocities']

    mass_select = (M_rock>subhalo_limit) #This is hard coded in

    id_rock_select = id_rock[mass_select]
    mass_rock_select = M_rock[mass_select]
    vmax_rock_select = Vmax_rock[mass_select]
    rvir_rock_select = Rvir_rock[mass_select]
    rmax_rock_select = Rmax_rock[mass_select]
    center_rock_select = centers_rock[mass_select]
    velocities_rock_select = velocities_rock[mass_select]

    particle_coords = PD_dict['halo']['coords'] #grab hi and low res particles
    low_res_coords = PD_dict['disk']['coords']

    M_hi_res = 0.0
    total_hi_res = []
    closest_low_res_list = []
    hi_res_masses = []
    hi_res_X, hi_res_Y, hi_res_Z = [],[],[]
    hi_res_rvir = []
    hi_res_rmax = []
    hi_res_vmax = []

    for j in range(len(mass_rock_select)):
        center = center_rock_select[j]
        R_gal = rvir_rock_select[j]
        id_gal = id_rock_select[j]
        mass_gal = mass_rock_select[j]
        vel_gal = velocities_rock_select[j]
        vmax_gal = vmax_rock_select[j]
        rmax_gal = rmax_rock_select[j]

        low_res_dist = np.sqrt((center[0]-low_res_coords[:,0])**2.0+(center[1]-low_res_coords[:,1])**2.0+(center[2]-low_res_coords[:,2])**2.0)
        if min(low_res_dist) > R_gal:
            #count the total number of high res halos above 1e10
            total_hi_res.append(id_gal)
            hi_res_masses.append(mass_gal)
            hi_res_rvir.append(R_gal)
            hi_res_rmax.append(rmax_gal)
            hi_res_vmax.append(vmax_gal)
            hi_res_X.append(center[0])
            hi_res_Y.append(center[1])
            hi_res_Z.append(center[2])
            closest_low_res_list.append(min(low_res_dist)/R_gal)
            if mass_gal > M_hi_res:
                #if this is the most massive hi res halo record some stats
                M_hi_res = mass_gal
                id_hi_res = id_gal
                rvir_hi_res = R_gal
                center_hi_res = center
                closest_low_res = min(low_res_dist)/R_gal
                vel_hi_res = vel_gal
    if print_values==True:
        print '\n '
        print 'total number of large galaxies: {0:d}'.format(len(mass_rock_select))
        print 'total number of high res galaxies: {0:d}'.format(len(total_hi_res))
        print '\n'
        print '{:#^30} \n'.format('HOST PROPERTIES')
        print 'id_gal: '+str(id_hi_res)
        print 'mass: '+'{0:.2e} Msun'.format(M_hi_res)
        print 'radius: '+'{0:.2e} kpc'.format(rvir_hi_res)
        print 'center: '+str(center_hi_res)
        print 'halo_velocity: '+str(vel_hi_res)
        print 'closest low res particle: {0:.2e} Rvir '.format(int(closest_low_res))

    sub_distances = [np.sqrt((hi_res_X[xx]-center_hi_res[0])**2.0+(hi_res_Y[xx]-center_hi_res[1])**2.0+(hi_res_Z[xx]-center_hi_res[2])**2.0) for xx in range(len(hi_res_X))]
    
    halos_matrix = np.zeros((len(hi_res_X),13))
    halos_matrix[:,0] = total_hi_res
    halos_matrix[:,1] = hi_res_masses
    halos_matrix[:,2] = hi_res_rvir
    halos_matrix[:,3] = hi_res_rmax
    halos_matrix[:,4] = hi_res_vmax
    halos_matrix[:,5] = sub_distances
    halos_matrix[:,6] = hi_res_X
    halos_matrix[:,7] = hi_res_Y
    halos_matrix[:,8] = hi_res_Z

    return halos_matrix

def Load_hi_res_halo_dict(giz_hdf5,halo_file,Low_res_tolerance=1.0,save_particles=False,halo_lim=1.0e9):
    #The idea behind this module is to take the particles from the file
    #and identify the hi res halos up to some tolerance, which is by
    #default 100% of the particles have to be hi_res
    #
    #Then save the halos and their properties to a dictionary
    #
    #Eventually I will also add in the capability to save all the 
    #particles in each halo, and then maybe the positions of the
    #closest low res particles if there are mutliple within Rvir
    
    #okay, this takes up an absurd amount of memory
    #I'm going to make another module that does the same
    #thing but uses less memory

    import numpy as np

    halo_dict = Load_Halo_Data(halo_file)
    PD_dict = Load_Particle_Data(giz_hdf5,add_low_res=True,add_stars=True,add_gas=True)

    id_rock = halo_dict['ids'] 
    M_rock = halo_dict['masses']
    Vmax_rock = halo_dict['vmax']
    Rvir_rock = halo_dict['rvir']
    Rmax_rock = halo_dict['rmax']
    centers_rock = halo_dict['centers']
    velocities_rock = halo_dict['velocities']

    mass_select = (M_rock>halo_lim)
    
    id_rock_select = id_rock[mass_select]
    mass_rock_select = M_rock[mass_select]
    vmax_rock_select = Vmax_rock[mass_select]
    rvir_rock_select = Rvir_rock[mass_select]
    rmax_rock_select = Rmax_rock[mass_select]
    center_rock_select = centers_rock[mass_select]
    velocities_rock_select = velocities_rock[mass_select]

    rgal_rock_select = rvir_rock_select*0.1 #galaxy radius defined as just 10% Rvir
    
    hi_res_coords = PD_dict['halo']['coords']
    hi_res_mass = PD_dict['halo']['masses'][0]
    low_res_coords = PD_dict['disk']['coords']
    low_res_mass = PD_dict['disk']['masses'][0]

    star_coords = PD_dict['star']['coords']
    star_mass = PD_dict['star']['masses']

    gas_coords = PD_dict['gas']['coords']
    gas_mass = PD_dict['gas']['masses']
    gas_temp = PD_dict['gas']['temperature']

    #Now I want to find the distance between every halo and every particle
    #I can do this "easily" by using masked arrays
    #
    #The first step is to create a 3d array I can do this by taking a
    #2d array of the halo centers and doing center_rock_select[:,np.newaxis]
    #which makes a 3d array where each array of axis=0 is the halo center
    
    halo_particle_diff = center_rock_select[:,np.newaxis] - low_res_coords
    halo_hi_res_diff = center_rock_select[:,np.newaxis] - hi_res_coords
    star_diff = center_rock_select[:,np.newaxis] - star_coords
    gas_diff = center_rock_select[:,np.newaxis] - gas_coords

    #star and gas multiples
    star_mass_multiple = np.repeat(star_mass[np.newaxis,:],len(id_rock_select),axis=0)
    
    gas_mass_multiple = np.repeat(gas_mass[np.newaxis,:],len(id_rock_select),axis=0)

    #now I need to calculate the distance which can be accomplished by
    #linalg.norm along the 2nd axis which is the "row" of each
    #element in the 3d array

    low_res_dist_all = np.linalg.norm(halo_particle_diff,axis=2)
    hi_res_dist_all = np.linalg.norm(halo_hi_res_diff,axis=2)
    star_dist_all = np.linalg.norm(star_diff,axis=2)
    gas_dist_all = np.linalg.norm(gas_diff,axis=2)

    #now I want to do two things, calculate the lowest distance in each 
    #row of this array, and also mask this and calculate the total number
    #of particles in Rvir
    #
    #maybe just mask it and do the fraction calculation in a loop?

    #low_res_mask = (low_res_dist_all<rvir_rock_select[np.newaxis].T)
    low_res_masked_distances = np.ma.array(low_res_dist_all,mask=np.logical_not((low_res_dist_all<rvir_rock_select[np.newaxis].T)))
    hi_res_masked_distances = np.ma.array(hi_res_dist_all,mask=np.logical_not((hi_res_dist_all<rvir_rock_select[np.newaxis].T)))
    star_masked_distances = np.ma.array(star_dist_all,mask=np.logical_not((star_dist_all<rvir_rock_select[np.newaxis].T)))
    star_masked_distances_gal = np.ma.array(star_dist_all,mask=np.logical_not((star_dist_all<rgal_rock_select[np.newaxis].T)))
    star_masked_masses = np.ma.array(star_mass_multiple,mask=np.logical_not((star_dist_all<rvir_rock_select[np.newaxis].T)))
    star_masked_masses_gal = np.ma.array(star_mass_multiple,mask=np.logical_not((star_dist_all<rgal_rock_select[np.newaxis].T)))
    gas_masked_distances = np.ma.array(gas_dist_all,mask=np.logical_not((gas_dist_all<rvir_rock_select[np.newaxis].T)))
    gas_masked_distances_gal = np.ma.array(gas_dist_all,mask=np.logical_not((gas_dist_all<rgal_rock_select[np.newaxis].T)))
    gas_masked_masses = np.ma.array(gas_mass_multiple,mask=np.logical_not((gas_dist_all<rvir_rock_select[np.newaxis].T)))
    gas_masked_masses_gal = np.ma.array(gas_mass_multiple,mask=np.logical_not((gas_dist_all<rgal_rock_select[np.newaxis].T)))

    cold_gas_masked_masses = np.ma.array(gas_mass_multiple,mask=np.logical_not((gas_dist_all<rvir_rock_select[np.newaxis].T)&(gas_temp<np.ones_like(rvir_rock_select[np.newaxis]).T*1.0e4)))
    cold_gas_masked_masses_gal = np.ma.array(gas_mass_multiple,mask=np.logical_not((gas_dist_all<rgal_rock_select[np.newaxis].T)&(gas_temp<np.ones_like(rvir_rock_select[np.newaxis]).T*1.0e4)))

    #I can sum across each row to get the total mass
    #I can also maybe do np.min along each row to get the closest thing

    low_res_n_in_rvir = np.sum((low_res_dist_all<rvir_rock_select[np.newaxis].T),axis=1)
    hi_res_n_in_rvir = np.sum((hi_res_dist_all<rvir_rock_select[np.newaxis].T),axis=1)

    assert len(low_res_n_in_rvir)==len(id_rock_select)

    low_res_mass_in_rvir = low_res_n_in_rvir*low_res_mass
    hi_res_mass_in_rvir = hi_res_n_in_rvir*hi_res_mass

    star_mass_in_rvir = np.sum(star_masked_masses,axis=1)
    star_mass_in_rgal = np.sum(star_masked_masses_gal,axis=1)

    gas_mass_in_rvir = np.sum(gas_masked_masses,axis=1)
    gas_mass_in_rgal = np.sum(gas_masked_masses_gal,axis=1)    

    temp_mask = np.ones_like(gas_masked_masses)*1.0e4
    temp_mask_gal = np.ones_like(gas_masked_masses_gal)*1.0e4

    cold_gas_mass_in_rvir = np.sum(cold_gas_masked_masses,axis=1)
    cold_gas_mass_in_rgal = np.sum(cold_gas_masked_masses_gal,axis=1)

    hi_res_fraction = np.divide(hi_res_mass_in_rvir,np.add(hi_res_mass_in_rvir,low_res_mass_in_rvir))

    assert len(hi_res_fraction)==len(id_rock_select)

    hi_res_frac_mask = (hi_res_fraction>=Low_res_tolerance)

    print len(id_rock_select[hi_res_frac_mask])
    
    print 'fraction of hi res particles:'
    print hi_res_fraction

    print 'galaxy stellar masses:'
    print star_mass_in_rgal

    print 'galaxy gas masses:'
    print gas_mass_in_rgal

    #now I want to make a dictionary that has each halo and its properties in it

    halo_dict = {'halo_ids':id_rock_select[hi_res_frac_mask],'M_vir':mass_rock_select[hi_res_frac_mask],'V_max':vmax_rock_select[hi_res_frac_mask],'R_vir':rvir_rock_select[hi_res_frac_mask],'centers':center_rock_select[hi_res_frac_mask],'hi_res_fraction':hi_res_fraction[hi_res_frac_mask],'M_star_gal':np.ma.getdata(star_mass_in_rgal[hi_res_frac_mask]),'M_star_vir':np.ma.getdata(star_mass_in_rvir[hi_res_frac_mask]),'M_gas_rvir':np.ma.getdata(gas_mass_in_rvir[hi_res_frac_mask]),'M_gas_cold_rvir':np.ma.getdata(cold_gas_mass_in_rvir[hi_res_frac_mask])}

    return halo_dict
    
def Load_hi_res_halo_dict_mem_save(giz_hdf5,halo_file,Low_res_tolerance=1.0,save_particles=False,halo_lim=1.0e9):
    #This is the same program as above, but reorganized to hopefully
    #save enough memory to be practical    

    import numpy as np

    halo_dict = Load_Halo_Data(halo_file)
    PD_dict = Load_Particle_Data(giz_hdf5,add_low_res=True,add_stars=True,add_gas=True)

    id_rock = halo_dict['ids'] 
    M_rock = halo_dict['masses']
    Vmax_rock = halo_dict['vmax']
    Rvir_rock = halo_dict['rvir']
    Rmax_rock = halo_dict['rmax']
    centers_rock = halo_dict['centers']
    velocities_rock = halo_dict['velocities']

    mass_select = (M_rock>halo_lim)
    
    id_rock_select = id_rock[mass_select]
    mass_rock_select = M_rock[mass_select]
    vmax_rock_select = Vmax_rock[mass_select]
    rvir_rock_select = Rvir_rock[mass_select]
    rmax_rock_select = Rmax_rock[mass_select]
    center_rock_select = centers_rock[mass_select]
    velocities_rock_select = velocities_rock[mass_select]

    rgal_rock_select = rvir_rock_select*0.1 #galaxy radius defined as just 10% Rvir
    
    #hi_res_coords = PD_dict['halo']['coords']
    #hi_res_mass = PD_dict['halo']['masses'][0]
    low_res_coords = PD_dict['disk']['coords']
    low_res_mass = PD_dict['disk']['masses'][0]

    star_coords = PD_dict['star']['coords']
    star_mass = PD_dict['star']['masses']

    gas_coords = PD_dict['gas']['coords']
    gas_mass = PD_dict['gas']['masses']
    gas_temp = PD_dict['gas']['temperature']

    #Now I want to find the distance between every halo and every particle
    #I can do this "easily" by using masked arrays
    #
    #The first step is to create a 3d array I can do this by taking a
    #2d array of the halo centers and doing center_rock_select[:,np.newaxis]
    #which makes a 3d array where each array of axis=0 is the halo center
    
    hi_res_ids = []
    hi_res_mass = []
    hi_res_vmax = []
    hi_res_rvir = []
    hi_res_centers = []
    hi_res_fraction_list = []
    hi_res_M_star_rvir = []
    hi_res_M_star_gal = []
    hi_res_M_gas_rvir = []
    hi_res_M_gas_rvir_cold = []
    closest_low_res_particle = []

    #First do all the calculations associated with the low_res_particles

    #okay this just takes alot of memory, so I'm just going to use a loop
    print 'looping over halos'
    for halo_id in range(len(center_rock_select)):
        
        halo_particle_diff = center_rock_select[halo_id] - low_res_coords
        low_res_dist_all = np.linalg.norm(halo_particle_diff,axis=1)
        assert len(low_res_dist_all) == len(low_res_coords)
        low_res_n_in_rvir = np.sum((low_res_dist_all<rvir_rock_select[halo_id]))
        low_res_mass_in_rvir = low_res_n_in_rvir*low_res_mass

        hi_res_fraction = mass_rock_select[halo_id]/(mass_rock_select[halo_id]+low_res_mass_in_rvir)

        if hi_res_fraction>=Low_res_tolerance:

            #print 'calculate hi-res particles '+str(halo_id)+' :'

            star_diff = center_rock_select[halo_id] - star_coords
            star_dist_all = np.linalg.norm(star_diff,axis=1)

            star_masked_distances = (star_dist_all<rvir_rock_select[halo_id])
            star_masked_distances_gal = (star_dist_all<rgal_rock_select[halo_id])

            star_mass_in_rvir = np.sum(star_mass[star_masked_distances])
            star_mass_in_rgal = np.sum(star_mass[star_masked_distances_gal])

            #gas calculations
            gas_diff = center_rock_select[halo_id] - gas_coords
            gas_dist_all = np.linalg.norm(gas_diff,axis=1)

            gas_masked_distances = (gas_dist_all<rvir_rock_select[halo_id])
            gas_masked_distances_gal = (gas_dist_all<rgal_rock_select[halo_id])

            cold_gas_masked_distances = (gas_dist_all<rvir_rock_select[halo_id])
            cold_gas_masked_distances_gal = (gas_dist_all<rgal_rock_select[halo_id])

            gas_mass_in_rvir = np.sum(gas_mass[gas_masked_distances])
            gas_mass_in_gal = np.sum(gas_mass[gas_masked_distances_gal])

            cold_gas_in_rvir = np.sum(gas_mass[cold_gas_masked_distances])
            cold_gas_in_rgal = np.sum(gas_mass[cold_gas_masked_distances_gal])

            hi_res_ids.append(id_rock_select[halo_id])
            hi_res_mass.append(mass_rock_select[halo_id])
            hi_res_vmax.append(vmax_rock_select[halo_id])
            hi_res_rvir.append(rvir_rock_select[halo_id])
            hi_res_centers.append(center_rock_select[halo_id])
            hi_res_fraction_list.append(hi_res_fraction)
            hi_res_M_star_rvir.append(star_mass_in_rvir)
            hi_res_M_star_gal.append(star_mass_in_rgal)
            hi_res_M_gas_rvir.append(gas_mass_in_rvir)
            hi_res_M_gas_rvir_cold.append(cold_gas_in_rvir)
            closest_low_res_particle.append(np.min(low_res_dist_all))
            
    #now I want to make a dictionary that has each halo and its properties in it

    print 'total halos: '+str(len(id_rock_select))
    print 'hi res halos: '+str(len(hi_res_ids))

    halo_dict = {'halo_ids':hi_res_ids,'M_vir':hi_res_mass,'V_max':hi_res_vmax,'R_vir':hi_res_rvir,'centers':hi_res_centers,'hi_res_fraction':hi_res_fraction_list,'M_star_gal':hi_res_M_star_gal,'M_star_vir':hi_res_M_star_rvir,'M_gas_rvir':hi_res_M_gas_rvir,'M_gas_cold_rvir':hi_res_M_gas_rvir_cold,'Closest_low_res_particle':closest_low_res_particle}

    return halo_dict

def Load_hi_res_halo_particle_dict(giz_hdf5,halo_file,Low_res_tolerance=1.0,save_particles=False,halo_lim=1.0e9,h=0.702,hf_type='rockstar'):
    #This is the same program as above, but now it saves the particle data
    #for each halo above your halo limit
    import numpy as np

    if hf_type=='rockstar':
        halo_dict = Load_Halo_Data(halo_file,h=h)
    elif hf_type=='AHF':
        halo_dict = Load_Halo_Data_AHF(halo_file,h=h)
        
    PD_dict = Load_Particle_Data(giz_hdf5,add_low_res=True,add_stars=True,add_gas=True)

    id_rock = halo_dict['ids'] 
    M_rock = halo_dict['masses']
    Vmax_rock = halo_dict['vmax']
    Rvir_rock = halo_dict['rvir']
    Rmax_rock = halo_dict['rmax']
    centers_rock = halo_dict['centers']
    velocities_rock = halo_dict['velocities']
    if hf_type=='AHF':
        f_hires_rock = halo_dict['f_hires']
    else:
        f_hires_rock = np.ones_like(id_rock)

    mass_select = (M_rock>halo_lim)&(f_hires_rock>=Low_res_tolerance)
    
    id_rock_select = id_rock[mass_select]
    mass_rock_select = M_rock[mass_select]
    vmax_rock_select = Vmax_rock[mass_select]
    rvir_rock_select = Rvir_rock[mass_select]
    rmax_rock_select = Rmax_rock[mass_select]
    center_rock_select = centers_rock[mass_select]
    velocities_rock_select = velocities_rock[mass_select]
    f_hires_select = f_hires_rock[mass_select]

    print 'total halos: {}, selected halos: {}'.format(str(len(id_rock)),str(len(id_rock_select)))

    rgal_rock_select = rvir_rock_select*0.1 #galaxy radius defined as just 10% Rvir
    
    #hi_res_coords = PD_dict['halo']['coords']
    #hi_res_mass = PD_dict['halo']['masses'][0]

    low_res_coords = PD_dict['disk']['coords']
    low_res_mass = PD_dict['disk']['masses'][0]

    star_coords = PD_dict['star']['coords']
    star_mass = PD_dict['star']['masses']
    star_ages = PD_dict['star']['ages']    

    gas_coords = PD_dict['gas']['coords']
    gas_mass = PD_dict['gas']['masses']
    gas_temp = PD_dict['gas']['temperature']

    #Now I want to find the distance between every halo and every particle
    #I can do this "easily" by using masked arrays
    #
    #The first step is to create a 3d array I can do this by taking a
    #2d array of the halo centers and doing center_rock_select[:,np.newaxis]
    #which makes a 3d array where each array of axis=0 is the halo center
    
    hi_res_ids = []
    hi_res_mass = []
    hi_res_vmax = []
    hi_res_rvir = []
    hi_res_centers = []
    hi_res_fraction_list = []
    hi_res_M_star_rvir = []
    hi_res_M_star_gal = []
    hi_res_M_gas_rvir = []
    hi_res_M_gas_rvir_cold = []
    closest_low_res_particle = []

    #First do all the calculations associated with the low_res_particles

    #okay this just takes alot of memory, so I'm just going to use a loop
    print 'looping over halos'

    hi_res_particle_dict={}

    for halo_id in range(len(center_rock_select)):
        
        halo_particle_diff = center_rock_select[halo_id] - low_res_coords
        low_res_dist_all = np.linalg.norm(halo_particle_diff,axis=1)
        assert len(low_res_dist_all) == len(low_res_coords)
        low_res_n_in_rvir = np.sum((low_res_dist_all<rvir_rock_select[halo_id]))
        low_res_mass_in_rvir = low_res_n_in_rvir*low_res_mass

        hi_res_fraction = mass_rock_select[halo_id]/(mass_rock_select[halo_id]+low_res_mass_in_rvir)


        if hi_res_fraction>=Low_res_tolerance:
            if hi_res_fraction<f_hires_select[halo_id]:
                print "Warning: hali {} has an AHF fraction of {}, but a calculated fraction of {}".format(str(halo_id),str(f_hires_select[halo_id]),str(hi_res_fraction))

            #print 'calculate hi-res particles '+str(halo_id)+' :'

            star_diff = center_rock_select[halo_id] - star_coords
            star_dist_all = np.linalg.norm(star_diff,axis=1)

            star_masked_distances = (star_dist_all<rvir_rock_select[halo_id])
            star_masked_distances_gal = (star_dist_all<rgal_rock_select[halo_id])

            star_mass_in_rvir = np.sum(star_mass[star_masked_distances])
            star_mass_in_rgal = np.sum(star_mass[star_masked_distances_gal])

            #gas calculations
            gas_diff = center_rock_select[halo_id] - gas_coords
            gas_dist_all = np.linalg.norm(gas_diff,axis=1)

            gas_masked_distances = (gas_dist_all<rvir_rock_select[halo_id])
            gas_masked_distances_gal = (gas_dist_all<rgal_rock_select[halo_id])

            cold_gas_masked_distances = (gas_dist_all<rvir_rock_select[halo_id])
            cold_gas_masked_distances_gal = (gas_dist_all<rgal_rock_select[halo_id])

            gas_mass_in_rvir = np.sum(gas_mass[gas_masked_distances])
            gas_mass_in_gal = np.sum(gas_mass[gas_masked_distances_gal])

            cold_gas_in_rvir = np.sum(gas_mass[cold_gas_masked_distances])
            cold_gas_in_rgal = np.sum(gas_mass[cold_gas_masked_distances_gal])

            hi_res_ids.append(id_rock_select[halo_id])
            hi_res_mass.append(mass_rock_select[halo_id])
            hi_res_vmax.append(vmax_rock_select[halo_id])
            hi_res_rvir.append(rvir_rock_select[halo_id])
            hi_res_centers.append(center_rock_select[halo_id])
            hi_res_fraction_list.append(hi_res_fraction)
            hi_res_M_star_rvir.append(star_mass_in_rvir)
            hi_res_M_star_gal.append(star_mass_in_rgal)
            hi_res_M_gas_rvir.append(gas_mass_in_rvir)
            hi_res_M_gas_rvir_cold.append(cold_gas_in_rvir)
            closest_low_res_particle.append(np.min(low_res_dist_all))

            #now lets actually save some halo properties
            hi_res_particle_dict.update({str(halo_id):{'coords':star_diff[star_masked_distances],'ages':star_ages[star_masked_distances]}})
            
    #now I want to make a dictionary that has each halo and its properties in it

    print 'total halos: '+str(len(id_rock_select))
    print 'hi res halos: '+str(len(hi_res_ids))

    halo_dict = {'halo_ids':hi_res_ids,'M_vir':hi_res_mass,'V_max':hi_res_vmax,'R_vir':hi_res_rvir,'centers':hi_res_centers,'hi_res_fraction':hi_res_fraction_list,'M_star_gal':hi_res_M_star_gal,'M_star_vir':hi_res_M_star_rvir,'M_gas_rvir':hi_res_M_gas_rvir,'M_gas_cold_rvir':hi_res_M_gas_rvir_cold,'Closest_low_res_particle':closest_low_res_particle}

    return halo_dict, hi_res_particle_dict

def return_specific_halo(halo_file,halo_id):
    #given a rockstar halo file and an id 
    #it will return that halo's center, virial radius, and velocity
    import numpy as np
    import yt, h5py, re, os
    from astropy.cosmology import FlatLambdaCDM

    #This script identifies all the halos above 1.0e10 (that's hard coded) and calculates all of the things
    #that are hi-res and then returns the center of the highest mass halo that doesn't have a low res particle
    #within Rvir

    halo_dict = Load_Halo_Data(halo_file) #load in the rockstar filw

    id_rock = halo_dict['ids'] 
    M_rock = halo_dict['masses']
    Vmax_rock = halo_dict['vmax']
    Rvir_rock = halo_dict['rvir']
    Rmax_rock = halo_dict['rmax']
    centers_rock = halo_dict['centers']
    velocities_rock = halo_dict['velocities']

    id_select = (id_rock==halo_id) #This is hard coded in

    id_rock_select = id_rock[id_select]
    mass_rock_select = M_rock[id_select]
    vmax_rock_select = Vmax_rock[id_select]
    rvir_rock_select = Rvir_rock[id_select]
    rmax_rock_select = Rmax_rock[id_select]
    center_rock_select = centers_rock[id_select]
    velocities_rock_select = velocities_rock[id_select]

    print 'id_gal: '+str(id_rock_select[0])
    print 'mass: '+'{0:.2e} Msun'.format(mass_rock_select[0])
    print 'radius: '+'{0:.2e} kpc'.format(rvir_rock_select[0])

    return center_rock_select[0], rvir_rock_select[0], velocities_rock_select[0]

def galaxy_statistics(giz_hdf5,halo_file,print_values=False,halo_id=None):
    import numpy as np
    import yt, h5py, re, os
    from math import log10
    from astropy.cosmology import FlatLambdaCDM
    from andrew_hydro_sim_modules.simple_tools import get_distance
    
    #print some simple stats of the galaxy, and then return the star particles within Rvir
    #specifically their ages (to compute SFHs) and FeH.
    if halo_id == None:
        host_center, host_rvir = Identify_Host(giz_hdf5,halo_file)
    else:
        host_center, host_rvir, host_vel = return_specific_halo(halo_file,int(halo_id))
    
    print 'loading pd'

    PD_dict = Load_Particle_Data(giz_hdf5,add_dm=False,add_gas=True,add_stars=True)

    star_coords = PD_dict['star']['coords']
    star_masses = PD_dict['star']['masses']
    star_age = PD_dict['star']['ages']
    star_Fe_H = PD_dict['star']['Fe_H']

    gas_coords = PD_dict['gas']['coords']
    gas_masses = PD_dict['gas']['masses']
    gas_temps = PD_dict['gas']['temperature']

    star_dists = get_distance(star_coords,host_center)
    gas_dists = get_distance(gas_coords,host_center)

    star_dist_mask = (star_dists<=0.1*host_rvir)
    gas_dist_mask = (gas_dists<=0.1*host_rvir)

    star_dist_mask_rvir = (star_dists<=host_rvir)
    gas_dist_mask_rvir = (gas_dists<=host_rvir)

    cold_gas_mask = (gas_dists<=host_rvir)&(gas_temps<=1.0e4)

    galaxy_star_parts_mass = star_masses[star_dist_mask]
    galaxy_star_part_ages = star_age[star_dist_mask]
    galaxy_star_part_FeH = star_Fe_H[star_dist_mask]
    galaxy_gas_mass = gas_masses[gas_dist_mask]
    cold_gas_mass = gas_masses[cold_gas_mask]

    gas_mass_rvir = gas_masses[gas_dist_mask_rvir]
    star_mass_rvir = star_masses[star_dist_mask_rvir]
    
    if print_values==True:
        print '{:#^30} \n'.format('This returns a list of ages and FeH for stars in Rvir')

        print 'Galaxy stellar mass (in 0.1*Rvir) is: {0:.2e}'.format(sum(galaxy_star_parts_mass))
        print 'Galaxy gas mass (in 0.1*Rvir) is: {0:.2e}'.format(sum(galaxy_gas_mass))

        print 'Galaxy stellar mass (in Rvir) is: {0:.2e}'.format(sum(star_mass_rvir))
        print 'Galaxy gas mass (in Rvir) is: {0:.2e}'.format(sum(gas_mass_rvir))

        print 'Cold gas mass (in Rvir) is: {0:.2e}'.format(sum(cold_gas_mass))

    agebins = np.linspace(0.0,14.0,100)
    cosmo = FlatLambdaCDM(H0=70.2,Om0=0.266) #THIS SHOULDN'T BE HARD CODED
    
    #theres this annoying thing where the age is a float so .value doesn't work in older versions of astropy

    star_age_gal_T = [cosmo.age(1.0/xx - 1.0).value for xx in galaxy_star_part_ages] #converts scale factor to time given cosmology

    #star_age_gal_T = [cosmo.age(1.0/xx - 1.0).value for xx in galaxy_star_part_ages] #converts scale factor to time given cosmology

    return star_age_gal_T, galaxy_star_part_FeH, galaxy_star_parts_mass, galaxy_gas_mass

def principle_axes(coordinates,masses,center,rad):
    #This code calculates principle axes of a given star particle
    #This is basically a modified version of Andrew Wetzel's code
    #to do the same thing
    import numpy as np
    import yt, h5py, re, os
    from math import log10
    from astropy.cosmology import FlatLambdaCDM
    from andrew_hydro_sim_modules.simple_tools import get_distance_vector, get_distance

    dm_dist_val = get_distance(coordinates, center)

    dist_mask = (dm_dist_val<=rad)
    coord_mod = coordinates[dist_mask]
    mass_mod = masses[dist_mask]

    dm_dist = get_distance_vector(coord_mod, center)

    weights = mass_mod/np.median(mass_mod)

    xx = np.sum(weights * dm_dist[:,0]**2.0)
    yy = np.sum(weights * dm_dist[:,1]**2.0)
    zz = np.sum(weights * dm_dist[:,2]**2.0)
    
    xy = yx = np.sum(weights * dm_dist[:,0] * dm_dist[:,1])
    xz = zx = np.sum(weights * dm_dist[:,0] * dm_dist[:,2])
    yz = zy = np.sum(weights * dm_dist[:,1] * dm_dist[:,2])
    
    I_tensor = [[xx, xy, xz],[yx, yy, yz],[zx, zy, zz]]

    eigen_values, eigen_vectors = np.linalg.eig(I_tensor)

    # order eigen-vectors by eigen-values, from largest to smallest                                  
    eigen_indices_sorted = np.argsort(eigen_values)[::-1]
    eigen_values = eigen_values[eigen_indices_sorted]
    eigen_values /= eigen_values.max()  # renormalize to 1                                           
    # make eigen_vectors[0] corresponds to vector of eigen_values[0]                                 
    eigen_vectors = eigen_vectors.transpose()[eigen_indices_sorted]

    axis_ratios = np.sqrt(
        [eigen_values[2] / eigen_values[0],
         eigen_values[2] / eigen_values[1],
         eigen_values[1] / eigen_values[0]]
    )

    return eigen_values, eigen_vectors, axis_ratios

def report_eigen(giz_hdf5,halo_file,add_dm=True,add_gas=False,add_stars=False,add_low_res=False):
    
    print '\n{:#^30}'.format('CALCULATING MOMENT OF INTERTIA TENSOR')

    host_center, host_rvir = Identify_Host(giz_hdf5,halo_file)
    PD_dict = Load_Particle_Data(giz_hdf5,add_dm=add_dm,add_gas=add_gas,add_stars=add_stars)
    if add_dm==True:
        dm_part = PD_dict['halo']['coords']
        dm_mass = PD_dict['halo']['masses']
        
        dm_eigen_val, dm_eigen_vec, dm_axis_ratios = principle_axes(dm_part,dm_mass,host_center,host_rvir)
    
        print '\nThe eigen values of the dark matter halo: {0}, {1}, {2}'.format(dm_eigen_val)
        print 'Sorted eigen vectors: \neignevector one: {0} \neigenvector two: {1} \neigenvector three: {2}'.format(dm_eigen_vec)
        print 'axis ratios: {0}, {1}, {2}'.format(dm_axis_ratios)
    
    if add_gas==True:
        gas_part = PD_dict['gas']['coords']
        gas_mass = PD_dict['gas']['masses']

        gas_eigen_val, gas_eigen_vec, gas_axis_ratios = principle_axes(gas_part,gas_mass,host_center,host_rvir)

        print '\nThe eigen values of the gas: {0}, {1}, {2}'.format(gas_eigen_val)
        print 'Sorted eigen vectors: \neignevector one: {0} \neigenvector two: {1} \neigenvector three: {2}'.format(gas_eigen_vec)
        print 'axis ratios: {0}, {1}, {2}'.format(gas_axis_ratios)

    if add_stars==True:
        star_part = PD_dict['star']['coords']
        star_mass = PD_dict['star']['masses']

        star_eigen_val, star_eigen_vec, star_axis_ratios = principle_axes(star_part,star_mass,host_center,host_rvir)

        print '\nThe eigen values of the stars: {a[0]}, {a[1]}, {a[2]}'.format(a=star_eigen_val)
        print 'Sorted eigen vectors: \neignevector one: {a[0]} \neigenvector two: {a[1]} \neigenvector three: {a[2]}'.format(a=star_eigen_vec)
        print 'axis ratios: {a[0]}, {a[1]}, {a[2]}'.format(a=star_axis_ratios)
    

def Calc_average_L(coordinates,masses,velocities,rad,center=[0.0,0.0,0.0],host_vel=[0.0,0.0,0.0]):
    import numpy as np
    import yt, h5py, re, os
    from math import log10
    from astropy.cosmology import FlatLambdaCDM
    from andrew_hydro_sim_modules.simple_tools import get_distance_vector, get_distance

    dm_dist_val = get_distance(coordinates, center)

    dist_mask = (dm_dist_val<=rad)
    coord_mod = coordinates[dist_mask]
    mass_mod = masses[dist_mask]
    vel_mod = velocities[dist_mask]

    coord_shift = coord_mod-center
    vel_shift = vel_mod-host_vel

    rcrossv = np.cross(coord_shift,vel_shift)
    L_vec = [mass_mod[ii]*rcrossv[ii] for ii in range(len(rcrossv))]

    L_avg = np.mean(L_vec,axis=0)
    den = np.sqrt(L_avg[0]**2.0+L_avg[1]**2.0+L_avg[2]**2.0)
    return L_avg/den

def Calc_average_L_shift(coordinates,masses,velocities):
    #For already shifted and cut cordinates
    import numpy as np
    import yt, h5py, re, os
    from math import log10
    from astropy.cosmology import FlatLambdaCDM
    from andrew_hydro_sim_modules.simple_tools import get_distance_vector, get_distance

    coord_mod = coordinates
    mass_mod = masses
    vel_mod = velocities

    coord_shift = coord_mod
    vel_shift = vel_mod

    rcrossv = np.cross(coord_shift,vel_shift)
    L_vec = [mass_mod[ii]*rcrossv[ii] for ii in range(len(rcrossv))]

    L_avg = np.mean(L_vec,axis=0)
    den = np.sqrt(L_avg[0]**2.0+L_avg[1]**2.0+L_avg[2]**2.0)
    return L_avg/den

def report_average_L(giz_hdf5,halo_file,add_dm=True,add_gas=False,add_stars=False,add_low_res=False):
    import numpy as np
    import yt, h5py, re, os
    from astropy.cosmology import FlatLambdaCDM

    host_center, host_rvir, host_vel = Identify_Host(giz_hdf5,halo_file,add_velocity=True)
    PD_dict = Load_Particle_Data(giz_hdf5,add_dm=add_dm,add_gas=add_gas,add_stars=add_stars)
    
    L_vec = np.zeros((5,3))

    if add_dm == True:
        dm_part = PD_dict['halo']['coords']
        dm_mass = PD_dict['halo']['masses']
        dm_vel = PD_dict['halo']['velocities']

        L_dm = Calc_average_L(dm_part,dm_mass,dm_vel,host_rvir,host_center,host_vel)
        L_vec[1] = L_dm
        print 'dark matter <L> = {a[0]}, {a[1]}, {a[2]}'.format(a=L_vec[1])

    if add_gas == True:
        gas_part = PD_dict['gas']['coords']
        gas_mass = PD_dict['gas']['masses']
        gas_vel = PD_dict['gas']['velocities']

        L_gas = Calc_average_L(gas_part,gas_mass,gas_vel,host_rvir,host_center,host_vel)
        L_vec[0] = L_gas
        print 'gas <L> = {a[0]}, {a[1]}, {a[2]}'.format(a=L_vec[0])

    if add_stars == True:
        stars_part = PD_dict['star']['coords']
        stars_mass = PD_dict['star']['masses']
        stars_vel = PD_dict['star']['velocities']

        L_star = Calc_average_L(stars_part,stars_mass,stars_vel,host_rvir,host_center,host_vel)
        L_vec[3] = L_star
        print 'stars <L> = {a[0]}, {a[1]}, {a[2]}'.format(a=L_vec[3])

    return L_vec

def convert_to_cylindrical(coordinate,velocity,halo_center=[0.0,0.0,0.0],halo_vel=[0.0,0.0,0.0]):
    #convert a 3 element matrix into a matrix of spherical coordinates
    #I could calculate V_tan and then project it into the new xy plane

    #it also could just be cylindrical and just use V_theta
    #just do cylindrical

    #the basic formula for converting velocities from cartesian to cylindrical is
    #
    # v = (x*dx/dt+y*dy/dt)/r rhat + (x*dy/dt-y*dx/dt)/r thetahat + dz/dt zhat

    import numpy as np

    halo_center = np.asarray(halo_center)
    halo_vel = np.asarray(halo_vel)

    X = np.asarray(coordinate[:,0])-halo_center[0]
    Y = np.asarray(coordinate[:,1])-halo_center[1]
    Z = np.asarray(coordinate[:,2])-halo_center[2]

    VX = np.asarray(velocity[:,0])-halo_vel[0]
    VY = np.asarray(velocity[:,1])-halo_vel[1]
    VZ = np.asarray(velocity[:,2])-halo_vel[2]

    R = np.sqrt(X**2.0+Y**2.0)
    theta = np.arctan2(Y,X)

    D_cylindrical = np.zeros((len(R),3))
    D_cylindrical[:,0] = R
    D_cylindrical[:,1] = theta
    D_cylindrical[:,2] = Z

    Vr = (X*VX+Y*VY)/R
    Vtheta = (X*VY-Y*VX)/R
    Vz = VZ

    V_cylindrical = np.zeros((len(Vr),3))
    V_cylindrical[:,0] = Vr
    V_cylindrical[:,1] = Vtheta
    V_cylindrical[:,2] = Vz

    return D_cylindrical, V_cylindrical

def Rotate_to_z_axis(coordinates,velocities,rotation_axis):
    import numpy as np
    import yt, h5py, re, os
    from math import log10
    from astropy.cosmology import FlatLambdaCDM
    from andrew_hydro_sim_modules.simple_tools import get_distance_vector, get_distance
    #Okay I want to take in a "z" axis, and then rotate the
    #coordinates so that that is the z axis
    #then calculate velocity vectors in that frame and
    #then decompose it into spherical coordinates IN THAT FRAME

    L = np.sqrt(rotation_axis[0]**2.0+rotation_axis[1]**2.0+rotation_axis[2]**2.0) #total length
    R = np.sqrt(rotation_axis[0]**2.0+rotation_axis[1]**2.0) #length in xy plane
    
    R1 = np.asarray([[rotation_axis[0]/R,rotation_axis[1]/R,0.0],[-rotation_axis[1]/R,rotation_axis[0],0.0],[0.0,0.0,1.0]]) #rotation about z axis to project into xz plane
    R2 = np.asarray([[rotation_axis[2]/L,0.0,-R/L],[0.0,1.0,0.0],[R/L,0.0,rotation_axis[2]/L]]) #rotation about y axis to make given axis the z axis

    #apply rotation to coordinates and velocities
    
    coord_rotate = np.asarray([R2.dot(R1.dot(xx)) for xx in coordinates])
    vel_rotate = np.asarray([R2.dot(R1.dot(xx)) for xx in velocities])

    return coord_rotate, vel_rotate
    
def report_velocities(giz_hdf5,halo_file,add_dm=True,add_gas=False,add_stars=False,add_low_res=False,vector=None,report_vel=True,halo_id=None):
    import numpy as np
    import yt, h5py, re, os
    from math import log10
    from astropy.cosmology import FlatLambdaCDM
    from andrew_hydro_sim_modules.simple_tools import get_distance_vector, get_distance
    #1) calculate host halo
    #2) calculate average L
    #3) rotate to L
    #4) convert to cylindrical

    #1) calculate host halo

    if halo_id == None:
        host_center, host_rvir, host_vel = Identify_Host(giz_hdf5,halo_file,add_velocity=True)
    else:
        host_center, host_rvir, host_vel = return_specific_halo(halo_file,int(halo_id))

    #2) calculate average L
    PD_dict = Load_Particle_Data(giz_hdf5,add_dm=add_dm,add_gas=add_gas,add_stars=add_stars)
    
    if add_dm == True:
        dm_part = PD_dict['halo']['coords']
        dm_mass = PD_dict['halo']['masses']
        dm_vel = PD_dict['halo']['velocities']

        dm_dist_val = get_distance(dm_part, host_center)
        dist_cut = (dm_dist_val<=host_rvir)
        coord_cut = dm_part[dist_cut]
        vel_cut = dm_vel[dist_cut]
        mass_cut = dm_mass[dist_cut]

        coord_shift = coord_cut-host_center
        vel_shift = vel_cut-host_vel
        
        if vector==None:
            L_dm = Calc_average_L_shift(coord_shift,mass_cut,vel_shift)
        else: 
            L_dm = vector

        print 'the vector shape is {}'.format(L_dm.shape)
        assert L_dm.shape == (1,3)
        #Rotate to L

        part_rotate, vel_rotate =  Rotate_to_z_axis(coord_shift,vel_shift,L_dm)

        #convert to cylindrical these are already cut and shifted
        
        part_cyl, vel_cyl = convert_to_cylindrical(part_rotate,vel_rotate,halo_center=[0.0,0.0,0.0],halo_vel=[0.0,0.0,0.0])
        
        V_theta = vel_cyl[:,1]

        print 'The rotational velocity is {} km/s'.format(np.mean(V_theta))

    if add_gas == True:
        gas_part = PD_dict['gas']['coords']
        gas_mass = PD_dict['gas']['masses']
        gas_vel = PD_dict['gas']['velocities']

        gas_dist_val = get_distance(gas_part, host_center)
        dist_cut = (gas_dist_val<=host_rvir)
        coord_cut = gas_part[dist_cut]
        vel_cut = gas_vel[dist_cut]
        mass_cut = gas_mass[dist_cut]

        coord_shift = coord_cut-host_center
        vel_shift = vel_cut-host_vel
        
        if vector==None:
            L_gas = Calc_average_L_shift(coord_shift,mass_cut,vel_shift)
        else: 
            L_gas = vector

        print 'the vector shape is {}'.format(L_gas.shape)
        assert L_gas.shape == (1,3)

        #Rotate to L

        part_rotate, vel_rotate =  Rotate_to_z_axis(coord_shift,vel_shift,L_gas)

        #convert to cylindrical these are already cut and shifted
        
        part_cyl, vel_cyl = convert_to_cylindrical(part_rotate,vel_rotate,halo_center=[0.0,0.0,0.0],halo_vel=[0.0,0.0,0.0])
        
        V_theta = vel_cyl[:,1]

        print 'The rotational velocity is {} km/s'.format(np.mean(V_theta))
    
    if add_stars == True:
        stars_part = PD_dict['star']['coords']
        stars_mass = PD_dict['star']['masses']
        stars_vel = PD_dict['star']['velocities']

        stars_dist_val = get_distance(stars_part, host_center)
        dist_cut = (stars_dist_val<=host_rvir)
        coord_cut = stars_part[dist_cut]
        vel_cut = stars_vel[dist_cut]
        mass_cut = stars_mass[dist_cut]

        vel_disp_one = np.var(np.sqrt(vel_cut[:,0]**2.0+vel_cut[:,1]**2.0+vel_cut[:,2]**2.0))

        coord_shift = coord_cut-host_center
        vel_shift = vel_cut-host_vel

        vel_disp_two = np.var(np.sqrt(vel_shift[:,0]**2.0+vel_shift[:,1]**2.0+vel_shift[:,2]**2.0))

        L_stars = Calc_average_L_shift(coord_shift,mass_cut,vel_shift)

        #Rotate to L

        if vector==None:
            L_stars = Calc_average_L_shift(coord_shift,mass_cut,vel_shift)
        else: 
            L_stars = vector

        assert L_stars.shape==tuple((3,))

        part_rotate, vel_rotate =  Rotate_to_z_axis(coord_shift,vel_shift,L_stars)

        vel_disp_three = np.var(np.sqrt(vel_rotate[:,0]**2.0+vel_rotate[:,1]**2.0+vel_rotate[:,2]**2.0))        

        #convert to cylindrical these are already cut and shifted
        
        part_cyl, vel_cyl = convert_to_cylindrical(part_rotate,vel_rotate,halo_center=[0.0,0.0,0.0],halo_vel=[0.0,0.0,0.0])
        
        vel_disp_four = np.var(np.sqrt(vel_cyl[:,0]**2.0+vel_cyl[:,1]**2.0+vel_cyl[:,2]**2.0))

        V_theta = vel_cyl[:,1]

        if report_vel==True:
            print 'The rotational velocity is {} km/s'.format(np.mean(V_theta))
            print 'Velocity dispersion is {} km/s'.format(np.sqrt(vel_disp_four))
            print 'V/sigma = {}'.format(np.mean(V_theta)/np.sqrt(vel_disp_four))

        return np.mean(V_theta), np.sqrt(vel_disp_four), np.mean(V_theta)/np.sqrt(vel_disp_four)
