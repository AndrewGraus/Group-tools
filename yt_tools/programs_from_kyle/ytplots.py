from ytstewart import SFR2, spin_parameters
#import yt.utilities.cosmology as cosmo
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.279)
import numpy as np
import matplotlib.pyplot as plt
from astropy.convolution import convolve, Gaussian1DKernel
import pylab as P

def SFR_vs_time():
    E_sfr, E_time = np.loadtxt('Enzo_SFRvstime.txt', skiprows=1, unpack=True)
    R_sfr, R_time = np.loadtxt('Ramses_SFRvstime.txt', skiprows=1, unpack=True)
    A_sfr, A_time = np.loadtxt('Arepo_SFRvstime.txt', skiprows=1, unpack=True)
    Art_sfr, Art_time = np.loadtxt('Art_SFRvstime.txt', skiprows=1, unpack=True)
    F_sfr, F_time = np.loadtxt('Fire_SFRvstime.txt', skiprows=1, unpack=True)

    z_tick_values = [6, 5, 4, 3, 2, 1.5, 1.25, 1.0 ]  #max and min of this are limits of the plot.  
    z_tick_value_labels = [' ', '5', '4', '3', '2', '1.5', '1.25', '1.0' ]  
    z_tick_locations = cosmo.lookback_time(z_tick_values).value
    LBT_max = cosmo.lookback_time(np.max(z_tick_values)).value
    LBT_min = cosmo.lookback_time(np.min(z_tick_values)).value

    fig = plt.figure()
    plt.ylabel('SFR  [M$_\odot$ yr$^{-1}$]')
    plt.ylim(0, 90)
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twiny()
    ax1.set_xlabel('Lookback Time [Gyr]')
    ax1.set_xlim(LBT_max, LBT_min)
    ax2.set_xlabel('Redshift')
    ax2.set_xlim(LBT_max, LBT_min)
    ax2.set_xticks(z_tick_locations)
    ax2.set_xticklabels(z_tick_value_labels)
    
    E_LBtime = cosmo.age(0).value-E_time/1e9
    ax1.plot(E_LBtime,E_sfr, linewidth=2)
    R_LBtime = cosmo.age(0).value-R_time/1e9
    ax1.plot(R_LBtime,R_sfr, ':', linewidth=3)
    Art_LBtime = cosmo.age(0).value-Art_time/1e9
    ax1.plot( Art_LBtime,convolve(Art_sfr,Gaussian1DKernel(1)), '-', linewidth=2 )
    A_LBtime = cosmo.age(0).value-A_time/1e9
    ax1.plot(A_LBtime,convolve(A_sfr,Gaussian1DKernel(0.5)), '-.', linewidth=2 )
    F_LBtime = cosmo.age(0).value-F_time/1e9
    ax1.plot(F_LBtime,F_sfr, linewidth=2)

    #str_z=str(abs(round(E_pf.current_redshift, 1)))
    ax1.legend(("Enzo","Ramses", "Art", "Arepo", "Gizmo-pSPH"), frameon=False)
    plt.savefig("ALL_SFR_vs_time.png", bbox_inches='tight')

#reads in sipn parameters from saved files
#order of saved files is: tot, dark, stars, particles, gas, cold, hot
def spinparam_2x2():
    #Ez=[0.00, 0.05, 0.10, 0.20, 0.25, 0.30, 0.40, 0.50, 0.75, 1.00, 3.00]
    #Rz=[0.00, 0.25, 0.31, 0.39, 0.50, 0.75, 1.00, 1.50, 2.01, 2.64, 3.00, 3.55, 4.00, 5.06]
    #Rz=[1.00, 1.50, 2.01, 2.64, 3.00, 3.55, 4.00, 5.06]
    #Az=[1.532, 1.82, 3.00, 4.94]
    #Artz=[2.03, 2.99, 4.0]
    #Fz=[0.80, 0.90, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.20, 2.50, 2.60, 2.80, 3.00, 3.40, 3.82, 3.90, 4.00, 4.50, 5.00]
    #Fz=[1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.20, 2.50, 2.60, 2.80, 3.00, 3.40, 3.82, 3.90, 4.00, 4.50, 5.00]
    #Fz=[1.70, 1.80, 1.90, 2.00, 2.20, 2.50, 2.60, 2.80, 3.00, 3.40, 3.82, 3.90, 4.00, 4.50, 5.00]
    Ez=[1.00, 3.00]
    Rz=[1.00, 1.50, 2.01, 2.64, 3.00]
    Az=[1.532, 1.82, 2.0, 3.00]
    Artz=[2.03, 2.33, 2.99]
    Fz=[1.00, 1.20, 1.40, 1.60, 1.80, 2.00, 2.20, 2.50, 2.60, 2.80, 3.00]
    E_LBtime=cosmo.lookback_time(Ez)
    R_LBtime=cosmo.lookback_time(Rz)
    A_LBtime=cosmo.lookback_time(Az)
    Art_LBtime=cosmo.lookback_time(Artz)
    F_LBtime=cosmo.lookback_time(Fz)

    endstring = "_spinparameters.txt" #this will average over everything in the halo
    
    #create 2-D numpy arrays.  First index is redshift, 2nd index is which lambda parameter:  Xz_array[redshift,lambda]
    simulation="Enzo"
    redshift_array=Ez
    Ez_array = np.zeros(E_LBtime.shape[0]*7).reshape(E_LBtime.shape[0], 7)
    for i in range(E_LBtime.shape[0]):
        fname="saved_data/"+simulation+"_z"+str(abs(round(redshift_array[i], 2)))+endstring
        names, trash, Ez_array[i][:] = np.loadtxt(fname, dtype={'names':('variable_name','equalsign','value'), 'formats':('|S15','|S8', np.float)}, unpack=True)

    simulation="Ramses"
    redshift_array=Rz
    Rz_array = np.zeros(R_LBtime.shape[0]*7).reshape(R_LBtime.shape[0], 7)
    for i in range(R_LBtime.shape[0]):
        fname="saved_data/"+simulation+"_z"+str(abs(round(redshift_array[i], 2)))+endstring
        names, trash, Rz_array[i][:] = np.loadtxt(fname, dtype={'names':('variable_name','equalsign','value'), 'formats':('|S15','|S8', np.float)}, unpack=True)

    simulation="Art"
    redshift_array=Artz
    Artz_array = np.zeros(Art_LBtime.shape[0]*7).reshape(Art_LBtime.shape[0], 7)
    for i in range(Art_LBtime.shape[0]):
        fname="saved_data/"+simulation+"_z"+str(abs(round(redshift_array[i], 2)))+endstring
        names, trash, Artz_array[i][:] = np.loadtxt(fname, dtype={'names':('variable_name','equalsign','value'), 'formats':('|S15','|S8', np.float)}, unpack=True)

    simulation="Arepo"
    redshift_array=Az
    Az_array = np.zeros(A_LBtime.shape[0]*7).reshape(A_LBtime.shape[0], 7)
    for i in range(A_LBtime.shape[0]):
        fname="saved_data/"+simulation+"_z"+str(abs(round(redshift_array[i], 2)))+endstring
        names, trash, Az_array[i][:] = np.loadtxt(fname, dtype={'names':('variable_name','equalsign','value'), 'formats':('|S15','|S8', np.float)}, unpack=True)

    simulation="Fire"
    redshift_array=Fz
    Fz_array = np.zeros(F_LBtime.shape[0]*7).reshape(F_LBtime.shape[0], 7)
    for i in range(F_LBtime.shape[0]):
        fname="saved_data/"+simulation+"_z"+str(abs(round(redshift_array[i], 2)))+endstring
        names, trash, Fz_array[i][:] = np.loadtxt(fname, dtype={'names':('variable_name','equalsign','value'), 'formats':('|S15','|S8', np.float)}, unpack=True)



    endstring = "_spinparameters_haloONLY.txt" #this will average over everything in the halo but exclude the innter 0.1Rvir of material

    #create 2-D numpy arrays.  First index is redshift, 2nd index is which lambda parameter:  Xz_array[redshift,lambda]
    simulation="Enzo"
    redshift_array=Ez
    Ez_array_H = np.zeros(E_LBtime.shape[0]*7).reshape(E_LBtime.shape[0], 7)
    for i in range(E_LBtime.shape[0]):
        fname="saved_data/"+simulation+"_z"+str(abs(round(redshift_array[i], 2)))+endstring
        names, trash, Ez_array_H[i][:] = np.loadtxt(fname, dtype={'names':('variable_name','equalsign','value'), 'formats':('|S15','|S8', np.float)}, unpack=True)

    simulation="Ramses"
    redshift_array=Rz
    Rz_array_H = np.zeros(R_LBtime.shape[0]*7).reshape(R_LBtime.shape[0], 7)
    for i in range(R_LBtime.shape[0]):
        fname="saved_data/"+simulation+"_z"+str(abs(round(redshift_array[i], 2)))+endstring
        names, trash, Rz_array_H[i][:] = np.loadtxt(fname, dtype={'names':('variable_name','equalsign','value'), 'formats':('|S15','|S8', np.float)}, unpack=True)

    simulation="Art"
    redshift_array=Artz
    Artz_array_H = np.zeros(Art_LBtime.shape[0]*7).reshape(Art_LBtime.shape[0], 7)
    for i in range(Art_LBtime.shape[0]):
        fname="saved_data/"+simulation+"_z"+str(abs(round(redshift_array[i], 2)))+endstring
        names, trash, Artz_array_H[i][:] = np.loadtxt(fname, dtype={'names':('variable_name','equalsign','value'), 'formats':('|S15','|S8', np.float)}, unpack=True)

    simulation="Arepo"
    redshift_array=Az
    Az_array_H = np.zeros(A_LBtime.shape[0]*7).reshape(A_LBtime.shape[0], 7)
    for i in range(A_LBtime.shape[0]):
        fname="saved_data/"+simulation+"_z"+str(abs(round(redshift_array[i], 2)))+endstring
        names, trash, Az_array_H[i][:] = np.loadtxt(fname, dtype={'names':('variable_name','equalsign','value'), 'formats':('|S15','|S8', np.float)}, unpack=True)

    simulation="Fire"
    redshift_array=Fz
    Fz_array_H = np.zeros(F_LBtime.shape[0]*7).reshape(F_LBtime.shape[0], 7)
    for i in range(F_LBtime.shape[0]):
        fname="saved_data/"+simulation+"_z"+str(abs(round(redshift_array[i], 2)))+endstring
        names, trash, Fz_array_H[i][:] = np.loadtxt(fname, dtype={'names':('variable_name','equalsign','value'), 'formats':('|S15','|S8', np.float)}, unpack=True)

    print Ez_array[:,5]
    print Ez_array_H[:,5]

    #FOR PLOTS RECALL: First index is redshift, 2nd index is which lambda parameter:  Xz_array[redshift,lambda]
    #FOR lambda variables, recall the order is: tot, dark, stars, particles, gas, cold, hot
        
    z_tick_values = [4, 3, 2, 1.5, 1.25, 1.0, 0.9 ]  #max and min of this are limits of the plot.  
    z_tick_value_labels = [' ', '3', '2', '1.5', '1.25', '1', ' ' ]  #max and min of this are limits of the plot.  
    z_tick_locations = cosmo.lookback_time(z_tick_values).value
    LBT_max = cosmo.lookback_time(np.max(z_tick_values)).value
    LBT_min = cosmo.lookback_time(np.min(z_tick_values)).value

    #make a 2x2 shared axes plot of the following 4 individual plots (shown afterwards)
    #fig = plt.figure(figsize=(10,10))
    fig, ((ax1, ax2), (ax3, ax4), (ax5,ax6)) = plt.subplots(3, 2, sharex='col', sharey='row', figsize=(10,15))
    fig.subplots_adjust(hspace=0)
    fig.subplots_adjust(wspace=0)
    ax1.set_ylabel('Spin Parameter $\lambda$')
    ax3.set_ylabel('Spin Parameter $\lambda$')
    ax5.set_ylabel('Spin Parameter $\lambda$')
    ax3.set_xlabel('Lookback Time [Gyr]')
    ax4.set_xlabel('Lookback Time [Gyr]')
    ax5.set_xlabel('Lookback Time [Gyr]')
    ax6.set_xlabel('Lookback Time [Gyr]')
    ax1.set_ylim(0.01, 0.3)
    ax3.set_ylim(0.00, 0.3)
    ax5.set_ylim(0.00, 0.3)
    ax3.set_xlim(LBT_max, LBT_min)
    ax4.set_xlim(LBT_max, LBT_min)
    ax1b = ax1.twiny()
    ax1b.set_xlabel('Redshift')
    ax1b.set_xlim(LBT_max, LBT_min)
    ax1b.set_xticks(z_tick_locations)
    ax1b.set_xticklabels(z_tick_value_labels)
    ax2b = ax2.twiny()
    ax2b.set_xlabel('Redshift')
    ax2b.set_xlim(LBT_max, LBT_min)
    ax2b.set_xticks(z_tick_locations)
    ax2b.set_xticklabels(z_tick_value_labels)
    ax1.plot(E_LBtime,Ez_array[:,1], color='k', linewidth=2)
    ax1.plot(E_LBtime,Ez_array_H[:,2], '--', color='y', linewidth=2)
    #ax1.plot(E_LBtime,Ez_array[:,4], color='m', linewidth=2)
    ax1.plot(E_LBtime,Ez_array_H[:,5], ':', color='b', linewidth=2)
    ax1.plot(E_LBtime,Ez_array_H[:,6], '-.', color='r', linewidth=2)
    ax2.plot(R_LBtime,Rz_array[:,1], color='k', linewidth=2)
    ax2.plot(R_LBtime,Rz_array_H[:,2], '--', color='y', linewidth=2)
    #ax2.plot(R_LBtime,Rz_array[:,4], color='m', linewidth=2)
    ax2.plot(R_LBtime,Rz_array_H[:,5], ':', color='b', linewidth=2)
    ax2.plot(R_LBtime,Rz_array_H[:,6], '-.', color='r', linewidth=2)
    ax3.plot(A_LBtime,Az_array[:,1], color='k', linewidth=2)
    ax3.plot(A_LBtime,Az_array_H[:,2], '--', color='y', linewidth=2)
    #ax3.plot(A_LBtime,Az_array[:,4], color='m', linewidth=2)
    ax3.plot(A_LBtime,Az_array_H[:,5], ':', color='b', linewidth=2)
    ax3.plot(A_LBtime,Az_array_H[:,6], '-.', color='r', linewidth=2)
    ax4.plot(F_LBtime,Fz_array[:,1], color='k', linewidth=2)
    ax4.plot(F_LBtime,Fz_array_H[:,2], '--', color='y', linewidth=2)
    #ax4.plot(F_LBtime,Fz_array[:,4], color='m', linewidth=2)
    ax4.plot(F_LBtime,Fz_array_H[:,5], ':', color='b', linewidth=2)
    ax4.plot(F_LBtime,Fz_array_H[:,6], '-.', color='r', linewidth=2)
    ax5.plot(Art_LBtime,Artz_array[:,1], color='k', linewidth=2)
    ax5.plot(Art_LBtime,Artz_array_H[:,2], '--', color='y', linewidth=2)
    #ax5.plot(Art_LBtime,Artz_array[:,4], color='m', linewidth=2)
    ax5.plot(Art_LBtime,Artz_array_H[:,5], ':', color='b', linewidth=2)
    ax5.plot(Art_LBtime,Artz_array_H[:,6], '-.', color='r', linewidth=2)
    #ax6.plot(Art_LBtime,Artz_array[:,6]*0, '-.', color='k', linewidth=2)
    ax2.legend(("Dark", "Stars", "Cold Gas", "Hot Gas"), frameon=False, loc=2)
    #ax6.legend(("Dark", "Stars", "Cold Gas", "Hot Gas"), frameon=False, loc=2)
    ax1.text(9.0, 0.27, 'Enzo', size='x-large')
    ax2.text(9.0, 0.27, 'Ramses', size='x-large')
    ax3.text(9.0, 0.27, 'Arepo', size='x-large')
    ax4.text(9.5, 0.27, 'Gizmo-PSPH', size='x-large')
    ax5.text(9.0, 0.27, 'Art', size='x-large')
    plt.savefig("ALL_2x2_allspins_vs_time.png", bbox_inches='tight')
    plt.close()


    Ecold_ratio = Ez_array_H[:,5]/Ez_array[:,1]
    Rcold_ratio = Rz_array_H[:,5]/Rz_array[:,1]
    Acold_ratio = Az_array_H[:,5]/Az_array[:,1]
    Fcold_ratio = Fz_array_H[:,5]/Fz_array[:,1]
    Artcold_ratio = Artz_array_H[:,5]/Artz_array[:,1]
    Ehot_ratio = Ez_array_H[:,6]/Ez_array[:,1]
    Rhot_ratio = Rz_array_H[:,6]/Rz_array[:,1]
    Ahot_ratio = Az_array_H[:,6]/Az_array[:,1]
    Fhot_ratio = Fz_array_H[:,6]/Fz_array[:,1]
    Arthot_ratio = Artz_array_H[:,6]/Artz_array[:,1]
    #Plot ratio of cold-gas to DM 
    fig = plt.figure()
    plt.ylabel('$\lambda_{cold}$ / $\lambda_{DM}$')
    plt.ylim(0, 10)
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel('Lookback Time [Gyr]')
    ax1.set_xlim(LBT_max, LBT_min)
    ax2 = ax1.twiny()
    ax2.set_xlabel('Redshift')
    ax2.set_xlim(LBT_max, LBT_min)
    ax2.set_xticks(z_tick_locations)
    ax2.set_xticklabels(z_tick_value_labels)
    ax1.plot(E_LBtime,Ez_array[:,5]/Ez_array[:,1], linewidth=2, color='k')
    ax1.plot(R_LBtime,Rz_array[:,5]/Rz_array[:,1], linewidth=2, color='g')
    ax1.plot(Art_LBtime,Artz_array[:,5]/Artz_array[:,1], linewidth=2, color='b')
    ax1.plot(A_LBtime,Az_array[:,5]/Az_array[:,1], linewidth=2, color='r')
    ax1.plot(F_LBtime,Fz_array[:,5]/Fz_array[:,1], linewidth=2, color='y')
    ax1.plot(E_LBtime,Ez_array[:,6]/Ez_array[:,1], ':', linewidth=2, color='k')
    ax1.plot(R_LBtime,Rz_array[:,6]/Rz_array[:,1], linewidth=2, color='g')
    ax1.plot(Art_LBtime,Artz_array[:,6]/Artz_array[:,1], ':', linewidth=2, color='b')
    ax1.plot(A_LBtime,Az_array[:,6]/Az_array[:,1], ':', linewidth=2, color='r')
    ax1.plot(F_LBtime,Fz_array[:,6]/Fz_array[:,1], ':', linewidth=2, color='y')
    ax1.plot([14,0], [1.0,1.0], linewidth=1, color='k')
    ax1.legend(("Enzo","Ramses", "Art", "Arepo", "Fire"), frameon=False)
    plt.savefig("ALL_spin_ratio_vs_time.png")
    plt.close()


    P.rc('font', family='serif', size=15)
    bins = [0, 0.5, 1, 1.5, 2,2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9]
    
    P.figure()
    P.xlabel("$j_{\\rm cold}(>0.1R_{\\rm vir})$ / $j_{\\rm DM}$", fontsize=20)
    P.ylabel("Frequency", fontsize=20)
    P.xlim(0,8)
    n, bins, patches = P.hist([Ecold_ratio,Artcold_ratio,Rcold_ratio, Acold_ratio, Fcold_ratio], stacked=True, histtype='bar', bins=bins, normed=1)
    P.legend(["Enzo","Art","Ramses","Arepo","Fire"], frameon=False)
    P.tight_layout()
    P.savefig("ALL_coldspin_histogram.png")
    P.close()

    P.figure()
    P.xlim(0,8)
    P.xlabel("$j_{\\rm hot}(>0.1R_{\\rm vir})$ / $j_{\\rm DM}$", fontsize=20)
    P.ylabel("Frequency", fontsize=20)
    n, bins, patches = P.hist([Ehot_ratio,Arthot_ratio,Rhot_ratio, Ahot_ratio, Fhot_ratio], stacked=True, histtype='bar', bins=bins, normed=1)
    P.legend(["Enzo","Art","Ramses","Arepo","Fire"], frameon=False)
    P.tight_layout()
    P.savefig("ALL_hotspin_histogram.png")
    P.close()

    cold=[]
    cold.extend(Ecold_ratio)
    cold.extend(Acold_ratio)
    cold.extend(Rcold_ratio)
    cold.extend(Fcold_ratio)
    cold.extend(Artcold_ratio)
    cold.sort()
    cold=np.array(cold)
    print "median, mean, stdev for cold ratios.....", cold[cold.shape[0]/2], cold.mean(), cold.std()
    hot=[]
    hot.extend(Ehot_ratio)
    hot.extend(Ahot_ratio)
    hot.extend(Rhot_ratio)
    hot.extend(Fhot_ratio)
    hot.extend(Arthot_ratio)
    hot.sort()
    hot=np.array(hot)
    print "median, mean, stdev for hot ratios......", hot[hot.shape[0]/2], hot.mean(), hot.std()


    #do same kinds of histograms, but just for spin param (not ratio)
    bins = np.array(range(21))/50.0    
    P.figure()
    P.xlabel("$\lambda_{cold gas}$")
    P.ylabel("Number")
    P.xlim(0,0.4)
    n, bins, patches = P.hist([Ez_array[:,5],Artz_array[:,5],Rz_array[:,5],Az_array[:,5],Fz_array[:,5]], stacked=True, histtype='bar', bins=bins)
    P.legend(["Enzo","Art","Ramses","Arepo","Fire"], frameon=False)
    P.savefig("ALL_coldspinparam_histogram.png")
    P.close()

    P.figure()
    P.xlabel("$\lambda_{hot gas}$")
    P.ylabel("Number")
    n, bins, patches = P.hist([Ez_array[:,6],Artz_array[:,6],Rz_array[:,6],Az_array[:,6],Fz_array[:,6]], stacked=True, histtype='bar', bins=bins)
    P.legend(["Enzo","Art","Ramses","Arepo","Fire"], frameon=False)
    P.savefig("ALL_hotspinparam_histogram.png")
    P.close()

def velocity_histogram(halo,pf):
    vz = halo['gas','velocity_z']
    vz_stars = halo['PartType4','particle_velocity_z']
    vz_dm = halo['PartType1','particle_velocity_z']

    P.figure()
    P.xlabel("velocity z")
    P.xlim(-1000,1000)
    P.ylabel("Number")
    n, bins, patches = P.hist([vz.in_units('km/s'), vz_stars.in_units('km/s'), vz_dm.in_units('km/s')], stacked=True, histtype='bar', bins=30)
    P.legend(["gas", "stars", "DM"], frameon=False)
    P.savefig("Arepo_test_gas_velocity_histogram.png")
    P.close()
    
    

#reads in sipn parameters from saved files
#order of saved files is: tot, dark, stars, particles, gas, cold, hot
def spinparam_vs_time():
    Ez=[0.00, 0.05, 0.10, 0.20, 0.25, 0.30, 0.40, 0.50, 0.75, 1.00, 3.00]
    Rz=[0.00, 0.25, 0.31, 0.39, 0.50, 0.75, 1.00, 1.50, 2.01, 2.64, 3.00, 3.55, 4.00, 5.06]
    Az=[1.82, 3.00, 4.94]
    Fz=[0.80, 0.90, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.20, 2.50, 2.60, 2.80, 3.00, 3.40, 3.82, 3.90, 4.00, 4.50, 5.00]
    #Fz=[1.70, 1.80, 1.90, 2.00, 2.20, 2.50, 2.60, 2.80, 3.00, 3.40, 3.82, 3.90, 4.00, 4.50, 5.00]
    E_LBtime=cosmo.lookback_time(Ez)
    R_LBtime=cosmo.lookback_time(Rz)
    A_LBtime=cosmo.lookback_time(Az)
    F_LBtime=cosmo.lookback_time(Fz)

    #create 2-D numpy arrays.  First index is redshift, 2nd index is which lambda parameter:  Xz_array[redshift,lambda]
    simulation="Enzo"
    redshift_array=Ez
    Ez_array = np.zeros(E_LBtime.shape[0]*7).reshape(E_LBtime.shape[0], 7)
    for i in range(E_LBtime.shape[0]):
        fname="saved_data/"+simulation+"_z"+str(abs(round(redshift_array[i], 2)))+"_spinparameters.txt"
        names, trash, Ez_array[i][:] = np.loadtxt(fname, dtype={'names':('variable_name','equalsign','value'), 'formats':('|S15','|S8', np.float)}, unpack=True)

    simulation="Ramses"
    redshift_array=Rz
    Rz_array = np.zeros(R_LBtime.shape[0]*7).reshape(R_LBtime.shape[0], 7)
    for i in range(R_LBtime.shape[0]):
        fname="saved_data/"+simulation+"_z"+str(abs(round(redshift_array[i], 2)))+"_spinparameters.txt"
        names, trash, Rz_array[i][:] = np.loadtxt(fname, dtype={'names':('variable_name','equalsign','value'), 'formats':('|S15','|S8', np.float)}, unpack=True)

    simulation="Arepo"
    redshift_array=Az
    Az_array = np.zeros(A_LBtime.shape[0]*7).reshape(A_LBtime.shape[0], 7)
    for i in range(A_LBtime.shape[0]):
        fname="saved_data/"+simulation+"_z"+str(abs(round(redshift_array[i], 2)))+"_spinparameters.txt"
        names, trash, Az_array[i][:] = np.loadtxt(fname, dtype={'names':('variable_name','equalsign','value'), 'formats':('|S15','|S8', np.float)}, unpack=True)

    simulation="Fire"
    redshift_array=Fz
    Fz_array = np.zeros(F_LBtime.shape[0]*7).reshape(F_LBtime.shape[0], 7)
    for i in range(F_LBtime.shape[0]):
        fname="saved_data/"+simulation+"_z"+str(abs(round(redshift_array[i], 2)))+"_spinparameters.txt"
        names, trash, Fz_array[i][:] = np.loadtxt(fname, dtype={'names':('variable_name','equalsign','value'), 'formats':('|S15','|S8', np.float)}, unpack=True)

    #FOR PLOTS RECALL: First index is redshift, 2nd index is which lambda parameter:  Xz_array[redshift,lambda]
    #FOR lambda variables, recall the order is: tot, dark, stars, particles, gas, cold, hot

    z_tick_values = [6, 5, 4, 3, 2, 1.5, 1.3 ]  #max and min of this are limits of the plot.  
    z_tick_value_labels = [' ', '5', '4', '3', '2', '1.5', ' ' ]  #max and min of this are limits of the plot.  
    z_tick_locations = cosmo.lookback_time(z_tick_values).value
    LBT_max = cosmo.lookback_time(np.max(z_tick_values)).value
    LBT_min = cosmo.lookback_time(np.min(z_tick_values)).value

    #plot Lambda total
    fig = plt.figure()
    plt.ylabel('Total Spin Parameter $\lambda$')
    plt.ylim(0, 0.3)
    #plt.semilogy()
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twiny()
    ax1.set_xlabel('Lookback Time [Gyr]')
    ax1.set_xlim(cosmo.age(0).value-0.5, 0.0)
    z_tick_values = [7,5,4,3,2,1,0.5,0]
    z_tick_locations = np.array(cosmo.age(0).value-cosmo.age(z_tick_values).value)
    ax2.set_xlabel('Redshift')
    ax2.set_xlim(cosmo.age(0).value-0.5, 0.0)
    ax2.set_xticks(z_tick_locations)
    ax2.set_xticklabels(z_tick_values)
    ax1.plot(E_LBtime,Ez_array[:,0])
    ax1.plot(R_LBtime,Rz_array[:,0])
    ax1.plot(A_LBtime,Az_array[:,0])
    ax1.plot(F_LBtime,Fz_array[:,0])
    ax1.legend(("Enzo","Ramses", "Arepo", "Fire"), frameon=False)
    plt.savefig("ALL_spin_tot_vs_time.png")
    plt.close()

    #Lambda-dark plot 
    fig = plt.figure()
    plt.ylabel('Dark Matter Particles Spin Parameter $\lambda$')
    plt.ylim(0, 0.3)
    #plt.semilogy()
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twiny()
    ax1.set_xlabel('Lookback Time [Gyr]')
    ax1.set_xlim(cosmo.age(0).value-0.5, 0.0)
    z_tick_values = [7,5,4,3,2,1,0.5,0]
    z_tick_locations = np.array(cosmo.age(0).value-cosmo.age(z_tick_values).value)
    ax2.set_xlabel('Redshift')
    ax2.set_xlim(cosmo.age(0).value-0.5, 0.0)
    ax2.set_xticks(z_tick_locations)
    ax2.set_xticklabels(z_tick_values)
    ax1.plot(E_LBtime,Ez_array[:,1])
    ax1.plot(R_LBtime,Rz_array[:,1])
    ax1.plot(A_LBtime,Az_array[:,1])
    ax1.plot(F_LBtime,Fz_array[:,1])
    ax1.legend(("Enzo","Ramses", "Arepo", "Fire"), frameon=False)
    plt.savefig("ALL_spin_dark_vs_time.png")
    plt.close()

    #Lambda-stars plot 
    fig = plt.figure()
    plt.ylabel('Star Particles Spin Parameter $\lambda$')
    plt.ylim(0, 0.3)
    #plt.semilogy()
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twiny()
    ax1.set_xlabel('Lookback Time [Gyr]')
    ax1.set_xlim(cosmo.age(0).value-0.5, 0.0)
    z_tick_values = [7,5,4,3,2,1,0.5,0]
    z_tick_locations = np.array(cosmo.age(0).value-cosmo.age(z_tick_values).value)
    ax2.set_xlabel('Redshift')
    ax2.set_xlim(cosmo.age(0).value-0.5, 0.0)
    ax2.set_xticks(z_tick_locations)
    ax2.set_xticklabels(z_tick_values)
    ax1.plot(E_LBtime,Ez_array[:,2])
    ax1.plot(R_LBtime,Rz_array[:,2])
    ax1.plot(A_LBtime,Az_array[:,2])
    ax1.plot(F_LBtime,Fz_array[:,2])
    ax1.legend(("Enzo","Ramses", "Arepo", "Fire"), frameon=False)
    plt.savefig("ALL_spin_stars_vs_time.png")
    plt.close()

    #Lambda-gas plot 
    fig = plt.figure()
    plt.ylabel('Gas Spin Parameter $\lambda$')
    plt.ylim(0, 0.3)
    #plt.semilogy()
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twiny()
    ax1.set_xlabel('Lookback Time [Gyr]')
    ax1.set_xlim(cosmo.age(0).value-0.5, 0.0)
    z_tick_values = [7,5,4,3,2,1,0.5,0]
    z_tick_locations = np.array(cosmo.age(0).value-cosmo.age(z_tick_values).value)
    ax2.set_xlabel('Redshift')
    ax2.set_xlim(cosmo.age(0).value-0.5, 0.0)
    ax2.set_xticks(z_tick_locations)
    ax2.set_xticklabels(z_tick_values)
    ax1.plot(E_LBtime,Ez_array[:,4])
    ax1.plot(R_LBtime,Rz_array[:,4])
    ax1.plot(A_LBtime,Az_array[:,4])
    ax1.plot(F_LBtime,Fz_array[:,4])
    ax1.legend(("Enzo","Ramses", "Arepo", "Fire"), frameon=False)
    plt.savefig("ALL_spin_gas_vs_time.png")
    plt.close()

    #Lambda-cold gas plot 
    fig = plt.figure()
    plt.ylabel('Cold Gas Spin Parameter $\lambda$')
    plt.ylim(0, 0.3)
    #plt.semilogy()
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twiny()
    ax1.set_xlabel('Lookback Time [Gyr]')
    ax1.set_xlim(cosmo.age(0).value-0.5, 0.0)
    z_tick_values = [7,5,4,3,2,1,0.5,0]
    z_tick_locations = np.array(cosmo.age(0).value-cosmo.age(z_tick_values).value)
    ax2.set_xlabel('Redshift')
    ax2.set_xlim(cosmo.age(0).value-0.5, 0.0)
    ax2.set_xticks(z_tick_locations)
    ax2.set_xticklabels(z_tick_values)
    ax1.plot(E_LBtime,Ez_array[:,5])
    ax1.plot(R_LBtime,Rz_array[:,5])
    ax1.plot(A_LBtime,Az_array[:,5])
    ax1.plot(F_LBtime,Fz_array[:,5])
    ax1.legend(("Enzo","Ramses", "Arepo", "Fire"), frameon=False)
    plt.savefig("ALL_spin_coldgas_vs_time.png")
    plt.close()

    #Lambda-hot gas plot 
    fig = plt.figure()
    plt.ylabel('Hot Gas Spin Parameter $\lambda$')
    plt.ylim(0, 0.3)
    #plt.semilogy()
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twiny()
    ax1.set_xlabel('Lookback Time [Gyr]')
    ax1.set_xlim(cosmo.age(0).value-0.5, 0.0)
    z_tick_values = [7,5,4,3,2,1,0.5,0]
    z_tick_locations = np.array(cosmo.age(0).value-cosmo.age(z_tick_values).value)
    ax2.set_xlabel('Redshift')
    ax2.set_xlim(cosmo.age(0).value-0.5, 0.0)
    ax2.set_xticks(z_tick_locations)
    ax2.set_xticklabels(z_tick_values)
    ax1.plot(E_LBtime,Ez_array[:,6])
    ax1.plot(R_LBtime,Rz_array[:,6])
    ax1.plot(A_LBtime,Az_array[:,6])
    ax1.plot(F_LBtime,Fz_array[:,6])
    ax1.legend(("Enzo","Ramses", "Arepo", "Fire"), frameon=False)
    plt.savefig("ALL_spin_hotgas_vs_time.png")
    plt.close()


    #Instead of plotting all 4 sims on same axis, one lambda at a time, also look at simulations one at a time and plot all lambdas    
    fig = plt.figure()
    plt.ylabel('Spin Parameter $\lambda$')
    plt.ylim(0, 0.3)
    #plt.semilogy()
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twiny()
    ax1.set_xlabel('Lookback Time [Gyr]')
    ax1.set_xlim(LBT_max, LBT_min)
    ax2.set_xlabel('Redshift')
    ax2.set_xlim(LBT_max, LBT_min)
    ax2.set_xticks(z_tick_locations)
    ax2.set_xticklabels(z_tick_values)
    ax1.plot(E_LBtime,Ez_array[:,1], color='k')
    ax1.plot(E_LBtime,Ez_array[:,2], color='y')
    ax1.plot(E_LBtime,Ez_array[:,4], color='m')
    ax1.plot(E_LBtime,Ez_array[:,5], color='b')
    ax1.plot(E_LBtime,Ez_array[:,6], color='r')
    ax1.legend(("Dark", "Stars", "All Gas", "Cold Gas", "Hot Gas"), frameon=False)
    #ax1.legend(("Dark", "Stars", "Cold Gas", "Hot Gas"), frameon=False)
    plt.text(12.7, 0.27, 'Enzo', size='x-large')
    plt.savefig("Enzo_allspins_vs_time.png")
    plt.close()

    fig = plt.figure()
    plt.ylabel('Spin Parameter $\lambda$')
    plt.ylim(0, 0.3)
    #plt.semilogy()
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twiny()
    ax1.set_xlabel('Lookback Time [Gyr]')
    ax1.set_xlim(LBT_max, LBT_min)
    #z_tick_values = [7,5,4,3,2,1,0.5,0]
    #z_tick_locations = np.array(cosmo.age(0).value-cosmo.age(z_tick_values).value)
    ax2.set_xlabel('Redshift')
    ax2.set_xlim(LBT_max, LBT_min)
    ax2.set_xticks(z_tick_locations)
    ax2.set_xticklabels(z_tick_values)
    ax1.plot(R_LBtime,Rz_array[:,1], color='k')
    ax1.plot(R_LBtime,Rz_array[:,2], color='y')
    ax1.plot(R_LBtime,Rz_array[:,4], color='m')
    ax1.plot(R_LBtime,Rz_array[:,5], color='b')
    ax1.plot(R_LBtime,Rz_array[:,6], color='r')
    ax1.legend(("Dark", "Stars", "All Gas", "Cold Gas", "Hot Gas"), frameon=False)
    plt.text(12.7, 0.27, 'Ramses', size='x-large')
    plt.savefig("Ramses_allspins_vs_time.png")
    plt.close()

    fig = plt.figure()
    plt.ylabel('Spin Parameter $\lambda$')
    plt.ylim(0, 0.3)
    #plt.semilogy()
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twiny()
    ax1.set_xlabel('Lookback Time [Gyr]')
    ax1.set_xlim(LBT_max, LBT_min)
    #z_tick_values = [7,5,4,3,2,1,0.5,0]
    #z_tick_locations = np.array(cosmo.age(0).value-cosmo.age(z_tick_values).value)
    ax2.set_xlabel('Redshift')
    ax2.set_xlim(LBT_max, LBT_min)
    ax2.set_xticks(z_tick_locations)
    ax2.set_xticklabels(z_tick_values)
    ax1.plot(A_LBtime,Az_array[:,1], color='k')
    ax1.plot(A_LBtime,Az_array[:,2], color='y')
    ax1.plot(A_LBtime,Az_array[:,4], color='m')
    ax1.plot(A_LBtime,Az_array[:,5], color='b')
    ax1.plot(A_LBtime,Az_array[:,6], color='r')
    ax1.legend(("Dark", "Stars", "All Gas", "Cold Gas", "Hot Gas"), frameon=False)
    plt.text(12.7, 0.27, 'Arepo', size='x-large')
    plt.savefig("Arepo_allspins_vs_time.png")
    plt.close()

    fig = plt.figure()
    plt.ylabel('Spin Parameter $\lambda$')
    plt.ylim(0, 0.3)
    #plt.semilogy()
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twiny()
    ax1.set_xlabel('Lookback Time [Gyr]')
    ax1.set_xlim(LBT_max, LBT_min)
    #z_tick_values = [7,5,4,3,2,1,0.5,0]
    #z_tick_locations = np.array(cosmo.age(0).value-cosmo.age(z_tick_values).value)
    ax2.set_xlabel('Redshift')
    ax2.set_xlim(LBT_max, LBT_min)
    ax2.set_xticks(z_tick_locations)
    ax2.set_xticklabels(z_tick_values)
    ax1.plot(F_LBtime,Fz_array[:,1], color='k')
    ax1.plot(F_LBtime,Fz_array[:,2], color='y')
    ax1.plot(F_LBtime,Fz_array[:,4], color='m')
    ax1.plot(F_LBtime,Fz_array[:,5], color='b')
    ax1.plot(F_LBtime,Fz_array[:,6], color='r')
    ax1.legend(("Dark", "Stars", "All Gas", "Cold Gas", "Hot Gas"), frameon=False)
    plt.text(12.7, 0.27, 'Fire', size='x-large')
    plt.savefig("Fire_allspins_vs_time.png")
    plt.close()


#reads in sipn parameters from saved files
#order of saved files is: Mtot, Mdark, Mstars, Mgas, Mcold, Mhot, Mstars(gal only), MHI (ENZO ONLY!)
def mass_vs_time():
    Ez=[0.00, 0.05, 0.10, 0.20, 0.25, 0.30, 0.40, 0.50, 0.75, 1.00, 3.00]
    Rz=[0.00, 0.25, 0.31, 0.39, 0.50, 0.75, 1.00, 1.50, 2.01, 2.64, 3.00, 3.55, 4.00, 5.06]
    Az=[1.532, 1.82, 3.00, 4.94]
    Artz=[2.03, 2.99, 4.0]
    Fz=[0.80, 0.90, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.20, 2.50, 2.60, 2.80, 3.00, 3.40, 3.82, 3.90, 4.00, 4.50, 5.00]
    #Fz=[1.70, 1.80, 1.90, 2.00, 2.20, 2.50, 2.60, 2.80, 3.00, 3.40, 3.82, 3.90, 4.00, 4.50, 5.00]
    E_LBtime=cosmo.lookback_time(Ez)
    R_LBtime=cosmo.lookback_time(Rz)
    A_LBtime=cosmo.lookback_time(Az)
    Art_LBtime=cosmo.lookback_time(Artz)
    F_LBtime=cosmo.lookback_time(Fz)

    #create 2-D numpy arrays.  First index is redshift, 2nd index is which lambda parameter:  Xz_array[redshift,Mass]
    simulation="Enzo"
    redshift_array=Ez
    Ez_array = np.zeros(E_LBtime.shape[0]*8).reshape(E_LBtime.shape[0], 8)
    for i in range(E_LBtime.shape[0]):
        fname="saved_data/"+simulation+"_z"+str(abs(round(redshift_array[i], 2)))+"_basicstats.txt"
        names, trash, Ez_array[i][:], trash2 = np.loadtxt(fname, dtype={'names':('variable_name','equalsign','value','units'), 'formats':('|S15','|S8', np.float, '|S15')}, unpack=True)

    simulation="Ramses"
    redshift_array=Rz
    Rz_array = np.zeros(R_LBtime.shape[0]*7).reshape(R_LBtime.shape[0], 7)
    for i in range(R_LBtime.shape[0]):
        fname="saved_data/"+simulation+"_z"+str(abs(round(redshift_array[i], 2)))+"_basicstats.txt"
        names, trash, Rz_array[i][:], trash2 = np.loadtxt(fname, dtype={'names':('variable_name','equalsign','value','units'), 'formats':('|S15','|S8', np.float, '|S15')}, unpack=True)

    simulation="Art"
    redshift_array=Artz
    Artz_array = np.zeros(Art_LBtime.shape[0]*7).reshape(Art_LBtime.shape[0], 7)
    for i in range(Art_LBtime.shape[0]):
        fname="saved_data/"+simulation+"_z"+str(abs(round(redshift_array[i], 2)))+"_basicstats.txt"
        names, trash, Artz_array[i][:], trash2 = np.loadtxt(fname, dtype={'names':('variable_name','equalsign','value','units'), 'formats':('|S15','|S8', np.float, '|S15')}, unpack=True)

    simulation="Arepo"
    redshift_array=Az
    Az_array = np.zeros(A_LBtime.shape[0]*7).reshape(A_LBtime.shape[0], 7)
    for i in range(A_LBtime.shape[0]):
        fname="saved_data/"+simulation+"_z"+str(abs(round(redshift_array[i], 2)))+"_basicstats.txt"
        names, trash, Az_array[i][:], trash2 = np.loadtxt(fname, dtype={'names':('variable_name','equalsign','value','units'), 'formats':('|S15','|S8', np.float, '|S15')}, unpack=True)

    simulation="Fire"
    redshift_array=Fz
    Fz_array = np.zeros(F_LBtime.shape[0]*7).reshape(F_LBtime.shape[0], 7)
    for i in range(F_LBtime.shape[0]):
        fname="saved_data/"+simulation+"_z"+str(abs(round(redshift_array[i], 2)))+"_basicstats.txt"
        names, trash, Fz_array[i][:], trash2 = np.loadtxt(fname, dtype={'names':('variable_name','equalsign','value','units'), 'formats':('|S15','|S8', np.float, '|S15')}, unpack=True)

    #FOR PLOTS RECALL: First index is redshift, 2nd index is which lambda parameter:  Xz_array[redshift,Mass]
    #FOR lambda variables, recall the order is:  Mtot, Mdark, Mstars, Mgas, Mcold, Mhot, Mstars(gal only), MHI (ENZO ONLY!)

    z_tick_values = [4, 3, 2, 1.5, 1.25, 1.0 ]  #max and min of this are limits of the plot.  
    z_tick_value_labels = [' ', '3', '2', '1.5', '1.25', '1.0' ]  
    z_tick_locations = cosmo.lookback_time(z_tick_values).value
    LBT_max = cosmo.lookback_time(np.max(z_tick_values)).value
    LBT_min = cosmo.lookback_time(np.min(z_tick_values)).value

    #make a SINGLE plot for all sims showing Mvir, Mgas, Mhotgas
    plt.rc('font', family='serif', size=15)
    linestyles=[ '-', '-', '-', '.', '.', '.']
    fig = plt.figure()
    plt.ylabel('log Mass ($ < R_{vir}$)', fontsize=20)
    plt.ylim(9, 13.0)
    #plt.semilogy()
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twiny()
    ax1.set_xlabel('Lookback Time [Gyr]')
    ax1.set_xlim(LBT_max, LBT_min)
    ax2.set_xlabel('Redshift')
    ax2.set_xlim(LBT_max, LBT_min)
    ax2.set_xticks(z_tick_locations)
    ax2.set_xticklabels(z_tick_value_labels)
    E1,=ax1.plot(E_LBtime,np.log10(Ez_array[:,0]), '-', color='k', linewidth=2)
    E2,=ax1.plot(E_LBtime,np.log10(Ez_array[:,3]), '-', color='cyan', linewidth=2)
    E3,=ax1.plot(E_LBtime,np.log10(Ez_array[:,5]), '-', color='m', linewidth=2)
    R1,=ax1.plot(R_LBtime,np.log10(Rz_array[:,0]), ':', color='k', linewidth=2)
    ax1.plot(R_LBtime,np.log10(Rz_array[:,3]), ':', color='cyan', linewidth=2)
    ax1.plot(R_LBtime,np.log10(Rz_array[:,5]), ':', color='m', linewidth=2)
    Art1,=ax1.plot(Art_LBtime,np.log10(Artz_array[:,0]), ':.', color='k', linewidth=2)
    ax1.plot(Art_LBtime,np.log10(Artz_array[:,3]), ':.', color='cyan', linewidth=2)
    ax1.plot(Art_LBtime,np.log10(Artz_array[:,5]), ':.', color='m', linewidth=2)
    A1,=ax1.plot(A_LBtime,np.log10(Az_array[:,0]), '--', color='k', linewidth=2)
    ax1.plot(A_LBtime,np.log10(Az_array[:,3]), '--', color='cyan', linewidth=2)
    ax1.plot(A_LBtime,np.log10(Az_array[:,5]), '--', color='m', linewidth=2)
    F1,=ax1.plot(F_LBtime,np.log10(Fz_array[:,0]), '-.', color='k', linewidth=2)
    ax1.plot(F_LBtime,np.log10(Fz_array[:,3]), '-.', color='cyan', linewidth=2)
    ax1.plot(F_LBtime,np.log10(Fz_array[:,5]), '-.', color='m', linewidth=2)
    ax1.text(9.0, 11.7, 'Total', color='k', fontsize=18)
    ax1.text(9.0, 11.0, 'All Gas', color='cyan', fontsize=18)
    ax1.text(9.0, 10.0, 'Hot Gas', color='m', fontsize=18)
    leg1=ax1.legend([E1,Art1], ["Enzo", "Art"], frameon=False, loc=2, fontsize=18)
    leg2=ax1.legend([R1,F1],   ["Ramses", "Gizmo-PSPH"], frameon=False, loc=9, fontsize=18)
    leg3=ax1.legend([A1],      ["Arepo"], frameon=False, loc=1, fontsize=18)
    plt.gca().add_artist(leg1)
    plt.gca().add_artist(leg2)
    plt.tight_layout()
    plt.savefig("ALL_onepanel_mass_vs_time.png")
    plt.close()
    #ax1.plot(E_LBtime,np.log10(Ez_array[:,6]), '-', color='y', linewidth=2)
    #ax1.plot(E_LBtime,np.log10(Ez_array[:,4]), '-', color='b', linewidth=2)
    #ax1.plot(R_LBtime,np.log10(Rz_array[:,6]), ':', color='y', linewidth=2)
    #ax1.plot(R_LBtime,np.log10(Rz_array[:,4]), ':', color='b', linewidth=2)
    #ax1.plot(Art_LBtime,np.log10(Artz_array[:,6]), ':.', color='y', linewidth=2)
    #ax1.plot(Art_LBtime,np.log10(Artz_array[:,4]), ':.', color='b', linewidth=2)
    #ax1.plot(A_LBtime,np.log10(Az_array[:,6]), '--', color='y', linewidth=2)
    #ax1.plot(A_LBtime,np.log10(Az_array[:,4]), '--', color='b', linewidth=2)
    #ax1.plot(F_LBtime,np.log10(Fz_array[:,6]), '-.', color='y', linewidth=2)
    #ax1.plot(F_LBtime,np.log10(Fz_array[:,4]), '-.', color='b', linewidth=2)



    #make a 4x4 shared axes plot of the following 4 individual plots (shown afterwards)
    #fig = plt.figure(figsize=(10,10))
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row', figsize=(10,10))
    fig.subplots_adjust(hspace=0)
    fig.subplots_adjust(wspace=0)
    ax1.set_ylabel('log Mass (< 300 ckpc)')
    ax3.set_ylabel('log Mass (< 300 ckpc)')
    ax3.set_xlabel('Lookback Time [Gyr]')
    ax4.set_xlabel('Lookback Time [Gyr]')
    ax1.set_ylim(9.01, 13)
    ax3.set_ylim(9.00, 13)
    ax3.set_xlim(LBT_max, LBT_min)
    ax4.set_xlim(LBT_max, LBT_min)
    ax1b = ax1.twiny()
    ax1b.set_xlabel('Redshift')
    ax1b.set_xlim(LBT_max, LBT_min)
    ax1b.set_xticks(z_tick_locations)
    ax1b.set_xticklabels(z_tick_value_labels)
    ax2b = ax2.twiny()
    ax2b.set_xlabel('Redshift')
    ax2b.set_xlim(LBT_max, LBT_min)
    ax2b.set_xticks(z_tick_locations)
    ax2b.set_xticklabels(z_tick_value_labels)
    ax1.plot(E_LBtime,np.log10(Ez_array[:,0]), color='k')
    ax1.plot(E_LBtime,np.log10(Ez_array[:,1]), color='k')
    ax1.plot(E_LBtime,np.log10(Ez_array[:,6]), color='y')
    ax1.plot(E_LBtime,np.log10(Ez_array[:,3]), color='m')
    ax1.plot(E_LBtime,np.log10(Ez_array[:,4]), color='b')
    ax1.plot(E_LBtime,np.log10(Ez_array[:,5]), color='r')
    ax2.plot(R_LBtime,np.log10(Rz_array[:,0]), color='k')
    ax2.plot(R_LBtime,np.log10(Rz_array[:,1]), color='k')
    ax2.plot(R_LBtime,np.log10(Rz_array[:,6]), color='y')
    ax2.plot(R_LBtime,np.log10(Rz_array[:,3]), color='m')
    ax2.plot(R_LBtime,np.log10(Rz_array[:,4]), color='b')
    ax2.plot(R_LBtime,np.log10(Rz_array[:,5]), color='r')
    ax3.plot(A_LBtime,np.log10(Az_array[:,0]), color='k')
    ax3.plot(A_LBtime,np.log10(Az_array[:,1]), color='k')
    ax3.plot(A_LBtime,np.log10(Az_array[:,6]), color='y')
    ax3.plot(A_LBtime,np.log10(Az_array[:,3]), color='m')
    ax3.plot(A_LBtime,np.log10(Az_array[:,4]), color='b')
    ax3.plot(A_LBtime,np.log10(Az_array[:,5]), color='r')
    ax4.plot(F_LBtime,np.log10(Fz_array[:,0]), color='k')
    ax4.plot(F_LBtime,np.log10(Fz_array[:,1]), color='k')
    ax4.plot(F_LBtime,np.log10(Fz_array[:,6]), color='y')
    ax4.plot(F_LBtime,np.log10(Fz_array[:,3]), color='m')
    ax4.plot(F_LBtime,np.log10(Fz_array[:,4]), color='b')
    ax4.plot(F_LBtime,np.log10(Fz_array[:,5]), color='r')
    #ax2.legend(("Dark", "Stars", "All Gas", "Cold Gas", "Hot Gas"), frameon=False)
    ax2.legend(("Total", "Dark", "Stars (Gal)", "All Gas", "Cold Gas", "Hot Gas"), frameon=False)
    ax1.text(12.5, 12.5, 'Enzo', size='x-large')
    ax2.text(12.5, 12.5, 'Ramses', size='x-large')
    ax3.text(12.5, 12.5, 'Arepo', size='x-large')
    ax4.text(12.5, 12.5, 'Fire', size='x-large')
    #plt.savefig("ALL_4x4_allspins_vs_time.png")
    plt.savefig("ALL_4x4_mass_vs_time.png")
    plt.close()


    #Enzo
    fig = plt.figure()
    plt.ylabel('log Mass (< 300 ckpc) -- Enzo')
    plt.ylim(9, 13)
    #plt.semilogy()
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twiny()
    ax1.set_xlabel('Lookback Time [Gyr]')
    ax1.set_xlim(cosmo.age(0).value-0.5, 0.0)
    z_tick_values = [7,5,4,3,2,1,0.5,0]
    z_tick_locations = np.array(cosmo.age(0).value-cosmo.age(z_tick_values).value)
    ax2.set_xlabel('Redshift')
    ax2.set_xlim(cosmo.age(0).value-0.5, 0.0)
    ax2.set_xticks(z_tick_locations)
    ax2.set_xticklabels(z_tick_values)
    ax1.plot(E_LBtime,np.log10(Ez_array[:,0]), color='k')
    ax1.plot(E_LBtime,np.log10(Ez_array[:,1]), color='k')
    ax1.plot(E_LBtime,np.log10(Ez_array[:,6]), color='y')
    ax1.plot(E_LBtime,np.log10(Ez_array[:,3]), color='m')
    ax1.plot(E_LBtime,np.log10(Ez_array[:,4]), color='b')
    ax1.plot(E_LBtime,np.log10(Ez_array[:,5]), color='r')
    ax1.legend(("Total", "Dark", "Stars (Gal)", "All Gas", "Cold Gas", "Hot Gas"), frameon=False)
    plt.savefig("Enzo_mass_vs_time.png")
    plt.close()

    #Ramses
    fig = plt.figure()
    plt.ylabel('log Mass (< 300 ckpc) -- Ramses')
    plt.ylim(9, 13)
    #plt.semilogy()
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twiny()
    ax1.set_xlabel('Lookback Time [Gyr]')
    ax1.set_xlim(cosmo.age(0).value-0.5, 0.0)
    z_tick_values = [7,5,4,3,2,1,0.5,0]
    z_tick_locations = np.array(cosmo.age(0).value-cosmo.age(z_tick_values).value)
    ax2.set_xlabel('Redshift')
    ax2.set_xlim(cosmo.age(0).value-0.5, 0.0)
    ax2.set_xticks(z_tick_locations)
    ax2.set_xticklabels(z_tick_values)
    ax1.plot(R_LBtime,np.log10(Rz_array[:,0]), color='k')
    ax1.plot(R_LBtime,np.log10(Rz_array[:,1]), color='k')
    ax1.plot(R_LBtime,np.log10(Rz_array[:,6]), color='y')
    ax1.plot(R_LBtime,np.log10(Rz_array[:,3]), color='m')
    ax1.plot(R_LBtime,np.log10(Rz_array[:,4]), color='b')
    ax1.plot(R_LBtime,np.log10(Rz_array[:,5]), color='r')
    ax1.legend(("Total", "Dark", "Stars (Gal)", "All Gas", "Cold Gas", "Hot Gas"), frameon=False)
    plt.savefig("Ramses_mass_vs_time.png")
    plt.close()

    #Arepo
    fig = plt.figure()
    plt.ylabel('log Mass (< 300 ckpc) -- Arepo')
    plt.ylim(9, 13)
    #plt.semilogy()
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twiny()
    ax1.set_xlabel('Lookback Time [Gyr]')
    ax1.set_xlim(cosmo.age(0).value-0.5, 0.0)
    z_tick_values = [7,5,4,3,2,1,0.5,0]
    z_tick_locations = np.array(cosmo.age(0).value-cosmo.age(z_tick_values).value)
    ax2.set_xlabel('Redshift')
    ax2.set_xlim(cosmo.age(0).value-0.5, 0.0)
    ax2.set_xticks(z_tick_locations)
    ax2.set_xticklabels(z_tick_values)
    ax1.plot(A_LBtime,np.log10(Az_array[:,0]), color='k')
    ax1.plot(A_LBtime,np.log10(Az_array[:,1]), color='k')
    ax1.plot(A_LBtime,np.log10(Az_array[:,6]), color='y')
    ax1.plot(A_LBtime,np.log10(Az_array[:,3]), color='m')
    ax1.plot(A_LBtime,np.log10(Az_array[:,4]), color='b')
    ax1.plot(A_LBtime,np.log10(Az_array[:,5]), color='r')
    ax1.legend(("Total", "Dark", "Stars (Gal)", "All Gas", "Cold Gas", "Hot Gas"), frameon=False)
    plt.savefig("Arepo_mass_vs_time.png")
    plt.close()

    #Ramses
    fig = plt.figure()
    plt.ylabel('log Mass (< 300 ckpc) -- Fire')
    plt.ylim(9, 13)
    #plt.semilogy()
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twiny()
    ax1.set_xlabel('Lookback Time [Gyr]')
    ax1.set_xlim(cosmo.age(0).value-0.5, 0.0)
    z_tick_values = [7,5,4,3,2,1,0.5,0]
    z_tick_locations = np.array(cosmo.age(0).value-cosmo.age(z_tick_values).value)
    ax2.set_xlabel('Redshift')
    ax2.set_xlim(cosmo.age(0).value-0.5, 0.0)
    ax2.set_xticks(z_tick_locations)
    ax2.set_xticklabels(z_tick_values)
    ax1.plot(F_LBtime,np.log10(Fz_array[:,0]), color='k')
    ax1.plot(F_LBtime,np.log10(Fz_array[:,1]), color='k')
    ax1.plot(F_LBtime,np.log10(Fz_array[:,6]), color='y')
    ax1.plot(F_LBtime,np.log10(Fz_array[:,3]), color='m')
    ax1.plot(F_LBtime,np.log10(Fz_array[:,4]), color='b')
    ax1.plot(F_LBtime,np.log10(Fz_array[:,5]), color='r')
    ax1.legend(("Total", "Dark", "Stars (Gal)", "All Gas", "Cold Gas", "Hot Gas"), frameon=False)
    plt.savefig("Fire_mass_vs_time.png")
    plt.close()
