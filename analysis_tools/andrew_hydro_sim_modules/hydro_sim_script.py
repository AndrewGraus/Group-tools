#This is just a script to run the basic analysis code on my sims

#plan
#1) give inputs (an hdf5 file and a rockstar file)
#2) print host statistics
#3) print galaxy statistics
#4) generate a SFH and a MDF for Fe/H
#5) print v/sigma

import sys
import matplotlib.pyplot as plt
import numpy as np

#NOTE: THIS ASSUMES YOU'VE ADDED THIS FOLDER TO YOUR PYTHON PATH
import andrew_hydro_sim_modules

if len(sys.argv) < 4:
    print 'usage: python hydro_sim_script <path to hdf5> <path to rockstar> <path to save figures>'

particle_file = sys.argv[1]
halo_file = sys.argv[2]
plot_dir = sys.argv[3]
#First thing is to print the statistics for the host

andrew_hydro_sim_modules.Identify_Host(particle_file, halo_file, print_values=True)

#print galaxy statistics

star_age_T, star_FeH = andrew_hydro_sim_modules.galaxy_statistics(particle_file, halo_file, print_values=True)

# plot SFH and MDF

agebins = np.linspace(0.0,14.0,100)
m_weight = np.ones_like(gal_FeH)/float(len(gal_FeH))
Z_bins = np.linspace(-6.0,0.5,60)
Z_bins_plot = [(Z_bins[xx]+Z_bins[xx])/2.0 for xx in range(len(Z_bins)-1)]

plt.figure(1,(10,10))
rc('axes',linewidth=3)
rc('axes',linewidth=3)
plt.yticks(fontsize = 25)
plt.xticks(fontsize = 25)
plt.xlim([0.0,13.7])
plt.xlabel('Time',fontsize = 25)
plt.ylabel('SFH',fontsize=25)
plt.hist(gal_T,bins=agebins,color='b',alpha=0.5,histtype='step',cumulative=True,normed=True,label='SFH (in Rvir)',linewidth=3)
plt.savefig(plot_dir+'galaxy_SFH.py')

plt.figure(2,(10,10))
rc('axes',linewidth=3)
rc('axes',linewidth=3)
plt.yticks(fontsize = 25)
plt.xticks(fontsize = 25)
plt.xlim([-4.0,0.0])
plt.xlabel('[Fe/H]',fontsize=25)
plt.ylabel('df/d[Fe/H]',fontsize=25)
plt.hist(gal_FeH,bins=Z_bins,color='b',weights=m_weight,alpha=0.5,histtype='step',linewidth=3)
plt.savefig(plot_dir+'galaxy_MDF.py')

#print V_sigma

andrew_hydro_sim_modules.report_velocities(hdf5_file,halo_file,add_dm=False,add_gas=False,add_stars=True,add_low_res=False)
