#This is just a script to run the basic analysis code on my sims

#plan
#1) give inputs (an hdf5 file and a rockstar file)
#2) print host statistics
#3) print galaxy statistics
#4) generate a SFH and a MDF for Fe/H
#5) print v/sigma

import sys, re
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

#NOTE: THIS ASSUMES YOU'VE ADDED THIS FOLDER TO YOUR PYTHON PATH
import andrew_hydro_sim_modules

if len(sys.argv) < 4:
    print 'usage: python hydro_sim_script <path to hdf5> <path to rockstar> <path to save figures>'

particle_file = sys.argv[1]
halo_file = sys.argv[2]
plot_dir = sys.argv[3]

print re.split('_',particle_file)

halo_num = re.split('_',particle_file)[1]

print 'analysis for halo {} \n'.format(halo_num)

#First thing is to print the statistics for the host

andrew_hydro_sim_modules.hydro_sim_analysis.Identify_Host(particle_file, halo_file, print_values=True)

#print galaxy statistics

star_age_T, star_FeH, star_parts_mass, gas_mass = andrew_hydro_sim_modules.hydro_sim_analysis.galaxy_statistics(particle_file, halo_file,print_values=True)

np.savetxt(plot_dir+'/'+str(halo_num)+'_star_ages.txt',np.asarray(star_age_T))
np.savetxt(plot_dir+'/'+str(halo_num)+'_star_FeH.txt',np.asarray(star_FeH)) 

# plot SFH and MDF

agebins = np.linspace(0.0,14.0,100)
m_weight = np.ones_like(star_FeH)/float(len(star_FeH))
Z_bins = np.linspace(-6.0,0.5,60)
Z_bins_plot = [(Z_bins[xx]+Z_bins[xx])/2.0 for xx in range(len(Z_bins)-1)]
'''
plt.figure(1,(10,10))
mpl.rcParams['axes.linewidth']=3
plt.yticks(fontsize = 25)
plt.xticks(fontsize = 25)
plt.xlim([0.0,13.7])
plt.xlabel('Time',fontsize = 25)
plt.ylabel('SFH',fontsize=25)
plt.hist(star_age_T,bins=agebins,color='b',alpha=0.5,histtype='step',cumulative=True,normed=True,label='SFH (in Rvir)',linewidth=3)
plt.savefig(plot_dir+'/'+str(halo_num)+'_galaxy_SFH.png')

plt.figure(2,(10,10))
mpl.rcParams['axes.linewidth']=3
plt.yticks(fontsize = 25)
plt.xticks(fontsize = 25)
plt.xlim([-4.0,0.0])
plt.xlabel('[Fe/H]',fontsize=25)
plt.ylabel('df/d[Fe/H]',fontsize=25)
plt.hist(star_FeH,bins=Z_bins,color='b',weights=m_weight,alpha=0.5,histtype='step',linewidth=3)
plt.savefig(plot_dir+'/'+str(halo_num)+'_galaxy_MDF.png')
'''
andrew_hydro_sim_modules.hydro_sim_analysis.report_velocities(particle_file,halo_file,add_dm=False,add_gas=False,add_stars=True,add_low_res=False)
