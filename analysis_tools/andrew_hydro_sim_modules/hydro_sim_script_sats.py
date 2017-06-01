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

#find sats 

sat_ids, clear_rad = andrew_hydro_sim_modules.hydro_sim_analysis.Identify_Host(particle_file,halo_file,print_hi_res_halos=True)

print sat_ids, clear_rad

#print galaxy statistics

for selected_sat_id in sat_ids:
    print 'for satellite listed as '.format(selected_sat_id)

    star_age_T, star_FeH = andrew_hydro_sim_modules.hydro_sim_analysis.galaxy_statistics(particle_file, halo_file, print_values=True,halo_id=selected_sat_id)

    np.savetxt(plot_dir+'/'+str(halo_num)+'_'+str(selected_sat_id)+'_star_ages.txt',np.asarray(star_age_T))
    np.savetxt(plot_dir+'/'+str(halo_num)+'_'+str(selected_sat_id)+'_star_FeH.txt',np.asarray(star_FeH)) 

    andrew_hydro_sim_modules.hydro_sim_analysis.report_velocities(particle_file,halo_file,add_dm=False,add_gas=False,add_stars=True,add_low_res=False,halo_id=selected_sat_id)
