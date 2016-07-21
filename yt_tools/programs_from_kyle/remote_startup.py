screen

qsub -I
ipython


import matplotlib
matplotlib.use('Agg')  #so you don't need an X-server
from matplotlib import pyplot as plt
import yt
from ytstewart import *


sim='arepo'
#for redshift in [1.53, 1.82, 2,3,4]:
for redshift in [1.53, 1.82, 2,3]:
    halo,pf=loadsim(sim,redshift)
    #get_halo_parameters(halo,pf)
    #paper1plots(halo,pf)
    #print_basic_stats(halo,pf)
    spin_parameters_outerhalo(halo,pf)
    
sim='art'
#for redshift in [2.03, 2.33, 2.99, 4.0]:
for redshift in [2.03, 2.33, 2.99]:
    #value, center = pf.h.find_max('all_density')
    #rvir=300.0*kpc*(1./(1.+pf.current_redshift))
    #halo    = pf.h.sphere(center, rvir)
    halo,pf=loadsim(sim,redshift)
    #center = halo.quantities.center_of_mass(use_particles=True)
    #halo    = pf.h.sphere(center, rvir)
    density_projections(halo,pf,'Halo')
    get_halo_parameters(halo,pf)
    print_basic_stats(halo,pf)
    spin_parameters(halo,pf)
    spin_parameters_outerhalo(halo,pf)
    paper1plots(halo,pf)

sim='ramses'
#for redshift in [0,0.25,0.32,0.39,0.5,0.75,1,1.5,2,2.64,3,3.55,4,5]:
for redshift in [1,1.5,2,2.64,3]:
    halo,pf=loadsim(sim,redshift)
    spin_parameters_outerhalo(halo,pf)
    
sim='fire'
#for redshift in [0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.2,2.5,2.6,2.8,3,3.4,3.82,3.9,4,4.5,5]:
for redshift in [1.0,1.2,1.4,1.6,1.8,2,2.2,2.5,2.6,2.8,3]:
    halo,pf=loadsim(sim,redshift)
    spin_parameters_outerhalo(halo,pf)
    

sim='enzo'
#for redshift in [0,0.05,0.1,0.2,0.25,0.3,0.4,0.5,0.75,1,3]:
for redshift in [1,3]:
    halo,pf=loadsim(sim,redshift)
    paper1plots(halo,pf)

#redo calculation after fixing code
qsub -I
ipython
import matplotlib
matplotlib.use('Agg')  #so you don't need an X-server
from matplotlib import pyplot as plt
import yt
import ytstewart
from ytstewart import *

import ytstewart
reload(ytstewart)
from ytstewart import *
if sim=='enzo': redshift_list=[1,3]
if sim=='ramses': redshift_list=[1,2,3]
if sim=='arepo': redshift_list=[2,3]
if sim=='fire': redshift_list=[1,2,3]
for redshift in redshift_list:
    halo,pf=loadsim(sim,redshift)
    ytstewart.paper1plots(halo,pf)
    


get_halo_parameters(halo,pf)
print_basic_stats(halo,pf)
spin_parameters(halo,pf)


halo2 = pf.h.sphere(halo.center, halo.radius/2)
center = halo2.quantities.center_of_mass(use_particles=True)
halo2 = pf.h.sphere(center, halo.radius/4)
center = halo2.quantities.center_of_mass(use_particles=True)
halo2 = pf.h.sphere(center, halo.radius/6)
center = halo2.quantities.center_of_mass(use_particles=True)
halo2 = pf.h.sphere(center, halo.radius/8)
center = halo2.quantities.center_of_mass(use_particles=True)
halo2 = pf.h.sphere(center, halo.radius/10)
center = halo2.quantities.center_of_mass(use_particles=True)
halo = pf.h.sphere(center, halo.radius)

    density_projections(halo,pf,'Halo')
    print_basic_stats(halo,pf)
    halo,pf=loadsim(sim,redshift)
    get_halo_parameters(halo,pf)
    paper1plots(halo,pf)
    spin_parameters(halo,pf)
    spin_parameters_outerhalo(halo,pf)
