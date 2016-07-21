from yt.mods import *
from yt.units import kpc, Msun, km, s, cm, g, Gyr, yr
from yt.visualization.plot_window import *
from yt.analysis_modules.halo_finding.api import *
import numpy as np
import scipy
import scipy.integrate
import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.279)
from cosmic_age_estimate import *
import os.path
from get_snapshot_name import get_snapshot_name as getsnap
from load_simulation import load_simulation as loadsim
from scylla_derived_fields import *

def basic_stats_vs_z(simulation,redshifts):
    for z in redshifts:
        halo,pf=loadsim(simulation,z)
        print_basic_stats(halo,pf)
                        
def make_all_plots(halo,pf):
    density_slices_xyz(halo,pf)
    density_projections(halo,pf)
    phase_plots(halo, pf)
    velocity_slices(halo,pf)
    cold_gas_plots(halo,pf)
    cold_dense_gas_plots(halo,pf)
    SFR2(halo,pf)

#figure out what kind of simulation you're looking at and return the simulation name (as a string)
def get_simulation_type(pf):
    if "Ramses" in pf.directory or "ramses" in pf.directory:
        simulation="Ramses"
    elif "Enzo" in pf.directory or "enzo" in pf.directory:
        simulation="Enzo"
    elif "Arepo" in pf.directory or "arepo" in pf.directory:
        simulation="Arepo"
    elif "Fire" in pf.directory or "fire" in pf.directory or "Gadget" in pf.directory or "gadget" in pf.directory:
        simulation="Fire"
    elif "Art" in pf.directory or "art" in pf.directory:
        simulation="Art"
    else:
        simulation="UNDETERMINED"
    return simulation

def get_halo_parameters(halo,pf):
    simulation=get_simulation_type(pf)
    simulation2=simulation+"_z"+str(abs(round(pf.current_redshift, 2)))
    rvir = pf.arr(300.0, 'kpccm')
    #halo    = pf.h.sphere(center, rvir)
    center = halo.get_field_parameter("center").in_units('code_length')
    gal = pf.h.sphere(center, rvir/10)
    if simulation=="Arepo":
        bv = gal.quantities.bulk_velocity(use_particles=False).in_units('cm/s') #Arepo gives weird errors if you try using particles....
        Lgal = gal.quantities.angular_momentum_vector(use_particles=False)
    else:
        bv = gal.quantities.bulk_velocity(use_particles=True).in_units('cm/s')
        Lgal = gal.quantities.angular_momentum_vector(use_particles=True)
    halo.set_field_parameter("bulk_velocity", bv)
    Lgal = YTArray( [Lgal[0].convert_to_cgs().value, Lgal[1].convert_to_cgs().value, Lgal[2].convert_to_cgs().value], 'cm**2/s')
    halo.set_field_parameter("angular_momentum_vector", Lgal)
    print simulation2 + "  bulk velocity: bv=", bv
    print simulation2 + "  Lgalaxy: Lgal=", Lgal

    f = open( "saved_data/"+simulation2+"_haloparameters.txt", "w")
    f.write(  "center(code_length)              = %.6g" % center[0].value + " %.6g" % center[1].value + " %.6g" % center[2].value)
    f.write("\nbulk_velocity(cm/s)              = %.6g" % bv[0].value     + " %.6g" % bv[1].value     + " %.6g" % bv[2].value)
    f.write("\nangular_momentum_vector(cm**2/s) = %.6g" % Lgal[0].value   + " %.6g" % Lgal[1].value   + " %.6g" % Lgal[2].value)
    f.close()

#Use equation from Bullock '05 to calculate spin parameters
# lambda = j / sqrt(2)*Mvir*Vvir, where Vvir is sqrt(GM/R) at virial radius
def spin_parameters(halo, pf):
    from yt.utilities.physical_constants import G
    simulation=get_simulation_type(pf)
    simulation2=simulation+"_z"+str(abs(round(pf.current_redshift, 2)))
    rvir=halo.get_field_parameter('radius').in_units('kpc')
    #halo    = pf.h.sphere(center, rvir)
    center = halo.get_field_parameter("center")
    
    #determine the (mass-weighted) specific angular momentum of various things (ALL, DM, stars, allgas, hotgas, coldgas)
    #aboslutely everything (built in routine)
    j_vector = halo.quantities.angular_momentum_vector()
    jall = np.sqrt(j_vector[0]**2+j_vector[1]**2 + j_vector[2]**2)

    #all gas
    Mgas  = (halo['gas','cell_mass']).in_units('Msun')
    jx = ( (Mgas*halo['gas','specific_angular_momentum_x']).sum() )/( Mgas.sum() )
    jy = ( (Mgas*halo['gas','specific_angular_momentum_y']).sum() )/( Mgas.sum() )
    jz = ( (Mgas*halo['gas','specific_angular_momentum_z']).sum() )/( Mgas.sum() )
    jgas = np.sqrt(jx*jx+jy*jy+jz*jz)

    #cold gas and hot gas
    Temp  = halo['gas','temperature'].in_units('K')
    cold = Temp < 2.5e5
    hot  = Temp > 2.5e5
    Mcold = ( halo['gas','cell_mass'] )[cold].in_units('Msun')
    jx = ( (Mcold*(halo['gas','specific_angular_momentum_x'])[cold]).sum() )/( Mcold.sum() )
    jy = ( (Mcold*(halo['gas','specific_angular_momentum_y'])[cold]).sum() )/( Mcold.sum() )
    jz = ( (Mcold*(halo['gas','specific_angular_momentum_z'])[cold]).sum() )/( Mcold.sum() )
    jcold = np.sqrt(jx*jx+jy*jy+jz*jz)
    Mhot = ( halo['gas','cell_mass'] )[hot].in_units('Msun')
    jx = ( (Mhot*(halo['gas','specific_angular_momentum_x'])[hot]).sum() )/( Mhot.sum() )
    jy = ( (Mhot*(halo['gas','specific_angular_momentum_y'])[hot]).sum() )/( Mhot.sum() )
    jz = ( (Mhot*(halo['gas','specific_angular_momentum_z'])[hot]).sum() )/( Mhot.sum() )
    jhot = np.sqrt(jx*jx+jy*jy+jz*jz)

    #stars, DM, and all particles (stars+DM)
    if simulation=="Enzo":
        ct    = halo["io","creation_time"]
        stars = (ct > 0)
        DM    = (ct <= 0)
    if simulation=="Ramses":
        age = halo["io","particle_age"]
        stars = (age < 0)
        DM = (age == 0)
    if simulation=="Enzo" or simulation=="Ramses":
        Mstars = (halo['all','particle_mass'])[stars].in_units('Msun')
        jx = ( (Mstars*(halo['all','particle_specific_angular_momentum_x'])[stars]).sum() )/( Mstars.sum() )
        jy = ( (Mstars*(halo['all','particle_specific_angular_momentum_y'])[stars]).sum() )/( Mstars.sum() )
        jz = ( (Mstars*(halo['all','particle_specific_angular_momentum_z'])[stars]).sum() )/( Mstars.sum() )
        jstars = np.sqrt(jx*jx+jy*jy+jz*jz)

        Mdark = (halo['all','particle_mass'])[DM].in_units('Msun')
        jx = ( (Mdark*(halo['all','particle_specific_angular_momentum_x'])[DM]).sum() )/( Mdark.sum() )
        jy = ( (Mdark*(halo['all','particle_specific_angular_momentum_y'])[DM]).sum() )/( Mdark.sum() )
        jz = ( (Mdark*(halo['all','particle_specific_angular_momentum_z'])[DM]).sum() )/( Mdark.sum() )
        jdark = np.sqrt(jx*jx+jy*jy+jz*jz)
        
    if simulation=="Arepo" or simulation=="Fire":
        #in these simulations "PartType4" is stars and "PartType1" is DM (in halo region at least)
        Mstars = ( halo['PartType4', 'particle_mass'] ).in_units('Msun')
        jx = ( (Mstars*halo['PartType4', 'particle_specific_angular_momentum_x'] ).sum() )/( Mstars.sum() )
        jy = ( (Mstars*halo['PartType4', 'particle_specific_angular_momentum_y'] ).sum() )/( Mstars.sum() )
        jz = ( (Mstars*halo['PartType4', 'particle_specific_angular_momentum_z'] ).sum() )/( Mstars.sum() )
        jstars = np.sqrt(jx*jx+jy*jy+jz*jz)

        Mdark = ( halo[('PartType1', 'particle_mass')] ).in_units('Msun')
        jx = ( (Mdark* halo['PartType1','particle_specific_angular_momentum_x'] ).sum() )/( Mdark.sum() )
        jy = ( (Mdark* halo['PartType1','particle_specific_angular_momentum_y'] ).sum() )/( Mdark.sum() )
        jz = ( (Mdark* halo['PartType1','particle_specific_angular_momentum_z'] ).sum() )/( Mdark.sum() )
        jdark = np.sqrt(jx*jx+jy*jy+jz*jz)

    if simulation=="Art":
        #in these simulations "PartType4" is stars and "PartType1" is DM (in halo region at least)
        Mstars = ( halo['stars', 'particle_mass'] ).in_units('Msun')
        jx = ( (Mstars*halo['stars', 'particle_specific_angular_momentum_x'] ).sum() )/( Mstars.sum() )
        jy = ( (Mstars*halo['stars', 'particle_specific_angular_momentum_y'] ).sum() )/( Mstars.sum() )
        jz = ( (Mstars*halo['stars', 'particle_specific_angular_momentum_z'] ).sum() )/( Mstars.sum() )
        jstars = np.sqrt(jx*jx+jy*jy+jz*jz)

        Mdark = ( halo[('darkmatter', 'particle_mass')] ).in_units('Msun')
        jx = ( (Mdark* halo['darkmatter','particle_specific_angular_momentum_x'] ).sum() )/( Mdark.sum() )
        jy = ( (Mdark* halo['darkmatter','particle_specific_angular_momentum_y'] ).sum() )/( Mdark.sum() )
        jz = ( (Mdark* halo['darkmatter','particle_specific_angular_momentum_z'] ).sum() )/( Mdark.sum() )
        jdark = np.sqrt(jx*jx+jy*jy+jz*jz)

    #combine info from stars and dark matter to get total "particle" angmom
    jpart = (jstars*Mstars.sum() + jdark*Mdark.sum()) / (Mstars.sum()+Mdark.sum())

    #determine the common denominator (for all different calculations)
    Mvir = Mdark.sum() + Mstars.sum() + Mgas.sum()
    Rvir = halo.get_field_parameter('radius').in_units('kpc')
    Vvir = np.sqrt(G*Mvir/Rvir).in_units('km/s')
    denom = (np.sqrt(2)*Rvir*Vvir).in_units('cm**2/s')

    lambda_all   = jall  /denom
    lambda_dark  = jdark /denom
    lambda_gas   = jgas  /denom
    lambda_cold  = jcold /denom
    lambda_hot   = jhot  /denom
    lambda_stars = jstars/denom
    lambda_part  = jpart /denom
    
    f = open( "saved_data/"+simulation2+"_spinparameters.txt", "w")
    f.write(  "Lambda_tot   = " "%.4g " % lambda_all)
    f.write("\nLambda_dark  = " "%.4g " % lambda_dark)
    f.write("\nLambda_stars = " "%.4g " % lambda_stars)
    f.write("\nLambda_part. = " "%.4g " % lambda_part)
    f.write("\nLambda_gas   = " "%.4g " % lambda_gas)
    f.write("\nLambda_cold  = " "%.4g " % lambda_cold)
    f.write("\nLambda_hot   = " "%.4g " % lambda_hot)
    f.close()


#Use equation from Bullock '05 to calculate spin parameters
# lambda = j / sqrt(2)*Mvir*Vvir, where Vvir is sqrt(GM/R) at virial radius
#same as previous routine, but EXCLUDE central galaxy (<Rvir/10)
def spin_parameters_outerhalo(halo, pf):
    from yt.utilities.physical_constants import G
    simulation=get_simulation_type(pf)
    simulation2=simulation+"_z"+str(abs(round(pf.current_redshift, 2)))
    rvir=halo.get_field_parameter('radius').in_units('kpc')
    #halo    = pf.h.sphere(center, rvir)
    center = halo.get_field_parameter("center")
    bv = halo.get_field_parameter('bulk_velocity')

    #sub-select only gas outside the galactic region
    gashalo=halo.cut_region(["(obj[('gas','radius')]).in_units('kpccm') > 30.0"])  #cut out the central 0.1Rvir
    
    #determine the (mass-weighted) specific angular momentum of various things (ALL, DM, stars, allgas, hotgas, coldgas)
    #aboslutely everything (built in routine)
    #j_vector = halo.quantities.angular_momentum_vector()
    #jall = np.sqrt(j_vector[0]**2+j_vector[1]**2 + j_vector[2]**2)

    #all gas
    Mgas  = (gashalo['gas','cell_mass']).in_units('Msun')
    jx = ( (Mgas*gashalo['gas','specific_angular_momentum_x']).sum() )/( Mgas.sum() )
    jy = ( (Mgas*gashalo['gas','specific_angular_momentum_y']).sum() )/( Mgas.sum() )
    jz = ( (Mgas*gashalo['gas','specific_angular_momentum_z']).sum() )/( Mgas.sum() )
    jgas = np.sqrt(jx*jx+jy*jy+jz*jz)

    #cold gas and hot gas
    Temp  = gashalo['gas','temperature'].in_units('K')
    cold = Temp < 2.5e5
    hot  = Temp > 2.5e5
    Mcold = ( gashalo['gas','cell_mass'] )[cold].in_units('Msun')
    jx = ( (Mcold*(gashalo['gas','specific_angular_momentum_x'])[cold]).sum() )/( Mcold.sum() )
    jy = ( (Mcold*(gashalo['gas','specific_angular_momentum_y'])[cold]).sum() )/( Mcold.sum() )
    jz = ( (Mcold*(gashalo['gas','specific_angular_momentum_z'])[cold]).sum() )/( Mcold.sum() )
    jcold = np.sqrt(jx*jx+jy*jy+jz*jz)
    Mhot = ( gashalo['gas','cell_mass'] )[hot].in_units('Msun')
    jx = ( (Mhot*(gashalo['gas','specific_angular_momentum_x'])[hot]).sum() )/( Mhot.sum() )
    jy = ( (Mhot*(gashalo['gas','specific_angular_momentum_y'])[hot]).sum() )/( Mhot.sum() )
    jz = ( (Mhot*(gashalo['gas','specific_angular_momentum_z'])[hot]).sum() )/( Mhot.sum() )
    jhot = np.sqrt(jx*jx+jy*jy+jz*jz)

    #stars, DM, and all particles (stars+DM)
    if simulation=="Enzo" or simulation=="Ramses":
        parthalo=halo.cut_region(["(obj[('all','particle_radius')]).in_units('kpccm') > 30.0"])  #cut out the central 0.1Rvir
        if simulation=="Enzo":
            ct    = parthalo['io','creation_time']
            stars = (ct > 0)
            DM    = (ct <= 0)
        if simulation=="Ramses":
            age = parthalo['io','particle_age']
            stars = (age < 0)
            DM = (age == 0)
        Mstars = (parthalo['all','particle_mass'])[stars].in_units('Msun')
        jx = ( (Mstars*(parthalo['all','particle_specific_angular_momentum_x'])[stars]).sum() )/( Mstars.sum() )
        jy = ( (Mstars*(parthalo['all','particle_specific_angular_momentum_y'])[stars]).sum() )/( Mstars.sum() )
        jz = ( (Mstars*(parthalo['all','particle_specific_angular_momentum_z'])[stars]).sum() )/( Mstars.sum() )
        jstars = np.sqrt(jx*jx+jy*jy+jz*jz)

        Mdark = (parthalo['all','particle_mass'])[DM].in_units('Msun')
        jx = ( (Mdark*(parthalo['all','particle_specific_angular_momentum_x'])[DM]).sum() )/( Mdark.sum() )
        jy = ( (Mdark*(parthalo['all','particle_specific_angular_momentum_y'])[DM]).sum() )/( Mdark.sum() )
        jz = ( (Mdark*(parthalo['all','particle_specific_angular_momentum_z'])[DM]).sum() )/( Mdark.sum() )
        jdark = np.sqrt(jx*jx+jy*jy+jz*jz)
        
    if simulation=="Arepo" or simulation=="Fire":
        #in these simulations "PartType4" is stars and "PartType1" is DM (in halo region at least)
        #parthalo=halo.cut_region(["np.logical_or( (obj[('PartType4','particle_radius')]).in_units('kpccm') > 30, obj[('PartType1','particle_radius')].in_units('kpccm') > 30)"])
        darkhalo=halo.cut_region(["(obj[('PartType1','particle_radius')].in_units('kpccm') > 30)"])
        starhalo=halo.cut_region(["(obj[('PartType4','particle_radius')]).in_units('kpccm') > 30"])
        
        Mstars = ( starhalo['PartType4', 'particle_mass'] ).in_units('Msun')
        jx = ( (Mstars*starhalo['PartType4', 'particle_specific_angular_momentum_x'] ).sum() )/( Mstars.sum() )
        jy = ( (Mstars*starhalo['PartType4', 'particle_specific_angular_momentum_y'] ).sum() )/( Mstars.sum() )
        jz = ( (Mstars*starhalo['PartType4', 'particle_specific_angular_momentum_z'] ).sum() )/( Mstars.sum() )
        jstars = np.sqrt(jx*jx+jy*jy+jz*jz)

        Mdark = ( darkhalo[('PartType1', 'particle_mass')] ).in_units('Msun')
        jx = ( (Mdark* darkhalo['PartType1','particle_specific_angular_momentum_x'] ).sum() )/( Mdark.sum() )
        jy = ( (Mdark* darkhalo['PartType1','particle_specific_angular_momentum_y'] ).sum() )/( Mdark.sum() )
        jz = ( (Mdark* darkhalo['PartType1','particle_specific_angular_momentum_z'] ).sum() )/( Mdark.sum() )
        jdark = np.sqrt(jx*jx+jy*jy+jz*jz)

    if simulation=="Art":
        #in these simulations "stars" is stars and "darkmatter" is DM 
        darkhalo=halo.cut_region(["(obj[('darkmatter','particle_radius')].in_units('kpccm') > 30)"])
        starhalo=halo.cut_region(["(obj[('stars','particle_radius')]).in_units('kpccm') > 30"])
        Mstars = ( starhalo['stars', 'particle_mass'] ).in_units('Msun')
        jx = ( (Mstars*starhalo['stars', 'particle_specific_angular_momentum_x'] ).sum() )/( Mstars.sum() )
        jy = ( (Mstars*starhalo['stars', 'particle_specific_angular_momentum_y'] ).sum() )/( Mstars.sum() )
        jz = ( (Mstars*starhalo['stars', 'particle_specific_angular_momentum_z'] ).sum() )/( Mstars.sum() )
        jstars = np.sqrt(jx*jx+jy*jy+jz*jz)

        Mdark = ( darkhalo[('darkmatter', 'particle_mass')] ).in_units('Msun')
        jx = ( (Mdark* darkhalo['darkmatter','particle_specific_angular_momentum_x'] ).sum() )/( Mdark.sum() )
        jy = ( (Mdark* darkhalo['darkmatter','particle_specific_angular_momentum_y'] ).sum() )/( Mdark.sum() )
        jz = ( (Mdark* darkhalo['darkmatter','particle_specific_angular_momentum_z'] ).sum() )/( Mdark.sum() )
        jdark = np.sqrt(jx*jx+jy*jy+jz*jz)

    #combine info from stars and dark matter to get total "particle" angmom
    jpart = (jstars*Mstars.sum() + jdark*Mdark.sum()) / (Mstars.sum()+Mdark.sum())
    jall =  (jstars*Mstars.sum() + jdark*Mdark.sum() + jgas*Mgas.sum()) / (Mstars.sum()+Mdark.sum()+Mgas.sum())

    #determine the common denominator (for all different calculations)
    Mvir = Mdark.sum() + Mstars.sum() + Mgas.sum()
    Rvir = halo.get_field_parameter('radius').in_units('kpc')
    Vvir = np.sqrt(G*Mvir/Rvir).in_units('km/s')
    denom = (np.sqrt(2)*Rvir*Vvir).in_units('cm**2/s')

    lambda_all   = jall  /denom
    lambda_dark  = jdark /denom
    lambda_gas   = jgas  /denom
    lambda_cold  = jcold /denom
    lambda_hot   = jhot  /denom
    lambda_stars = jstars/denom
    lambda_part  = jpart /denom
    
    f = open( "saved_data/"+simulation2+"_spinparameters_haloONLY.txt", "w")
    f.write(  "Lambda_tot   = " "%.4g " % lambda_all)
    f.write("\nLambda_dark  = " "%.4g " % lambda_dark)
    f.write("\nLambda_stars = " "%.4g " % lambda_stars)
    f.write("\nLambda_part. = " "%.4g " % lambda_part)
    f.write("\nLambda_gas   = " "%.4g " % lambda_gas)
    f.write("\nLambda_cold  = " "%.4g " % lambda_cold)
    f.write("\nLambda_hot   = " "%.4g " % lambda_hot)
    f.close()

def print_basic_stats(halo, pf):
    simulation=get_simulation_type(pf)
    simulation2=simulation+"_z"+str(abs(round(pf.current_redshift, 2)))
    rvir=halo.get_field_parameter('radius').in_units('kpc')
    #halo    = pf.h.sphere(center, rvir)
    center = halo.get_field_parameter("center")
    Rgal=20.0*kpc
    gal = pf.h.sphere(center, Rgal)

    #print out some basic statistisc for the halo within <300 comoving kpc.  Mvir, Mcold, Mhot, etc.
    #Mass  = (halo["particle_mass"].sum() + halo["cell_mass"].sum()).in_units('Msun')
    #Mgas, Mparticle  = halo.quantities.total_mass().in_units('Msun')
    Mgas  = (halo['gas','cell_mass'].sum()).in_units('Msun')
    Mparticles  = (halo['all','particle_mass'].sum()).in_units('Msun')
    Mass = Mgas + Mparticles
    if simulation=="Enzo":
        ct    = halo["creation_time"]
        stars = (ct > 0)
        DM    = (ct <= 0)
    if simulation=="Ramses":
        age = halo["io","particle_age"]
        stars = (age < 0)
        DM = (age == 0)

    #Mgas =  ( halo["cell_mass"].sum() ).in_units('Msun')
    #Mstardark = halo["particle_mass"].sum()
    if simulation=="Enzo":
        MHI=( halo['H_mass'].sum() ).in_units('Msun')

    if simulation=="Arepo" or simulation=="Fire":
        Mstar = ( halo[('PartType4', 'Masses')].sum() ).in_units('Msun')
        #Mdark = halo["particle_mass"].sum().in_units('Msun') - Mstar
        Mdark = Mparticles - Mstar 
    if simulation=="Enzo" or simulation=="Ramses":
        Mstar = ( (halo["particle_mass"])[stars].sum() ).in_units('Msun')
        #Mdark = ( (halo["particle_mass"])[DM].sum()    ).in_units('Msun')
        Mdark = Mparticles - Mstar 
    if simulation=="Art":
        Mstar = ( (halo["stars","particle_mass"]).sum() ).in_units('Msun')
        Mdark = ( (halo["darkmatter","particle_mass"]).sum() ).in_units('Msun')
    
    Temp  = halo[('gas','temperature')]
    Mcold=( (halo[('gas','cell_mass')])[Temp<2.5e5].sum() ).in_units('Msun')
    Mhot =( (halo[('gas','cell_mass')])[Temp>2.5e5].sum() ).in_units('Msun')
    #L     = halo.quantities["AngularMomentumVector"]()
    #Lstar = halo.quantities["StarAngularMomentumVector"]()
    #Lambda= halo.quantities["BaryonSpinParameter"]()
    #print L, Lstar, Lambda
    if simulation=="Enzo":
        ct = gal["creation_time"]
        stars=(ct > 0)
        Mstargal = ( (gal["particle_mass"])[stars].sum() ).in_units('Msun')
    if simulation=="Ramses":
        age = gal["particle_age"]
        stars=(age < 0)
        Mstargal = ( (gal["io","particle_mass"])[stars].sum() ).in_units('Msun')
    if simulation=="Arepo" or simulation=="Fire":
        Mstargal = ( gal[('PartType4', 'Masses')].sum() ).in_units('Msun')
    if simulation=="Art":
        Mstargal = ( gal[('stars', 'particle_mass')].sum() ).in_units('Msun')
            
    f = open( "saved_data/"+simulation2+"_basicstats.txt", "w")
    f.write("Mtot = " "%.4g Msun" % Mass)
    f.write("\nMdark = " "%.4g Msun" % Mdark)
    f.write("\nMstar = " "%.4g Msun" % Mstar)
    f.write("\nMgas = " "%.4g Msun" % Mgas)
    f.write("\nMcold = " "%.4g Msun" % Mcold)
    f.write("\nMhot = " "%.4g Msun" % Mhot)
    f.write("\nMstargal = " "%.4g Msun" % Mstargal)
    if simulation=="Enzo":
        f.write("\nMHI = " "%.4g Msun" % MHI)
    f.close()

def density_slices_xyz(halo, pf):
    simulation=get_simulation_type(pf)
    simulation2=simulation+"_z"+str(abs(round(pf.current_redshift, 1)))
    rvir = pf.arr(300.0, 'kpccm')
    center = halo.get_field_parameter("center")

    for ax in 'xyz':
        p=SlicePlot(pf, ax, [('gas','density'), ('gas','temperature')], center=center, width=(4*rvir))
        p.set_zlim(('gas','density'), 1e-30,1e-21)
        p.set_zlim(('gas','temperature'), 1e4, 1e8)
        p.save(simulation2+"_2Rvir")
    for ax in 'xyz':
        p=SlicePlot(pf, ax, [('gas','density'), ('gas','temperature')], center=center, width=(2*rvir))
        p.set_zlim(('gas','density'), 1e-30,1e-21)
        p.set_zlim(('gas','temperature'), 1e4, 1e8)
        p.save(simulation2+"_Rvir")



#NOTE FOR ALL OFFAXIS PLOTS: width HAS to be a pf.arr variable, or it will crash (needs to understand 'unitary' unit)
#quivers only work for Enzo/Ramses so far......
def density_slices_angmom(halo, pf):
    simulation=get_simulation_type(pf)
    simulation2=simulation+"_z"+str(abs(round(pf.current_redshift, 1)))
    rvir=pf.arr(300.0*(1./(1.+pf.current_redshift)), 'kpc') 
    #halo    = pf.h.sphere(center, rvir)
    center = halo.get_field_parameter("center")
    Lgal = halo.get_field_parameter("angular_momentum_vector")
    Lgal = YTArray( [Lgal[0].convert_to_cgs().value, Lgal[1].convert_to_cgs().value, Lgal[2].convert_to_cgs().value], 'cm**2/s')
    Rgal=pf.arr(20.0, 'kpc')
    gal = pf.h.sphere(center, Rgal)
    
    #slices in axis of galaxy angular momentum
    #if simulation=="Enzo":
    #    p=OffAxisSlicePlot(pf, Lgal, ["density", "H_density", "temperature"], center=center, width=(2*rvir))
    #if simulation=="Ramses":
    #    p=OffAxisSlicePlot(pf, Lgal, ["density", "temperature"], center=center, width=(2*rvir), axes_unit='kpc')
    #if simulation=="Arepo" or simulation=="Fire":
    p=OffAxisSlicePlot(pf, Lgal, [('gas','density'), ('gas','temperature')], center=center, width=(2*rvir))
    p.set_zlim(('gas','density'), 1e-30,1e-22)
    p.set_zlim(('gas','temperature'), 1e3, 1e8)
    p.annotate_velocity(factor=32)
    p.save(simulation2+"_L_velocityflow_Rvir")

    p=OffAxisSlicePlot(pf, Lgal, [('gas','density'), ('gas','temperature')], center=center, width=(2*rvir), north_vector=[0,0,1])
    p.set_zlim(('gas','density'), 1e-30,1e-22)
    p.set_zlim(('gas','temperature'), 1e3, 1e8)
    p.annotate_cquiver(('all','vEast'),('all','vNorth'), factor=32)
    p.save(simulation2+"_L_velocityflow_Rvir")
        
    #slices perpendicular to galaxy angular momentum direction, and defines Lgal as "up" on plot
    #if simulation=="Enzo":
    #    p=OffAxisSlicePlot(pf, np.cross(Lgal, [0,1,0]), ["density", "H_density", "temperature"], center=center, width=(2*rvir), north_vector=Lgal/(cm**2/s))
    #if simulation=="Ramses":
    #    p=OffAxisSlicePlot(pf, np.cross(Lgal, [0,1,0]), ["density", "temperature"], center=center, width=(2*rvir), north_vector=Lgal/(cm**2/s))
    #if simulation=="Arepo" or simulation=="Fire":
    p=OffAxisSlicePlot(pf, np.cross(Lgal, [0,1,0]), [('gas','density'), ('gas','temperature')], center=center, width=(2*rvir), north_vector=Lgal)
    p.annotate_velocity(factor=32)
    p.set_zlim(('gas','density'), 1e-30,1e-22)
    p.set_zlim(('gas','temperature'), 1e3, 1e8)
    p.save(simulation2+"_Lperp_velocityflow_Rvir")


def angmom_slices(halo,pf):
    simulation=get_simulation_type(pf)
    simulation2=simulation+"_z"+str(abs(round(pf.current_redshift, 1)))
    rvir=pf.arr(300.0*(1./(1.+pf.current_redshift)), 'kpc') 
    center = halo.get_field_parameter("center")

    #define unit vectors based on angular momentum vector, to use in defining cutting planes 
    bv=halo.get_field_parameter('bulk_velocity')
    Lgal = halo.get_field_parameter("angular_momentum_vector")
    Lgal = YTArray( [Lgal[0].convert_to_cgs().value, Lgal[1].convert_to_cgs().value, Lgal[2].convert_to_cgs().value], 'cm**2/s')
    face_on = Lgal/YTQuantity(np.linalg.norm(Lgal), 'cm**2/s')
    vec = np.array([0,0,1])
    north_vector = np.cross(vec, face_on)
    north_vector = north_vector / np.linalg.norm(north_vector)
    edge_on = np.cross(north_vector,face_on)
    edge_on = edge_on / np.linalg.norm(edge_on) 
    ptype='gas'

    #set density and temperature limits for plots:
    logTmin=3
    logTmax=8
    logrhomin=-30
    logrhomax=-22
    
    #define a region as our data source
    region=pf.h.region(center, [center[0]-1.5*rvir, center[1]-1.5*rvir, center[2]-1.5*rvir], [center[0]+1.5*rvir, center[1]+1.5*rvir, center[2]+1.5*rvir])
    region.set_field_parameter('bulk_velocity', bv)
    region.set_field_parameter('angular_mometnum_vector', Lgal)

    #Define Cutting Plane Velocities
    def _gas_offaxis_velocity_north(field, data):
        v = north_vector[0]*(data[ptype,'velocity_x']-bv[0]) + north_vector[1]*(data[ptype,'velocity_y']-bv[1]) + north_vector[2]*(data[ptype,'velocity_z']-bv[2])
        return v
    def _gas_offaxis_velocity_east(field, data):
        v =      edge_on[0]*(data[ptype,'velocity_x']-bv[0]) +      edge_on[1]*(data[ptype,'velocity_y']-bv[1]) +      edge_on[2]*(data[ptype,'velocity_z']-bv[2])
        return v
    def _gas_offaxis_velocity_normal(field, data):
        v =      face_on[0]*(data[ptype,'velocity_x']-bv[0]) +      face_on[1]*(data[ptype,'velocity_y']-bv[1]) +      face_on[2]*(data[ptype,'velocity_z']-bv[2])
        return v
    pf.add_field(('gas','offaxis_velocity_north'),  function=_gas_offaxis_velocity_north, units='km/s', particle_type=False, take_log=False, force_override=True)
    pf.add_field(('gas','offaxis_velocity_east'),   function=_gas_offaxis_velocity_east,  units='km/s', particle_type=False, take_log=False, force_override=True)
    pf.add_field(('gas','offaxis_velocity_normal'), function=_gas_offaxis_velocity_normal,units='km/s', particle_type=False, take_log=False, force_override=True)

    #look at galaxy face-on
    p=OffAxisSlicePlot(pf, face_on, [('gas','density'), ('gas','temperature')], center=center, width=(2*rvir), north_vector=north_vector, data_source=region)
    #p=OffAxisSlicePlot(pf, face_on, [('gas','density'), ('gas','temperature'), ('gas','offaxis_velocity_normal')], center=center, width=(2*rvir), north_vector=north_vector)
    #p.set_unit(('gas','offaxis_velocity_normal'),'km/s')
    #p.set_zlim(('gas','offaxis_velocity_normal'),-300,300)
    if simulation=="Enzo":
        p.set_zlim(('gas','density'), 1e-30,1e-22)
        p.set_zlim(('gas','temperature'), 1e3, 1e8)
        vx=('gas','offaxis_velocity_east')
        vy=('gas','offaxis_velocity_north')
        p.annotate_cquiver(vx, vy, factor=32)       #The "easy" way only works for grid codes.  Can't figure out how to do it for arepo/FIRE
        p.save(simulation2+"_L_velocityflow_Rvir")

    #MANUALLY read out data from the OffAxisSlicePlot
    T=p.frb[('gas','temperature')]
    den=p.frb[('gas','density')]
    vx=p.frb[('gas','offaxis_velocity_east')]
    vy=p.frb[('gas','offaxis_velocity_north')]
    
    #MANUALLY create velocity quivers
    res = 800
    factor = 32
    lim=float(round(rvir.value))
    nx = res / factor; ny = res / factor
    X,Y = np.meshgrid(np.linspace(-lim,lim,nx,endpoint=True), np.linspace(-lim,lim,ny,endpoint=True))
    pixX = vx[::factor,::factor]
    pixY = vy[::factor,::factor]

    #MANUALLY create Temperature Plot
    fig=plt.figure()
    plt.xlabel('x (kpc)')
    plt.ylabel('y (kpc)')
    im=plt.imshow(np.log10(T), origin='lower', extent=[-lim,lim,-lim,lim], vmin=logTmin, vmax=logTmax)
    plt.quiver(X, Y, pixX, pixY, scale=None, scale_units=None)
    cax =fig.add_axes([0.81, 0.1, 0.03, 0.8])
    cb=fig.colorbar(im, orientation='vertical', cax=cax)
    cb.set_label('log Temperature (K)')
    fig.savefig(simulation2+"_L_velocityflow_Rvir_OffAxisSlice_temperature2.png")

    #MANUALLY create Density Plot
    fig=plt.figure()
    plt.xlabel('x (kpc)')
    plt.ylabel('y (kpc)')
    im=plt.imshow(np.log10(den), origin='lower', extent=[-lim,lim,-lim,lim], vmin=logrhomin, vmax=logrhomax)
    plt.quiver(X, Y, pixX, pixY, scale=None, scale_units=None)
    cax = fig.add_axes([0.81, 0.1, 0.03, 0.8])
    cb=fig.colorbar(im, orientation='vertical', cax=cax)
    cb.set_label('log Density (g/cm^3)')
    fig.savefig(simulation2+"_L_velocityflow_Rvir_OffAxisSlice_density2.png")

    #rotate galaxy upwards, to look at galaxy edge-on --> 'north' is away from you, 'face-on' is up, 'east' is still to right
    p=OffAxisSlicePlot(pf, -1*north_vector, [('gas','density'), ('gas','temperature')], center=center, width=(2*rvir), north_vector=face_on)
    #p=OffAxisSlicePlot(pf, -1*north_vector, [('gas','density'), ('gas','temperature'), ('gas','offaxis_velocity_north')], center=center, width=(2*rvir), north_vector=face_on)
    #p.set_unit(('gas','offaxis_velocity_north'),'km/s')
    #p.set_zlim(('gas','offaxis_velocity_north'),-300,300)
    if simulation=="Enzo":
        p.set_zlim(('gas','density'), 1e-30,1e-22)
        p.set_zlim(('gas','temperature'), 1e3, 1e8)
        vx=('gas','offaxis_velocity_east')
        vy=('gas','offaxis_velocity_normal')
        p.annotate_cquiver(vx, vy, factor=32)        #you can do all this in Enzo with just:  p.annotate_velocity(factor=32)
        p.save(simulation2+"_Lperp_velocityflow_Rvir")

    #MANUALLY read out data from the OffAxisSlicePlot
    T=p.frb[('gas','temperature')]
    den=p.frb[('gas','density')]
    vx=p.frb[('gas','offaxis_velocity_east')]
    vy=p.frb[('gas','offaxis_velocity_normal')]
    
    #MANUALLY create velocity quivers
    res = 800
    factor = 32
    lim=float(round(rvir.value))
    nx = res / factor; ny = res / factor
    X,Y = np.meshgrid(np.linspace(-lim,lim,nx,endpoint=True), np.linspace(-lim,lim,ny,endpoint=True))
    pixX = vx[::factor,::factor]
    pixY = vy[::factor,::factor]

    #MANUALLY create Temperature Plot
    fig=plt.figure()
    plt.xlabel('x (kpc)')
    plt.ylabel('y (kpc)')
    im=plt.imshow(np.log10(T), origin='lower', extent=[-lim,lim,-lim,lim], vmin=logTmin, vmax=logTmax)
    plt.quiver(X, Y, pixX, pixY, scale=None, scale_units=None)
    cax = fig.add_axes([0.81, 0.1, 0.03, 0.8])
    cb=fig.colorbar(im, orientation='vertical', cax=cax)
    cb.set_label('log Temperature (K)')
    fig.savefig(simulation2+"_Lperp_velocityflow_Rvir_OffAxisSlice_temperature2.png")

    #MANUALLY create Density Plot
    fig=plt.figure()
    plt.xlabel('x (kpc)')
    plt.ylabel('y (kpc)')
    im=plt.imshow(np.log10(den), origin='lower', extent=[-lim,lim,-lim,lim], vmin=logrhomin, vmax=logrhomax)
    plt.quiver(X, Y, pixX, pixY, scale=None, scale_units=None)
    cax = fig.add_axes([0.81, 0.1, 0.03, 0.8])
    cb=fig.colorbar(im, orientation='vertical', cax=cax)
    cb.set_label('log Density (g/cm^3)')
    fig.savefig(simulation2+"_Lperp_velocityflow_Rvir_OffAxisSlice_density2.png")
    

#make density and temperature plots aligned with angmom direction
#note that projections cannot overlay velcity "quiver" arrows, the way slices can
def angmom_proj(halo,pf):
    simulation=get_simulation_type(pf)
    simulation2=simulation+"_z"+str(abs(round(pf.current_redshift, 1)))
    rvir=pf.arr(300.0*(1./(1.+pf.current_redshift)), 'kpc') 
    #halo    = pf.h.sphere(center, rvir)
    center = halo.get_field_parameter("center")
    Rgal=pf.arr(20.0, 'kpc')
    gal = pf.h.sphere(center, Rgal)

    #define unit vectors based on angular momentum vector, to use in defining cutting planes (same as in angmom_slices)
    bv=halo.get_field_parameter('bulk_velocity')
    Lgal = halo.get_field_parameter("angular_momentum_vector")
    Lgal = YTArray( [Lgal[0].convert_to_cgs().value, Lgal[1].convert_to_cgs().value, Lgal[2].convert_to_cgs().value], 'cm**2/s')
    face_on = Lgal/YTQuantity(np.linalg.norm(Lgal.in_units('cm**2/s')), "cm**2/s")
    vec = np.array([0,0,1])
    north_vector = np.cross(vec, face_on)
    north_vector = north_vector / np.linalg.norm(north_vector)
    edge_on = np.cross(north_vector,face_on)
    edge_on = edge_on / np.linalg.norm(edge_on) 

    #Define Cutting Plane Velocities
    ptype='gas'
    def _gas_offaxis_velocity_north(field, data):
        v = north_vector[0]*(data[ptype,'velocity_x']-bv[0]) + north_vector[1]*(data[ptype,'velocity_y']-bv[1]) + north_vector[2]*(data[ptype,'velocity_z']-bv[2])
        return v
    def _gas_offaxis_velocity_east(field, data):
        v =      edge_on[0]*(data[ptype,'velocity_x']-bv[0]) +      edge_on[1]*(data[ptype,'velocity_y']-bv[1]) +      edge_on[2]*(data[ptype,'velocity_z']-bv[2])
        return v
    def _gas_offaxis_velocity_normal(field, data):
        v =      face_on[0]*(data[ptype,'velocity_x']-bv[0]) +      face_on[1]*(data[ptype,'velocity_y']-bv[1]) +      face_on[2]*(data[ptype,'velocity_z']-bv[2])
        return v
    pf.add_field(('gas','offaxis_velocity_north'),  function=_gas_offaxis_velocity_north, units='km/s', particle_type=False, take_log=False, force_override=True)
    pf.add_field(('gas','offaxis_velocity_east'),   function=_gas_offaxis_velocity_east,  units='km/s', particle_type=False, take_log=False, force_override=True)
    pf.add_field(('gas','offaxis_velocity_normal'), function=_gas_offaxis_velocity_normal,units='km/s', particle_type=False, take_log=False, force_override=True)

    #look at galaxy face-on
    p=OffAxisProjectionPlot(pf, face_on, [('gas','density'), ('gas','temperature'), ('gas','offaxis_velocity_normal')], weight_field=('gas','density'), center=center, width=(2*rvir), depth=(2*rvir), north_vector=north_vector)
    p.set_zlim(('gas','density'), 1e-30,1e-22)
    p.set_zlim(('gas','temperature'), 1e3, 1e8)
    p.set_unit(('gas','offaxis_velocity_normal'),'km/s')
    p.set_zlim(('gas','offaxis_velocity_normal'),-250,250)
    p.save(simulation2+"_L_proj_Rvir")

    #rotate galaxy upwards, to look at galaxy edge-on --> 'north' is away from you, 'face-on' is up, 'east' is still to right
    p=OffAxisProjectionPlot(pf, -1*north_vector, [('gas','density'), ('gas','temperature'), ('gas','offaxis_velocity_north')], weight_field=('gas','density'), center=center, width=(2*rvir), depth=(2*rvir), north_vector=face_on)
    p.set_zlim(('gas','density'), 1e-30,1e-22)
    p.set_zlim(('gas','temperature'), 1e3, 1e8)
    p.set_unit(('gas','offaxis_velocity_north'),'km/s')
    p.set_zlim(('gas','offaxis_velocity_north'),-250,250)
    p.save(simulation2+"_Lperp_proj_Rvir")


    #look at galaxy face-on but zoom out to Environment
    p=OffAxisProjectionPlot(pf, face_on, [('gas','density'), ('gas','temperature'), ('gas','offaxis_velocity_normal')], weight_field=('gas','density'), center=center, width=(8*rvir), depth=(8*rvir), north_vector=north_vector)
    p.set_zlim(('gas','density'), 1e-30,1e-22)
    p.set_zlim(('gas','temperature'), 1e3, 1e8)
    p.set_unit(('gas','offaxis_velocity_normal'),'km/s')
    p.set_zlim(('gas','offaxis_velocity_normal'),-200,200)
    p.save(simulation2+"_L_proj_Environment")

    #rotate galaxy upwards, to look at galaxy edge-on --> 'north' is away from you, 'face-on' is up, 'east' is still to right, but zoom out to Environment
    p=OffAxisProjectionPlot(pf, -1*north_vector, [('gas','density'), ('gas','temperature'), ('gas','offaxis_velocity_north')], weight_field=('gas','density'), center=center, width=(8*rvir), depth=(8*rvir), north_vector=face_on)
    p.set_zlim(('gas','density'), 1e-30,1e-22)
    p.set_zlim(('gas','temperature'), 1e3, 1e8)
    p.set_unit(('gas','offaxis_velocity_north'),'km/s')
    p.set_zlim(('gas','offaxis_velocity_north'),-200,200)
    p.save(simulation2+"_Lperp_proj_Environment")

def density_projections(halo, pf, plottype):
    simulation=get_simulation_type(pf)
    simulation2=simulation+"_z"+str(abs(round(pf.current_redshift, 1)))
    rvir = pf.arr(300.0, 'kpccm')
    center = halo.get_field_parameter("center")

    if plottype=="Environment" or plottype=="all" or plottype=="All":
        #make a projetion of only the region around the halo (4Rvir radius box)
        region=pf.h.region(center, [center[0]-4*rvir, center[1]-4*rvir, center[2]-4*rvir], [center[0]+4*rvir, center[1]+4*rvir, center[2]+4*rvir])
        for ax in 'xyz':
            p1 = pf.h.proj(('gas','density'), ax, weight_field=('gas','density'), data_source = region).to_pw(fields=('gas','density'), center=center, width=8*rvir)
            p1.set_unit(('gas','density'), "Msun/pc**3")
            p1.set_zlim(('gas','density'), 1e-8, 1)
            p1.set_axes_unit('kpc')
            p1.set_width(8*rvir.in_units('kpc')-1.0*kpc)
            p1.save(simulation2+"_Environment")

    if plottype=="2Rvir" or plottype=="all" or plottype=="All":
        #projected quantities through the entire halo (Rvir radius box)
        region=pf.h.region(center, [center[0]-2*rvir, center[1]-2*rvir, center[2]-2*rvir], [center[0]+2*rvir, center[1]+2*rvir, center[2]+2*rvir])
        for ax in 'x':
            p1 = pf.h.proj([('gas','density'), ('gas','temperature')], ax, weight_field=('gas','density'), data_source = region).to_pw(fields=('gas','density'), center=center, width=4*rvir)
            p1.set_unit(('gas','density'), "Msun/pc**3")
            p1.set_zlim(('gas','density'), 1e-8, 1)
            p1.set_zlim(('gas','temperature'), 1e3, 1e7)
            p1.set_axes_unit('kpc')
            p1.set_width(4*rvir.in_units('kpc')-1.0*kpc)
            p1.save(simulation2+"_2Rvir")

    if plottype=="Halo" or plottype=="all" or plottype=="All":
        #projected quantities through the entire halo (Rvir radius box)
        region=pf.h.region(center, [center[0]-rvir, center[1]-rvir, center[2]-rvir], [center[0]+rvir, center[1]+rvir, center[2]+rvir])
        for ax in 'xyz':
            p1 = pf.h.proj([('gas','density'), ('gas','temperature')], ax, weight_field=('gas','density'), data_source = region).to_pw(fields=('gas','density'), center=center, width=2*rvir)
            p1.set_unit(('gas','density'), "Msun/pc**3")
            p1.set_zlim(('gas','density'), 1e-8, 1)
            p1.set_axes_unit('kpc')
            p1.set_width(2*rvir.in_units('kpc'))
            p1.annotate_scale(coeff=30, unit='kpc')
            p1.set_font_size(30)
            p1.annotate_sphere(center, radius=(rvir.in_units('kpc').value, 'kpc'), circle_args={'color':'black'})
            p1.save(simulation2+"_Halo")


def relative_velocity_proj(halo, pf):
    simulation=get_simulation_type(pf)
    simulation2=simulation+"_z"+str(abs(round(pf.current_redshift, 1)))
    rvir = pf.arr(300.0, 'kpccm')
    #halo    = pf.h.sphere(center, rvir)
    center = halo.get_field_parameter("center")
    bv = halo.get_field_parameter('bulk_velocity')
    Rgal=20.0*kpc
    gal = pf.h.sphere(center, Rgal)

    #do relative velocity maps to virial radius 
    region=pf.h.region(center, [center[0]-rvir, center[1]-rvir, center[2]-rvir], [center[0]+rvir, center[1]+rvir, center[2]+rvir])
    region.set_field_parameter("bulk_velocity", bv)
    vel = "relative_velocity_%s"
    for ax in 'xyz':
        p1 = pf.h.proj(('gas',vel % ax), ax, weight_field=('gas','density'), data_source = region).to_pw(fields=('gas',vel % ax), center=center, width=2*rvir)
        p1.set_unit(('gas',vel % ax),'km/s')
        p1.set_zlim(('gas',vel % ax),-250,250)
        p1.set_axes_unit('kpc')
        p1.set_width(2*rvir.in_units('kpc')-1.0*kpc)
        #p1.annotate_contour(('gas','density'), ncont=5)
        p1.save(simulation2+"_Halo")

    #do relative velocity maps to virial radius, but only for somewhat dense gas
    region=pf.h.region(center, [center[0]-rvir, center[1]-rvir, center[2]-rvir], [center[0]+rvir, center[1]+rvir, center[2]+rvir])
    region.set_field_parameter("bulk_velocity", bv)
    gas_ds = region.cut_region(["(obj[('gas','density')]).in_units('Msun/pc**3') > 7.5e-5"])
    vel = "relative_velocity_%s"
    for ax in 'xyz':
        p1 = pf.h.proj(('gas',vel % ax), ax, weight_field=('gas','density'), data_source = gas_ds).to_pw(fields=('gas',vel % ax), center=center, width=2*rvir)
        p1.set_unit(('gas',vel % ax),'km/s')
        p1.set_zlim(('gas',vel % ax),-250,250)
        p1.set_axes_unit('kpc')
        p1.set_width(2*rvir.in_units('kpc')-1.0*kpc)
        #p1.annotate_contour(('gas','density'), ncont=5)
        p1.save(simulation2+"_Dense_Halo")

    #do relative velocity maps of environment (R<4Rvir), but only for somewhat dense gas
    region=pf.h.region(center, [center[0]-4*rvir, center[1]-4*rvir, center[2]-4*rvir], [center[0]+4*rvir, center[1]+4*rvir, center[2]+4*rvir])
    region.set_field_parameter("bulk_velocity", bv)
    #gas_ds = region.cut_region(["np.logical_and( (obj[('gas','density')]).in_units('Msun/pc**3') > 5e-4, obj[('gas','temperature')] < 2.5e5)"])
    gas_ds = region.cut_region(["(obj[('gas','density')]).in_units('Msun/pc**3') > 7.5e-6"])
    vel = "relative_velocity_%s"
    for ax in 'xyz':
        p1 = pf.h.proj(('gas',vel % ax), ax, weight_field=('gas','density'), data_source = gas_ds).to_pw(fields=('gas',vel % ax), center=center, width=8*rvir)
        p1.set_unit(('gas',vel % ax),'km/s')
        p1.set_zlim(('gas',vel % ax),-200,200)
        p1.set_axes_unit('kpc')
        p1.set_width(8*rvir.in_units('kpc')-2.0*kpc)
        #p1.annotate_contour(('gas','density'), ncont=3)
        p1.save(simulation2+"_Dense_Environment")

    #do relative velocity maps of environment (R<4Rvir), but only for cold and somewhat dense gas
    #region=pf.h.region(center, [center[0]-4*rvir, center[1]-4*rvir, center[2]-4*rvir], [center[0]+4*rvir, center[1]+4*rvir, center[2]+4*rvir])
    #region.set_field_parameter("bulk_velocity", bv)
    #gas_ds = region.cut_region(["np.logical_and( (obj[('gas','density')]).in_units('Msun/pc**3') > 7.5e-6, obj[('gas','temperature')] < 2.5e5)"])
    #vel = "relative_velocity_%s"
    #for ax in 'xyz':
    #    p1 = pf.h.proj(('gas',vel % ax), ax, weight_field=('gas','density'), data_source = gas_ds).to_pw(fields=('gas',vel % ax), center=center, width=8*rvir)
    #    p1.set_unit(('gas',vel % ax),'km/s')
    #    p1.set_zlim(('gas',vel % ax),-200,200)
    #    #p1.annotate_contour(('gas','density'), ncont=5)
    #    p1.save(simulation2+"_Cold_Dense_Environment")


def phase_plots(halo, pf):
    simulation=get_simulation_type(pf)
    simulation2=simulation+"_z"+str(abs(round(pf.current_redshift, 1)))
    rvir=299.0*kpc*(1./(1.+pf.current_redshift))
    #halo    = pf.h.sphere(center, rvir)
    center = halo.get_field_parameter("center")
    Rgal=20.0*kpc
    gal = pf.h.sphere(center, Rgal)

    units  =dict(radial_velocity='km/s',      temperature='K')
    extrema=dict(radial_velocity=(-800,800), temperature=(1e3,1e8))
    logs   =dict(radial_velocity=False,       temperature=True)
    profile=create_profile(halo, [('gas','radial_velocity'), ('gas','temperature')], fields=["cell_mass"], n_bins=[128,128], weight_field=None, fractional=True, logs=logs, units=units, extrema=extrema)
    pp=PhasePlot.from_profile(profile)
    pp.set_zlim('cell_mass', 1e-7,1e-2)
    pp.save(simulation2)
    
    units  =dict(density='Msun/pc**3',     temperature='K')
    extrema=dict(density=(1e-8,1e3), temperature=(1e3,1e8))
    logs   =dict(density=True,          temperature=True)
    profile=create_profile(halo, [('gas','density'), ('gas','temperature')], fields=['cell_mass'], n_bins=[128,128], weight_field=None, fractional=True, logs=logs, units=units, extrema=extrema)
    pp=PhasePlot.from_profile(profile)
    pp.set_zlim('cell_mass', 1e-7,1e-2)
    pp.save(simulation2)

    units  =dict(radius='kpc',    temperature='K')
    extrema=dict(radius=(0,rvir), temperature=(1e3,1e8))
    logs   =dict(radius=False,    temperature=True)
    profile=create_profile(halo, [('gas','radius'), ('gas','temperature')], fields=['cell_mass'], n_bins=[128,128], weight_field=None, fractional=True, logs=logs, units=units, extrema=extrema)
    pp=PhasePlot.from_profile(profile)
    pp.set_zlim('cell_mass', 1e-7,1e-2)
    pp.save(simulation2)
        

def creationtime(halo, pf):
    simulation=get_simulation_type(pf)
    simulation2=simulation+"_z"+str(abs(round(pf.current_redshift, 1)))
    rvir=299.0*kpc*(1./(1.+pf.current_redshift))
    #halo    = pf.h.sphere(center, rvir)
    center = halo.get_field_parameter("center")
    Rgal=20.0*kpc
    gal = pf.h.sphere(center, Rgal)

    #density vs radius profile
    profile=create_profile(gal, 'radius', fields='density', n_bins=64, weight_field=None, fractional=False, logs={'radius':False, 'density':True}, units={'radius':'kpc', 'density':'g/cm**3'}, extrema={'radius':(0,20), 'density':(1e-30, 1e-22)})
    pp=ProfilePlot.from_profiles(profile)
    pp.save(simulation2)

    #creation time vs particle age profile
    profile=create_profile(gal, 'age', fields='particle_mass', n_bins=64, weight_field=None, fractional=False, logs={'age':False, 'particle_mass':True}, units={'age':'Gyr', 'particle_mass':'Msun'}, extrema={'age':(0,14), 'particle_mass':(1e2,1e8)})
    pp=ProfilePlot.from_profiles(profile)
    pp.save(simulation2)

    #plot=ProfilePlot(gal, 'age', 'particle_mass', weight_field=None)
    #plot.set_log('age', False)
    #plot.set_xlim(0,14.0)
    #plot.set_unit('age', 'Gyr')
    #plot.save(simulation2)
    


def SFR(halo, pf):
    simulation=get_simulation_type(pf)
    simulation2=simulation+"_z"+str(abs(round(pf.current_redshift, 1)))
    rvir=299.0*kpc*(1./(1.+pf.current_redshift))
    #halo    = pf.h.sphere(center, rvir)
    center = halo.get_field_parameter("center")
    Rgal=20.0*kpc
    gal = pf.h.sphere(center, Rgal)
    
    if simulation=="Ramses":
        import ramses_get_tform as ram
        ad=pf.all_data()                           #
        ad.set_field_parameter('center',center)    # select all star particles and get their radii from galaxy center, mass, etc.
        stars=ad['particle_age']<0                 #
        radius=(ad['particle_radius'])[stars]      #
        mass = (ad['particle_mass'])[stars]        #
        stars_radius_filter = radius.in_units('kpc') < Rgal  #create filter that goes from all stars --> only stars in "the galaxy"
        ct_stars = ram.load_creation_time(pf)                         # load actual creation time from routine for *all star particles*
        creation_time = ct_stars[stars_radius_filter].in_units('yr')  # and select only those *in radius filter*
        mass = mass[stars_radius_filter].in_units('Msun')
        #print "compare lengths from using sphere object: ", len(creation_time), len( (gal['particle_age'])[gal['particle_age']<0] )
                
    if simulation=="Enzo":
        ct_stars = gal['creation_time'].in_units('yr')
        creation_time = ct_stars[ct_stars > 0]
        mass = (gal['particle_mass'])[ct_stars > 0].in_units('Msun')

    if simulation=="Arepo":
        creation_scale = halo[('PartType4', 'GFM_StellarFormationTime')]         #Gadget-format files give SF time in terms of scale factor
        creation_time = (get_age_from_scale(creation_scale)*Gyr).in_units('yr')  #Convert scale factor to cosmic time using look-up table
    if simulation=="Fire":
        creation_scale = (halo[('PartType4', 'StellarFormationTime')])[stars_radius_filter]             #Gadget-format files give SF time in terms of scale factor
        creation_time = (get_age_from_scale(creation_scale)*Gyr).in_units('yr')  #Convert scale factor to cosmic time using look-up table
        
    time_range = [0, cosmo.age(0).value*1e9] # years
    n_bins = 1000
    hist, bins = np.histogram(creation_time, bins=n_bins, range=time_range,)
    inds = np.digitize(creation_time, bins=bins)
    time = (bins[:-1] + bins[1:])/2
    sfr = np.array([mass[inds == j].sum()/(bins[j+1]-bins[j])
                    for j in range(len(time))])
    sfr[sfr == 0] = np.nan

    plt.plot((time.max()-time)/1e9, sfr)
    plt.xlabel('Lookback Time [Gyr]')
    plt.ylabel('SFR  [M$_\odot$ yr$^{-1}$]')
    plt.ylim(0, 90)
    plt.xlim(cosmo.age(0).value-0.5, 0.0)
    plt.savefig(simulation2+"_SFR.png")


#this function returns numpy arrays for time and sfr, saving them to a file, to be used for plotting purposes.
def SFR2(halo, pf):
    simulation=get_simulation_type(pf)
    simulation2=simulation+"_z"+str(abs(round(pf.current_redshift, 1)))
    rvir=299.0*kpc*(1./(1.+pf.current_redshift))
    #halo    = pf.h.sphere(center, rvir)
    center = halo.get_field_parameter("center")
    #Rgal=20.0*kpc
    Rgal=rvir
    gal = pf.h.sphere(center, Rgal)
    #gal=halo
    
    if simulation=="Ramses":
        import ramses_get_tform as ram
        ad=pf.all_data()                           #
        ad.set_field_parameter('center',center)    # select all star particles and get their radii from galaxy center, mass, etc.
        stars=ad['particle_age']<0                 #
        radius=(ad['particle_radius'])[stars]      #
        mass = (ad['particle_mass'])[stars]        #
        stars_radius_filter = radius.in_units('kpc') < Rgal  #create filter that goes from all stars --> only stars in "the galaxy"
        ct_stars = ram.load_creation_time(pf)                         # load actual creation time from routine for *all star particles*
        creation_time = ct_stars[stars_radius_filter].in_units('yr')  # and select only those *in radius filter*
        mass = mass[stars_radius_filter].in_units('Msun')
        #print "compare lengths from using sphere object: ", len(creation_time), len( (gal['particle_age'])[gal['particle_age']<0] )
              
    if simulation=="Enzo":
        ct_stars = gal['creation_time'].in_units('yr')
        creation_time = ct_stars[ct_stars > 0]
        mass = (gal['particle_mass'])[ct_stars > 0].in_units('Msun')

    if simulation=="Arepo": 
        creation_scale = gal[('PartType4', 'GFM_StellarFormationTime')]          #Gadget-format files give SF time in terms of scale factor
        creation_time = (get_age_from_scale(creation_scale)*Gyr).in_units('yr')  #Convert scale factor to cosmic time using look-up table
        mass = (gal['PartType4','Masses']).in_units('Msun')

    if simulation=="Fire": 
        creation_scale = gal[('PartType4', 'StellarFormationTime')]              #Gadget-format files give SF time in terms of scale factor
        creation_time = (get_age_from_scale(creation_scale)*Gyr).in_units('yr')  #Convert scale factor to cosmic time using look-up table
        mass = (gal['PartType4','Masses']).in_units('Msun')

    if simulation=="Art": 
        creation_time = gal[('stars', 'particle_creation_time')].in_units('yr')  
        #mass = (gal['stars','particle_mass_initial']).in_units('Msun')
        mass = (gal['stars','particle_mass']).in_units('Msun')

    time_range = [0, cosmo.age(0).value*1e9] # years
    n_bins = 1000
    hist, bins = np.histogram(creation_time, bins=n_bins, range=time_range,)
    inds = np.digitize(creation_time, bins=bins)
    time = (bins[:-1] + bins[1:])/2

    sfr = np.array([mass[inds == j].sum()/(bins[j+1]-bins[j])
                    for j in range(len(time))])
    sfr[sfr == 0] = np.nan

    time=np.array(time)
    f = open( simulation+"_SFRvstime.txt", "w")
    f.write("sfr [Msun/yr] \t time [yr]\n")
    for i in range(len(time)):
        f.write("%10.10g \t" % sfr[i])
        f.write("%10.10g \n" % time[i])
    f.close()

    return sfr, time

def SFR_ALL(Enzo_halo, Enzo_pf, Ramses_halo, Ramses_pf, Arepo_halo, Arepo_pf, Fire_halo, Fire_pf):
    if not os.path.isfile("Enzo_SFRvstime.txt"):
        print "creating Enzo SFR vs time table\n"
        E_sfr, E_time = SFR2(Enzo_halo, Enzo_pf)
    else:
        E_sfr, E_time = np.loadtxt('Enzo_SFRvstime.txt', skiprows=1, unpack=True)

    if not os.path.isfile("Ramses_SFRvstime.txt"):
        print "creating Ramses SFR vs time table\n"
        R_sfr, R_time = SFR2(Ramses_halo, Ramses_pf)
    else:
        R_sfr, R_time = np.loadtxt('Ramses_SFRvstime.txt', skiprows=1, unpack=True)

    if not os.path.isfile("Arepo_SFRvstime.txt"):
        print "creating Arepo SFR vs time table\n"
        A_sfr, A_time = SFR2(Arepo_halo, Arepo_pf)
    else:
        A_sfr, A_time = np.loadtxt('Arepo_SFRvstime.txt', skiprows=1, unpack=True)

    if not os.path.isfile("Fire_SFRvstime.txt"):
        print "creating Fire SFR vs time table\n"
        F_sfr, F_time = SFR2(Fire_halo, Fire_pf)
    else:
        F_sfr, F_time = np.loadtxt('Fire_SFRvstime.txt', skiprows=1, unpack=True)

    #E_time, E_sfr = SFR2(Enzo_halo, Enzo_pf)
    #R_time, R_sfr = SFR2(Ramses_halo, Ramses_pf)
    #A_time, A_sfr = SFR2(Arepo_halo, Arepo_pf)
    #F_time, F_sfr = SFR2(Fire_halo, Fire_pf)
    
    fig = plt.figure()
    plt.ylabel('SFR  [M$_\odot$ yr$^{-1}$]')
    plt.ylim(0, 90)
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
    
    E_LBtime = cosmo.age(0).value-E_time/1e9
    ax1.plot(E_LBtime,E_sfr)
    R_LBtime = cosmo.age(0).value-R_time/1e9
    ax1.plot(R_LBtime,R_sfr)
    A_LBtime = cosmo.age(0).value-A_time/1e9
    ax1.plot(A_LBtime,A_sfr)
    F_LBtime = cosmo.age(0).value-F_time/1e9
    ax1.plot(F_LBtime,F_sfr)

    str_z=str(abs(round(Enzo_pf.current_redshift, 1)))
    ax1.legend(("Enzo","Ramses", "Arepo", "Fire"), frameon=False)
    plt.savefig("ALL_SFR_z="+str_z+".png")

#make density and relative velocity plots for various cuts on cold gas
#plottype: choose from 'cold', 'cold_dense', 'cold_lowz', 'cold_inflow'
def cold_gas_plots(halo,pf,plottype):
    simulation=get_simulation_type(pf)
    simulation2=simulation+"_z"+str(abs(round(pf.current_redshift, 1)))
    rvir=299.0*kpc*(1./(1.+pf.current_redshift))
    center = halo.get_field_parameter("center")
    bv = halo.get_field_parameter("bulk_velocity")
    region=pf.h.region(center, [center[0]-rvir, center[1]-rvir, center[2]-rvir], [center[0]+rvir, center[1]+rvir, center[2]+rvir])
    region.set_field_parameter("bulk_velocity", bv)

    #set the data_source for projections based on desired cuts, defined based on each keyword 'plottype'
    if plottype=="cold":
        gas_ds = region.cut_region(["np.logical_and( (obj[('gas','density')]).in_units('Msun/pc**3') > 1e-4, obj[('gas','temperature')] < 2.5e5)"])
    elif plottype=="cold_dense":
        gas_ds = region.cut_region(["np.logical_and( (obj[('gas','density')]).in_units('Msun/pc**3') > 5e-4, obj[('gas','temperature')] < 2.5e5)"])
    elif plottype=="cold_lowz":
        gas_ds = region.cut_region(["np.logical_and( (obj[('gas','density')]).in_units('Msun/pc**3') > 1e-4, obj[('gas','temperature')] < 2.5e5, (obj[ ('gas','metallicity')]).in_units('Zsun') < 0.2)"])
    elif plottype=="cold_inflow":
        gas_ds = region.cut_region(["np.logical_and( (obj[('gas','density')]).in_units('Msun/pc**3') > 1e-4, obj[('gas','temperature')] < 2.5e5, obj[ ('gas','radial_velocity')] < 0.0)"])
    else:
        print "WARNING: plottype not recognized. Please choose from 'cold_gas', 'cold_dense_gas', 'cold_lowz_gas', 'cold_inflow_gas'.\n"
        print "WARNING: reverting to default type: 'cold_gas'."
        plottype='cold_gas'
        gas_ds = region.cut_region(["obj['temperature'] < 2.5e5"])

    #create projections
    savename = simulation2+"_"+plottype+"_gas_halo"
    vel = "relative_velocity_%s"
    for ax in 'xyz':
        p1 = pf.h.proj(('gas','density'), ax, weight_field=('gas','density'), data_source = gas_ds).to_pw(fields=('gas','density'), center=center, width=2*rvir)
        p1.set_unit(('gas','density'), "Msun/pc**3")
        p1.set_zlim(('gas','density'), 1e-4, 1)
        p1.save(savename)
    for ax in 'xyz':
        proj_object = pf.h.proj(('gas',vel % ax), ax, weight_field=('gas','density'), data_source = gas_ds)
        p1=proj_object.to_pw(fields=('gas',vel % ax), center=center, width=(2*rvir))
        p1.set_zlim(('gas',vel % ax),-300,300)
        p1.annotate_contour(('gas','density'), ncont=5)
        p1.save(savename)

def gas_flux(halo, pf):
    simulation=get_simulation_type(pf)
    simulation2=simulation+"_z"+str(abs(round(pf.current_redshift, 1)))
    rvir=299.0*kpc*(1./(1.+pf.current_redshift))
    #halo    = pf.h.sphere(center, rvir)
    center = halo.get_field_parameter("center")
    Rgal=20.0*kpc
    gal = pf.h.sphere(center, Rgal)
    

#take a snapshot and make a pretty 3D plot, uploaded to sketchfab website
#plot is created by making an isocontour of field1 (up to limit of field_limit)
#coloration of the plot is based on the value of field2
#e.g. a plot showing the temperature of dense gas > 5e-27 would use "density", 5e-27, "temperature"
# examples: density=5e-27, coldgas--> temperature=2.5e5
def sketchfab_density(halo, pf, field1, field_limit, field2):
    simulation=get_simulation_type(pf)
    simulation2=simulation+"_z"+str(abs(round(pf.current_redshift, 1)))
    rvir=299.0*kpc*(1./(1.+pf.current_redshift))
    #halo    = pf.h.sphere(center, rvir)
    center = halo.get_field_parameter("center")
    Rgal=20.0*kpc
    gal = pf.h.sphere(center, Rgal)

    apikey = "4234914637944751b2528ca55506d6d2"
    bounds = [(center[i] - rvir, center[i] + rvir) for i in range(3)]

    surf = pf.surface(halo, field1, field_limit)
    upload_id = surf.export_sketchfab(
        title = simulation2+" "+field1 + " contour showing " + field2 + " variations",
        description = "Color map shows changes in " + field2 + " in " + simulation + " simulation at z=" + str(abs(round(pf.current_redshift, 1))),
        color_field = field2,
        color_map = "hot",
        color_log = True,
        bounds = bounds,
        api_key = apikey )

#do all the projection/slice/LOS velocity plots needed for paper 1 (but only the projections/sizes we're actually using in the paper) 
def paper1plots(halo,pf):
    #LOS velocity plots
    simulation=get_simulation_type(pf)
    simulation2=simulation+"_z"+str(abs(round(pf.current_redshift, 1)))
    rvir = pf.arr(300.0, 'kpccm')
    #halo    = pf.h.sphere(center, rvir)
    center = halo.get_field_parameter("center")
    bv = halo.get_field_parameter('bulk_velocity')
    Rgal=20.0*kpc
    gal = pf.h.sphere(center, Rgal)

    do_axes='xyz'
    rounded_redshift = abs(round(pf.current_redshift, 1))
    if rounded_redshift==1.0:
        print "Redshift 1 --> Only creating x-axis plots"
        do_axes='x'
    if rounded_redshift==2.0:
        print "Redshift 2 --> Only creating x-axis plots"
        do_axes='xy'
    if rounded_redshift==3.0:
        print "Redshift 3 --> Only creating z-axis plots"
        do_axes='z'
    if rounded_redshift==4.0 or rounded_redshift==5.0:
        print "Redshift 4 --> Only creating y-axis plots"
        do_axes='y'
        
    #do relative velocity maps to virial radius, but only for somewhat dense gas
    region=pf.h.region(center, [center[0]-rvir, center[1]-rvir, center[2]-rvir], [center[0]+rvir, center[1]+rvir, center[2]+rvir])
    region.set_field_parameter("bulk_velocity", bv)
    gas_ds = region.cut_region(["np.logical_and( (obj[('gas','density')]).in_units('Msun/pc**3') > 7.5e-5, obj[('gas','temperature')] < 2.5e5)"])
    #gas_ds = region.cut_region(["(obj[('gas','density')]).in_units('Msun/pc**3') > 7.5e-5"])
    vel = "relative_velocity_%s"
    for ax in do_axes:
        p1 = pf.h.proj(('gas',vel % ax), ax, weight_field=('gas','density'), data_source = gas_ds).to_pw(fields=('gas',vel % ax), center=center, width=2*rvir)
        p1.set_unit(('gas',vel % ax),'km/s')
        p1.set_zlim(('gas',vel % ax),-250,250)
        p1.set_axes_unit('kpc')
        p1.set_width(2*rvir.in_units('kpc'))
        p1.annotate_scale(coeff=25, unit='kpc', size_bar_args={'color':'black'}, corner='lower_left')
        p1.set_font_size(30)
        p1.annotate_sphere(center, radius=(rvir.in_units('kpc').value, 'kpc'), circle_args={'color':'black'})
        #p1.annotate_contour(('gas','density'), ncont=5)
        p1.save(simulation2+"_Dense_Halo")

    #do relative velocity maps of environment (R<4Rvir), but only for somewhat dense gas
    region=pf.h.region(center, [center[0]-4*rvir, center[1]-4*rvir, center[2]-4*rvir], [center[0]+4*rvir, center[1]+4*rvir, center[2]+4*rvir])
    region.set_field_parameter("bulk_velocity", bv)
    gas_ds = region.cut_region(["np.logical_and( (obj[('gas','density')]).in_units('Msun/pc**3') > 7.5e-6, obj[('gas','temperature')] < 2.5e5)"])
    #gas_ds = region.cut_region(["(obj[('gas','density')]).in_units('Msun/pc**3') > 7.5e-6"])
    vel = "relative_velocity_%s"
    for ax in do_axes:
        p1 = pf.h.proj(('gas',vel % ax), ax, weight_field=('gas','density'), data_source = gas_ds).to_pw(fields=('gas',vel % ax), center=center, width=8*rvir)
        p1.set_unit(('gas',vel % ax),'km/s')
        p1.set_zlim(('gas',vel % ax),-200,200)
        p1.set_axes_unit('kpc')
        p1.set_width(8*rvir.in_units('kpc'))
        p1.annotate_scale(coeff=100, unit='kpc', size_bar_args={'color':'black'}, corner='lower_left')
        p1.set_font_size(30)
        p1.annotate_sphere(center, radius=(rvir.in_units('kpc').value, 'kpc'), circle_args={'color':'black'})
        #p1.annotate_contour(('gas','density'), ncont=3)
        p1.save(simulation2+"_Dense_Environment")

    if True:
    #if rounded_redshift==3:
        for ax in do_axes:
            #make a projetion of only the region around the halo (4Rvir radius box)
            region=pf.h.region(center, [center[0]-4*rvir, center[1]-4*rvir, center[2]-4*rvir], [center[0]+4*rvir, center[1]+4*rvir, center[2]+4*rvir])
            p1 = pf.h.proj(('gas','density'), ax, weight_field=('gas','density'), data_source = region).to_pw(fields=('gas','density'), center=center, width=8*rvir)
            p1.set_unit(('gas','density'), "Msun/pc**3")
            p1.set_zlim(('gas','density'), 1e-8, 1)
            p1.set_axes_unit('kpc')
            p1.set_width(8*rvir.in_units('kpc'))
            p1.annotate_scale(coeff=100, unit='kpc', size_bar_args={'color':'w'}, corner='lower_left')
            p1.set_font_size(30)
            p1.annotate_sphere(center, radius=(rvir.in_units('kpc').value, 'kpc'), circle_args={'color':'w'})
            p1.save(simulation2+"_Environment")

            #make a projetion of only the halo (Rvir radius box)
            region=pf.h.region(center, [center[0]-rvir, center[1]-rvir, center[2]-rvir], [center[0]+rvir, center[1]+rvir, center[2]+rvir])
            p1 = pf.h.proj(('gas','density'), ax, weight_field=('gas','density'), data_source = region).to_pw(fields=('gas','density'), center=center, width=2*rvir)
            p1.set_unit(('gas','density'), "Msun/pc**3")
            p1.set_zlim(('gas','density'), 1e-8, 1)
            p1.set_axes_unit('kpc')
            p1.set_width(2*rvir.in_units('kpc'))
            p1.annotate_scale(coeff=25, unit='kpc', size_bar_args={'color':'w'}, corner='lower_left')
            p1.set_font_size(30)
            p1.annotate_sphere(center, radius=(rvir.in_units('kpc').value, 'kpc'), circle_args={'color':'w'})
            p1.save(simulation2+"_Halo")

        for ax in 'x':
            #projected quantities through 2Rvir radius box
            region=pf.h.region(center, [center[0]-2*rvir, center[1]-2*rvir, center[2]-2*rvir], [center[0]+2*rvir, center[1]+2*rvir, center[2]+2*rvir])
            p1 = pf.h.proj([('gas','density'), ('gas','temperature')], ax, weight_field=('gas','density'), data_source = region).to_pw(fields=('gas','density'), center=center, width=4*rvir)
            p1.set_unit(('gas','density'), "Msun/pc**3")
            p1.set_zlim(('gas','density'), 1e-8, 1)
            p1.set_zlim(('gas','temperature'), 1e3, 1e7)
            p1.set_axes_unit('kpc')
            p1.set_width(4*rvir.in_units('kpc'))
            p1.annotate_scale(coeff=50, unit='kpc', size_bar_args={'color':'w'}, corner='lower_left')  #this crashes yt here for some reason...
            p1.set_font_size(30)
            p1.annotate_sphere(center, radius=(rvir.in_units('kpc').value, 'kpc'), circle_args={'color':'w'})
            p1.save(simulation2+"_2Rvir")
        

def firemovie():
    from load_simulation import load_saved_halo_parameters
    import matplotlib.pyplot as plt
    from scipy.interpolate import interp1d

    simulation="Fire"
    Fz=[1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.20, 2.50, 2.60, 2.80, 3.00, 3.40, 3.82, 3.90, 4.00, 4.50, 5.00, 10.0]
    F_time=cosmo.lookback_time(Fz)

    F_center = np.zeros(F_time.shape[0]*3).reshape(F_time.shape[0], 3)
    for i in range(F_time.shape[0]):
        redshift = Fz[i]
        center,bv,Lgal = load_saved_halo_parameters(simulation,redshift)
        F_center[i,:]=center

    Fz_long, trash= np.loadtxt('Fire_redshift.txt', dtype={'names':('redshift','fname'), 'formats':(np.float,'|S8')}, unpack=True)
    Fz_long = Fz_long[Fz_long<=10.0]
    Fz_long = Fz_long[Fz_long>=1.0]
    F_time_long = cosmo.lookback_time(Fz_long)

    fx = interp1d(F_time, F_center[:,0])
    fy = interp1d(F_time, F_center[:,1])
    fz = interp1d(F_time, F_center[:,2])

    x_long = fx(F_time_long)
    y_long = fy(F_time_long)
    z_long = fz(F_time_long)
    
    #plot the interpolations to make sure they work
    z_tick_values = [6, 5, 4, 3, 2, 1.5, 1.25, 1.0 ]  #max and min of this are limits of the plot.  
    z_tick_value_labels = [' ', '5', '4', '3', '2', '1.5', '1.25', '1.0' ]  
    z_tick_locations = cosmo.lookback_time(z_tick_values).value
    LBT_max = cosmo.lookback_time(np.max(z_tick_values)).value
    LBT_min = cosmo.lookback_time(np.min(z_tick_values)).value

    fig = plt.figure()
    plt.ylabel('x')
    plt.ylim(9000, 14000)
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twiny()
    ax1.set_xlabel('Lookback Time [Gyr]')
    ax1.set_xlim(LBT_max, LBT_min)
    ax2.set_xlabel('Redshift')
    ax2.set_xlim(LBT_max, LBT_min)
    ax2.set_xticks(z_tick_locations)
    ax2.set_xticklabels(z_tick_value_labels)
    ax1.plot(F_time,F_center[:,0], 'o')
    ax1.plot(F_time,F_center[:,1], 'x')
    ax1.plot(F_time,F_center[:,2], '.')
    ax1.plot(F_time_long, x_long)
    ax1.plot(F_time_long, y_long)
    ax1.plot(F_time_long, z_long)
    plt.close()
        

    istart = 20 # start at arbirtary snapshot (must be >= 20 -- corresponding to z=10) 
    istop = 290 # stop at snap290 --> z=1
    for i in range(istop-istart+1):
        j=i+istart
        snap="/home/kstewart/RyanMW/gadget/snapdir_%03d/snapshot_%03d.0.hdf5"%(j,j)
        print "loading snapshot...%03d"%j
        #print [x_long[i-i0], y_long[i-i0], z_long[i-i0]]
        pf=load(snap)
        center=pf.arr(  [fx(cosmo.lookback_time(pf.current_redshift)), fy(cosmo.lookback_time(pf.current_redshift)), fz(cosmo.lookback_time(pf.current_redshift))], 'code_length')
        print "Center = ", center
        rvir = pf.arr(300.0, 'kpccm')
        halo=pf.h.sphere(center, rvir)
        simulation2="Fire_snap%03d_z"%j+str(abs(round(pf.current_redshift, 4)))
        ztext="z = %1.2f"%pf.current_redshift

        for ax in 'x':
            #make a projetion of only the region around the halo (4Rvir radius box)
            region=pf.h.region(center, [center[0]-4*rvir, center[1]-4*rvir, center[2]-4*rvir], [center[0]+4*rvir, center[1]+4*rvir, center[2]+4*rvir])
            p1 = pf.h.proj(('gas','density'), ax, weight_field=('gas','density'), data_source = region).to_pw(fields=('gas','density'), center=center, width=8*rvir)
            p1.set_unit(('gas','density'), "Msun/pc**3")
            p1.set_zlim(('gas','density'), 1e-8, 1)
            p1.set_axes_unit('kpc')
            p1.set_width(8*rvir.in_units('kpc')-2.0*kpc)
            p1.annotate_text((0.05, 0.05), ztext, text_args={'color':'white', 'size':'x-large'}, coord_system='axis')
            p1.save(simulation2+"_Environment")

            #make a projetion of only the halo (Rvir radius box)
            region=pf.h.region(center, [center[0]-rvir, center[1]-rvir, center[2]-rvir], [center[0]+rvir, center[1]+rvir, center[2]+rvir])
            p1 = pf.h.proj(('gas','density'), ax, weight_field=('gas','density'), data_source = region).to_pw(fields=('gas','density'), center=center, width=2*rvir)
            p1.set_unit(('gas','density'), "Msun/pc**3")
            p1.set_zlim(('gas','density'), 1e-8, 1)
            p1.set_axes_unit('kpc')
            p1.set_width(2*rvir.in_units('kpc')-2.0*kpc)
            p1.annotate_text((0.05, 0.05), ztext, text_args={'color':'white', 'size':'x-large'}, coord_system='axis')
            p1.save(simulation2+"_Halo")


        
