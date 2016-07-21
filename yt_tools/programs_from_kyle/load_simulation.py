from yt.mods import *
from yt.utilities.math_utils import get_sph_r_component, get_sph_theta, get_sph_phi, get_sph_theta_component, get_sph_phi_component
from yt.units import kpc, Msun, km, s, cm, g, Gyr, yr
import os as os
from get_snapshot_name import get_snapshot_name as getsnap

#loads the appropriate simulation based on simulation name and redshift (both as strings)
#value, center = pf.find_max('density') determines halo center (long computation so this is pre-calculated, not run every time)
#also pre-load bulk velocity, angular momentum, center into "halo" sphere object with halo.quantities.____() commands
#syntax:   halo, pf = load_simulation(simulation, redshift).  halo is a sphere object focused on main galaxy halo
#syntax of input variable "simulation" is NOT case-sensitive
def load_simulation(simulation, redshift):
    pf=False
    simulation = simulation.lower()[0].upper()+simulation[1:].lower()   #make input simulation-name case-insensitive

    #-------------------------------------------------------- ENZO -----------------------------------------------------------------#
    if simulation=="Enzo":
        if 'linux' in os.sys.platform:
            pf=load(getsnap('enzo',redshift))
        else:
            if redshift==3:
                pf=load("/Users/kstewart/Research/RyanMW/Enzo/RD0018/redshift0018")
            if redshift==2:
                pf=load("/Users/kstewart/Research/RyanMW/Enzo/RD00xx/redshift00xx")
            if redshift==1:
                pf=load("/Users/kstewart/Research/RyanMW/Enzo/RD0038/redshift0038")
            if redshift==0:
                pf=load("/Users/kstewart/Research/RyanMW/Enzo/RD0058/redshift0058")
    #-------------------------------------------------------- ENZO -----------------------------------------------------------------#

    #-------------------------------------------------------- RAMSES -----------------------------------------------------------------#
    if simulation=="Ramses":
        if 'linux' in os.sys.platform:
            pf=load(getsnap('ramses',redshift))
        else:
            if redshift==8:
                pf=load("/Users/kstewart/Research/RyanMW/Ramses/output_00018/info_00018.txt")
            if redshift==4:
                pf=load("/Users/kstewart/Research/RyanMW/Ramses/output_00039/info_00039.txt")
            if redshift==3:
                pf=load("/Users/kstewart/Research/RyanMW/Ramses/output_00060/info_00060.txt")
            if redshift==2:
                pf=load("/Users/kstewart/Research/RyanMW/Ramses/output_00067/info_00067.txt")
            if redshift==1:
                pf=load("/Users/kstewart/Research/RyanMW/Ramses/output_00102/info_00102.txt")
            if redshift=="0.5":
                pf=load("/Users/kstewart/Research/RyanMW/Ramses/output_00131/info_00131.txt")
            if redshift==0:
                pf=load("/Users/kstewart/Research/RyanMW/Ramses/output_00210/info_00210.txt")
    #-------------------------------------------------------- RAMSES -----------------------------------------------------------------#

    #-------------------------------------------------------- ART -----------------------------------------------------------------#
    if simulation=="Art":
        if 'linux' in os.sys.platform:
            pf=load(getsnap('art',redshift))
        else:
            if redshift==4:
                pf=load("")
            if redshift==3:
                pf=load("")
            if redshift==2:
                pf=load("")
    #-------------------------------------------------------- ART -----------------------------------------------------------------#

    #-------------------------------------------------------- AREPO -----------------------------------------------------------------#
    if simulation=="Arepo":
        #center location determined by command: value, center = pf.h.find_max('all_density')
        #units_override={'velocity_unit':(1.0,'kmcm/s')}
        if 'linux' in os.sys.platform:
            pf=load(getsnap('arepo',redshift))
        else:
            if redshift==4:
                pf=load("/Users/kstewart/Research/RyanMW/Arepo/snapdir_060/snap_060.0.hdf5")
            if redshift==3:
                pf=load("/Users/kstewart/Research/RyanMW/Arepo/snapdir_060/snap_060.0.hdf5")
            if redshift==2:
                pf=load("/Users/kstewart/Research/RyanMW/Arepo/snapdir_060/snap_070.0.hdf5")
    #-------------------------------------------------------- AREPO -----------------------------------------------------------------#

    #-------------------------------------------------------- FIRE -----------------------------------------------------------------#
    if simulation=="Fire":
        #center location determined by command: value, center = pf.h.find_max('all_density')
        if 'linux' in os.sys.platform:
            pf=load(getsnap('gadget',redshift))  #n_ref=8*64 smooth out (lower resolution) by factor of 8
        else:
            if redshift==4:
                pf=load("/Users/kstewart/Research/RyanMW/Fire/snapdir_090/snapshot_090.0.hdf5") #smooth out (lower resolution) by factor of 8
            if redshift==3:
                pf=load("/Users/kstewart/Research/RyanMW/Fire/snapdir_140/snapshot_140.0.hdf5") #smooth out (lower resolution) by factor of 8
            if redshift==2:
                pf=load("/Users/kstewart/Research/RyanMW/Fire/snapdir_190/snapshot_190.0.hdf5") #smooth out (lower resolution) by factor of 8
    #-------------------------------------------------------- FIRE -----------------------------------------------------------------#

    
    #-------------------------------------------------------- DERIVED FIELDS -----------------------------------------------------------------#
    #Grid based (Enzo/Ramses) specific derived fields
    if simulation=="Enzo" or simulation=="Ramses":
        pf.add_field( ('gas', 'radius'), function=_grid_alias_gas_radius, particle_type=False, units='kpc')

    #Art specific derived fields 
    if simulation=="Art":
        #pf.add_field( ('gas', 'temperature'), function=_art_temperature, particle_type=False, units='K', force_override=True)
        pf.add_field( ('gas', 'radius'), function=_grid_alias_gas_radius, particle_type=False, units='kpc')
        #pf.add_field( ('gas', 'velocity_x'), function=_art_gas_velocity_x, particle_type=False, units='km/s')
        #pf.add_field( ('gas', 'velocity_y'), function=_art_gas_velocity_y, particle_type=False, units='km/s')
        #pf.add_field( ('gas', 'velocity_z'), function=_art_gas_velocity_z, particle_type=False, units='km/s')
        pf.add_field( ('all', 'velocity_x'), function=_sph_alias_velocity_x, particle_type=True, units='km/s')
        pf.add_field( ('all', 'velocity_y'), function=_sph_alias_velocity_y, particle_type=True, units='km/s')
        pf.add_field( ('all', 'velocity_z'), function=_sph_alias_velocity_z, particle_type=True, units='km/s')

    #Arepo / Fire (basically, Gadget-format) specific derived fields
    if simulation=="Arepo" or simulation=="Fire":
        #create aliases to map sph quantities to use the same names as grid codes use.
        pf.add_field( ('all', 'velocity_x'), function=_sph_alias_velocity_x, particle_type=True, units='km/s')
        pf.add_field( ('all', 'velocity_y'), function=_sph_alias_velocity_y, particle_type=True, units='km/s')
        pf.add_field( ('all', 'velocity_z'), function=_sph_alias_velocity_z, particle_type=True, units='km/s')
        #pf.add_field( ('gas', 'velocity_x'), function=_sph_alias_gas_velocity_x, particle_type=False, units='km/s')
        #pf.add_field( ('gas', 'velocity_y'), function=_sph_alias_gas_velocity_y, particle_type=False, units='km/s')
        #pf.add_field( ('gas', 'velocity_z'), function=_sph_alias_gas_velocity_z, particle_type=False, units='km/s')
        pf.add_field( ('all', 'specific_angular_momentum_x'), function=_sph_alias_specific_angular_momentum_x, particle_type=True, units='cm**2/s')
        pf.add_field( ('all', 'specific_angular_momentum_y'), function=_sph_alias_specific_angular_momentum_x, particle_type=True, units='cm**2/s')
        pf.add_field( ('all', 'specific_angular_momentum_z'), function=_sph_alias_specific_angular_momentum_x, particle_type=True, units='cm**2/s')
        #pf.add_field( ('gas', 'specific_angular_momentum_x'), function=_sph_alias_gas_specific_angular_momentum_x, particle_type=False, units='cm**2/s')
        #pf.add_field( ('gas', 'specific_angular_momentum_y'), function=_sph_alias_gas_specific_angular_momentum_x, particle_type=False, units='cm**2/s')
        #pf.add_field( ('gas', 'specific_angular_momentum_z'), function=_sph_alias_gas_specific_angular_momentum_x, particle_type=False, units='cm**2/s')
        pf.add_field( ('all', 'radius'), function=_sph_alias_radius, particle_type=True, units='kpc')
        pf.add_field( ('gas', 'radius'), function=_grid_alias_gas_radius, particle_type=False, units='kpc')
        #pf.add_field(('gas', "radial_velocity"), function=_gas_radial_velocity, particle_type=False, units="km/s", force_override=True)
        pf.add_field(('all', "radial_velocity"), function=_particle_radial_velocity, particle_type=True, units="km/s", force_override=True)
        #Need to define 'unitary' units for these simulations as well
        #DW = pf.arr(pf.domain_right_edge - pf.domain_left_edge, "code_length")
        #pf.unit_registry.add("unitary", float(DW.max() * DW.units.cgs_value), DW.units.dimensions)

    #New derived fields for ALL simulation types 
    pf.add_field( ('gas', 'relative_velocity_x'), units='km/s', function=_gas_relative_velocity_x, particle_type=False, take_log=False)
    pf.add_field( ('gas', 'relative_velocity_y'), units='km/s', function=_gas_relative_velocity_y, particle_type=False, take_log=False)
    pf.add_field( ('gas', 'relative_velocity_z'), units='km/s', function=_gas_relative_velocity_z, particle_type=False, take_log=False)
    #pf.add_field(('gas','velocity_north'),  function=_gas_offaxis_velocity_north, units='km/s', particle_type=False, take_log=False)
    #pf.add_field(('gas','velocity_east'),   function=_gas_offaxis_velocity_east,  units='km/s', particle_type=False, take_log=False)
    #pf.add_field(('gas','velocity_normal'), function=_gas_offaxis_velocity_normal,units='km/s', particle_type=False, take_log=False)

    #-------------------------------------------------------- DERIVED FIELDS -----------------------------------------------------------------#

    #for ANY SIMULATION, load the halo center, bulk velocity, and galaxy angular momentum from pre-saved values, then define the halo
    center,bv,Lgal = load_saved_halo_parameters(simulation,redshift)
    center=pf.arr(center, 'code_length')
    rvir = pf.arr(300.0, 'kpccm')
    halo    = pf.h.sphere(center, rvir)
    halo.set_field_parameter("bulk_velocity", bv)
    halo.set_field_parameter("angular_momentum_vector", Lgal)
    halo.set_field_parameter("normal", YTArray([0,0,1]))
    return halo, pf


#returns saved halo paraemters based on redshift and simulation name
#format of simulation name is NOT case sensitive.  Options are -- Enzo, Ramses, Art, Fire, Arepo
#center = center of main halo -- determined by command: value, center = pf.h.find_max('all_density')
#bv = bulk velocity of halo
#Lgal = angular momentum of the central galaxy (R<20kpc) of the halo
#returns a value of 0 (and prints error message) for any parameter that is not yet defined.
def load_saved_halo_parameters(simulation,redshift):
    center = bv = Lgal = YTArray([0,0,0])
    simulation = simulation.lower()[0].upper()+simulation[1:].lower()   #make input simulation-name case-insensitive
    
    #-------------------------------------------------------- ENZO -----------------------------------------------------------------#
    if simulation=="Enzo":
        if redshift==3:
            #center=pf.arr([3.671374166e-01, 4.492826843e-01, 5.193484562e-01], 'code_length') # 0-1 units
            center = YTArray( [ 0.3670425415039062, 0.4494857788085938, 0.5195541381835888]) # 'code_length' # 0-1 units
            bv     = YTArray( [ 7000411.32979963,  5880285.02591697, -7374647.88594685], 'cm/s')
            Lgal   = YTArray( [ 9.55343663206e+28, 4.19386247426e+28, 1.78052400375e+28 ], 'cm**2/s' )
            rvir=7.915565054e-03 #unitary -> 70.67 kpc = 283 ckpc
            Mvir=3.079760837e+11
        #if redshift==2:
            #center=pf.arr([0.4032863361, 0.4717615279, 0.4613209662])#, 'code_length') # 0-1 units
            #rvir=8.796488424e-03 #unitary -> 314.16 kpc
            #Mvir=1.472159713e+12
        if redshift==1:
            center=YTArray([3.858510505e-01, 4.606776594e-01, 4.917943380e-01])#, 'code_length') # 0-1 units
            rvir=9.057392092e-03 #unitary -> 161.34 kpc = 323 ckpc
            Mvir=8.249904732e+11
        if redshift==0:
            center=YTArray([4.032168108e-01, 4.717957861e-01, 4.612480366e-01])#, 'code_length') # 0-1 units
            bv     = YTArray( [ 7182080.88471607,  4306191.46708033, -8124103.42890124], 'cm/s')
            Lgal   = YTArray( [ 5.59964164918e+28, -2.72588675515e+29, -1.59971211606e+29], 'cm**2/s' )
            rvir=8.796488424e-03 #unitary -> 314.16 kpc = 314 ckpc
            Mvir=1.472159713e+12
    #-------------------------------------------------------- ENZO -----------------------------------------------------------------#

    #-------------------------------------------------------- RAMSES -----------------------------------------------------------------#
    if simulation=="Ramses":
        if redshift==8:
            center=YTArray([0.35783386,  0.43812561,  0.52461243])#, 'code_length') # 0-1 units
        if redshift==4:
            center= YTArray([ 0.3652114868164062, 0.4456253051757812, 0.5204849243164062])#, 'code_length') # 0-1 units
            bv   =  YTArray([ 7541049.23334357,  4810140.57765677, -2958693.19227095], 'cm/s' )
            Lgal =  YTArray([ 2.60769695502e+28, -3.3472615748e+28, 5.59879117032e+27], 'cm**2/s')
        if redshift==3:
            center=YTArray([ 3.684333980e-01, 4.503621044e-01, 5.163866495e-01])#, 'code_length') # 0-1 units
            bv   = YTArray([ 5146727.94147222,  4304587.59521447, -6569302.99548824], 'cm/s' )
            Lgal = YTArray([  3.55172677e+28,   5.43672100e+28,  -4.07964768e+28], 'cm**2/s')
            rvir = 9.058849889e-03 #unitary -> 80.74 kpc = 323 ckpc
            Mvir=4.507163417e+11
        if redshift==2:
            #center=[3.767523944e-01, 4.542925497e-01, 5.051947237e-01], 'code_length') # 0-1 units. This is offset from galaxy due to merger
            center=YTArray([0.37549591,  0.45398712,  0.50614166])#, 'code_length')
            bv =   YTArray( [ 6684080.15489569,  4826907.63566979, -9797914.44517351], 'cm/s')
            Lgal = YTArray( [4.46832855915e+27, 1.3701426861e+28, -7.46806248209e+28], 'cm**2/s')
            rvir=1.254609909e-02 #unitary -> 149 kpc = 448 ckpc
            Mvir=6.630044935e+11
        if redshift==1.5:
            center=YTArray([3.804927702e-01,	4.567554563e-01,	4.966276947e-01])#, 'code_length') # 0-1 units
            rvir=9.645748667e-03 #unitary -> 137 kpc = 334 ckpc
            Mvir=8.249583087e+11
        if redshift==1:
            center=YTArray([3.853373685e-01, 4.601577714e-01, 4.867490185e-01])#, 'code_length') # 0-1 units
            rvir=8.189311630e-03 #unitary -> 146 kpc = 292 ckpc
            Mvir=9.601470413e+11
        if redshift=="0.5":
            center=YTArray([3.853373685e-01, 4.601577714e-01, 4.867490185e-01])#, 'code_length') # 0-1 units
        if redshift==0:
            center=YTArray([4.041181896e-01, 4.701136575e-01, 4.570915379e-01])#, 'code_length') # 0-1 units
            bv   = YTArray([ 6716207.18482715,  3058015.80276357, -9069901.28064168], 'cm/s' )
            Lgal = YTArray([ 8.66592104934e+28, -8.36982710857e+28, -1.0708525743e+29], 'cm**2/s')
            rvir=9.067457497e-03 #unitary -> 323 kpc = 323 ckpc
            Mvir=1.702014832e+12
    #-------------------------------------------------------- RAMSES -----------------------------------------------------------------#

    #-------------------------------------------------------- ART -----------------------------------------------------------------#
    if simulation=="Art":
        if redshift==4:
            center= YTArray([ 0.4473209381103516, 0.3720149993896484, 0.3946857452392578])#, 'code_length') # 0-1 units
            #bv   = YTArray([ 7541049.23334357,  4810140.57765677, -2958693.19227095], 'cm/s' )
            #Lgal = YTArray([ 2.60769695502e+28, -3.3472615748e+28, 5.59879117032e+27], 'cm**2/s')
        if redshift==3:
            center= YTArray([ 0.4919796 ,  0.41488075,  0.43264961])#, 'code_length') # 0-1 units
            #bv   = YTArray([ 7541049.23334357,  4810140.57765677, -2958693.19227095], 'cm/s' )
            #Lgal = YTArray([ 2.60769695502e+28, -3.3472615748e+28, 5.59879117032e+27], 'cm**2/s')
        if redshift==2:
            center= YTArray([ 0.4017314910888672, 0.5545940399169922, 0.5247402191162109])#, 'code_length') # 0-1 units
            #bv   = YTArray([ 7541049.23334357,  4810140.57765677, -2958693.19227095], 'cm/s' )
            #Lgal = YTArray([ 2.60769695502e+28, -3.3472615748e+28, 5.59879117032e+27], 'cm**2/s')
    #-------------------------------------------------------- ART -----------------------------------------------------------------#

    #-------------------------------------------------------- AREPO -----------------------------------------------------------------#
    if simulation=="Arepo":
        if redshift==4:
            #center=pf.arr([0.36784959,  0.44960427,  0.51697516], 'code_length') # 0-1 units
            center = YTArray( [ 9040.3437614440917969, 11053.3595085144042969, 13032.4244499206542969])#, 'code_length') # 1 kpccm/h units for Arepo
            bv     = YTArray( [ 118.21270788,  108.31640739,  -65.74191243], 'km/s' )  #bulk velocity of galaxy
            Lgal   = YTArray( [ -7.24897984735e+28, 1.09763521083e+27, -2.27694563177e+28], 'cm**2/s') #ang mom vector for galaxy
          
        if redshift==3:
            #center=pf.arr([0.36784959,  0.44960427,  0.51697516], 'code_length') # 0-1 units
            center = YTArray( [ 9196.23970985,  11240.10682106,  12924.37911034])#, 'code_length') # 1 kpccm/h units for Arepo
            bv     = YTArray( [ 67.35413895,  139.12740292, -151.01760116], 'km/s' )  #bulk velocity of galaxy
            Lgal   = YTArray( [ 1.23397273703e+29, 2.83824709682e+28, 4.18442793656e+27], 'cm**2/s') #ang mom vector for galaxy
          
        if redshift==2:
            #center=pf.arr([0.36784959,  0.44960427,  0.51697516], 'code_length') # 0-1 units
            center = YTArray( [ 9404.0215015411376953, 11365.7414913177490234, 12640.5656337738037109])#, 'code_length') # 1 kpccm/h units for Arepo
            bv     = YTArray( [ 181.36570404,   42.34693956, -144.31245317], 'km/s' )  #bulk velocity of galaxy
            Lgal   = YTArray( [ 9.10986183823e+28, 6.28365846327e+28, -6.73906653098e+28], 'cm**2/s') #ang mom vector for galaxy
    #-------------------------------------------------------- AREPO -----------------------------------------------------------------#

    #-------------------------------------------------------- FIRE -----------------------------------------------------------------#
    if simulation=="Fire":
        #center location determined by command: value, center = pf.h.find_max('all_density')
        if redshift==4:
            #center = pf.arr([0.36773559,  0.44975769,  0.51700393], 'code_length') # 0-1 units
            center = YTArray( [  9101.0749340057373047, 11110.8362674713134766, 12984.8897457122802734])#, 'code_length') # 1 kpccm/h units for Arepo
            ##bv   = YTArray( [ 7.50038669,  132.1431487 , -221.68740748], 'km/s' )  #bulk velocity of galaxy (everything)
            bv     = YTArray( [  122.28942968,  102.63906415,  -71.3656294], 'km/s' )  #bulk velocity of galaxy (only gas)            
            Lgal   = YTArray( [-6.12388547094e+27, -4.10178305753e+28, 8.35080636406e+27], 'cm**2/s') #ang mom vector for galaxy

        if redshift==3:
            #center = pf.arr([0.36773559,  0.44975769,  0.51700393], 'code_length') # 0-1 units
            center = YTArray( [  9193.38986742,  11243.94226565,  12925.09833421])#, 'code_length') # 1 kpccm/h units for Arepo
            #bv    = YTArray( [ 7.50038669,  132.1431487 , -221.68740748], 'km/s' )  #bulk velocity of galaxy (everything)
            bv     = YTArray( [  45.6186475 ,  221.67851557, -107.29292011], 'km/s' )  #bulk velocity of galaxy (only gas)            
            Lgal   = YTArray( [ 1.9462109626e+29, -1.00050846265e+29, -8.76505777129e+28], 'cm**2/s') #ang mom vector for galaxy

        if redshift==2:
            #center = pf.arr([0.36773559,  0.44975769,  0.51700393], 'code_length') # 0-1 units
            center = YTArray( [  9372.12303426,  11350.02978143,  12671.74072015])#, 'code_length') # 1 kpccm/h units for Arepo
            ##bv   = YTArray( [ ], 'km/s' )  #bulk velocity of galaxy (everything)
            bv     = YTArray( [  77.86907977,   52.82035534, -146.35911295], 'km/s' )  #bulk velocity of galaxy (only gas)            
            Lgal   = YTArray( [ -4.07187308717e+27, 1.04049839857e+29, -1.45969665202e+28], 'cm**2/s') #ang mom vector for galaxy
    #-------------------------------------------------------- FIRE -----------------------------------------------------------------#
    if not center.any() or not bv.any() or not Lgal.any():
        try:
            simulation2="saved_data/"+simulation+"_z"+str(abs(round(redshift, 2)))
            fname=simulation2+"_haloparameters.txt"
            data=np.zeros(9).reshape(3, 3)
            names, trash, data[:,0], data[:,1], data[:,2] = np.loadtxt(fname, dtype='S4,S4,f4,f4,f4', unpack=True)
            center = YTArray(data[0,:]) # 'code_length' # 0-1 units
            bv     = YTArray(data[1,:], 'cm/s')
            Lgal   = YTArray(data[2,:], 'cm**2/s' )
        except:
            print "IMPORTANT NOTE: halo parameter file DOES NOT EXIST"
            
    if not center.any():
        print "IMPORTANT NOTE: halo center not defined yet.\nPlease run value, center = pf.h.find_max('all_density') to determine halo center\n"
    if not bv.any():
        print "IMPORTANT NOTE: halo bulk velocity not defined yet.\nPlease run get_halo_parameters(halo,pf)\n"
    if not Lgal.any():
        print "IMPORTANT NOTE: galaxy angular momentum not defined yet.\nPlease run get_halo_parameters(halo,pf)\n"

    return center, bv, Lgal

#define relative velocities by the bulk velocity
def _gas_relative_velocity_x(field, data):
    bulk_velocity = data.get_field_parameter("bulk_velocity")
    if bulk_velocity == None:
        bulk_velocity = YTArray( [0., 0., 0.], 'cm/s')
    return data['gas','velocity_x'] - bulk_velocity[0]

def _gas_relative_velocity_y(field, data):
    bulk_velocity = data.get_field_parameter("bulk_velocity")
    if bulk_velocity == None:
        bulk_velocity = YTArray( [0., 0., 0.], 'cm/s')
    return data['gas','velocity_y'] - bulk_velocity[1]

def _gas_relative_velocity_z(field, data):
    bulk_velocity = data.get_field_parameter("bulk_velocity")
    if bulk_velocity == None:
        bulk_velocity = YTArray( [0., 0., 0.], 'cm/s')
    return data['gas','velocity_z'] - bulk_velocity[2]

#Gadget/Arepo derived field for Temperature based on gas internal energy and electron abundance
def _gadget_temperature(field, data):
    from yt.utilities.physical_constants import kb,mp
    gamma = 5.0/3.0
    H_MASSFRAC = 0.76
    mu = 4.0 / (3.0 * H_MASSFRAC + 1.0 + 4.0 * H_MASSFRAC * data[('gas','ElectronAbundance')])
    #yt error reads the internal energy per unit mass as dimensionless.  Add proper dimensions.
    u = data.ds.arr(data[('gas','InternalEnergy')], 'code_velocity**2')
    return u*(gamma-1)*mu*mp/kb

#Art aliases adn derived fields 
def _art_temperature(field, data):
    from yt.utilities.physical_constants import kb,mp
    gamma = 5.0/3.0
    H_MASSFRAC = 0.76
    mu = 4.0 / (3.0 * H_MASSFRAC + 1.0 )  #No ElectronAbundance variable defined so don't take it into account
    #mu = 4.0 / (3.0 * H_MASSFRAC + 1.0 + 4.0 * H_MASSFRAC * data[('gas','ElectronAbundance')])  
    #yt error reads the internal energy per unit mass as dimensionless.  Add proper dimensions.
    u = data.ds.arr(data[('art','GasEnergy')], 'code_velocity**2')
    return u*(gamma-1)*mu*mp/kb

def _art_gas_velocity_x(field, data):
    mom=data['gas','momentum_x']
    mass=data['gas','cell_mass']
    return mom/mass
def _art_gas_velocity_y(field, data):
    mom=data['gas','momentum_y']
    mass=data['gas','cell_mass']
    return mom/mass
def _art_gas_velocity_z(field, data):
    mom=data['gas','momentum_z']
    mass=data['gas','cell_mass']
    return mom/mass

#Radial velocity bug in code: re-define to fix.
def _gas_radial_velocity(field, data):
    ptype="gas"
    normal = data.get_field_parameter('normal')
    center = data.get_field_parameter('center')
    bv = data.get_field_parameter("bulk_velocity")
    pos = "particle_position_%s"
    pos = YTArray([data[ptype, pos % ax].in_units('cm') for ax in "xyz"], 'cm')
    vel = "particle_velocity_%s"
    vel = YTArray([data[ptype, vel % ax].in_units('km/s') for ax in "xyz"], 'km/s')
    pos = np.reshape(pos, (3, pos.size/3))
    pos = pos - np.reshape(center, (3, 1))
    vel = np.reshape(vel, (3, vel.size/3))
    vel = vel - np.reshape(bv, (3, 1))
    theta = get_sph_theta(pos, normal)
    phi = get_sph_phi(pos, normal)
    sphr = get_sph_r_component(vel, theta, phi, normal)
    return sphr

#Radial velocity bug in code: re-define to fix.
def _particle_radial_velocity(field, data):
    ptype="all"
    normal = data.get_field_parameter('normal')
    center = data.get_field_parameter('center')
    bv = data.get_field_parameter("bulk_velocity")
    pos = "particle_position_%s"
    pos = YTArray([data[ptype, pos % ax].in_units('cm') for ax in "xyz"], 'cm')
    vel = "particle_velocity_%s"
    vel = YTArray([data[ptype, vel % ax].in_units('km/s') for ax in "xyz"], 'km/s')
    pos = np.reshape(pos, (3, pos.size/3))
    pos = pos - np.reshape(center, (3, 1))
    vel = np.reshape(vel, (3, vel.size/3))
    vel = vel - np.reshape(bv, (3, 1))
    theta = get_sph_theta(pos, normal)
    phi = get_sph_phi(pos, normal)
    sphr = get_sph_r_component(vel, theta, phi, normal)
    return sphr

#Create aliases to map SPH simulation quantities to the "usual" names (for grid codes)
def _sph_alias_specific_angular_momentum_x(field,data):
    return data['all','particle_specific_angular_momentum_x']
def _sph_alias_specific_angular_momentum_y(field,data):
    return data['all','particle_specific_angular_momentum_y']
def _sph_alias_specific_angular_momentum_z(field,data):
    return data['all','particle_specific_angular_momentum_z']
def _sph_alias_gas_specific_angular_momentum_x(field,data):
    return data['gas','particle_specific_angular_momentum_x']
def _sph_alias_gas_specific_angular_momentum_y(field,data):
    return data['gas','particle_specific_angular_momentum_y']
def _sph_alias_gas_specific_angular_momentum_z(field,data):
    return data['gas','particle_specific_angular_momentum_z']
def _sph_alias_radius(field,data):
    return data['all','particle_radius']
def _sph_alias_gas_radius(field,data):
    return data['gas','particle_radius']
def _sph_alias_radial_velocity(field,data):
    return data['all','particle_radial_velocity']
def _sph_alias_velocity_x(field,data):
    return data['all','particle_velocity_x']
def _sph_alias_velocity_y(field,data):
    return data['all','particle_velocity_y']
def _sph_alias_velocity_z(field,data):
    return data['all','particle_velocity_z']
def _sph_alias_gas_velocity_x(field,data):
    return data['gas','particle_velocity_x']
def _sph_alias_gas_velocity_y(field,data):
    return data['gas','particle_velocity_y']
def _sph_alias_gas_velocity_z(field,data):
    return data['gas','particle_velocity_z']

def _grid_alias_gas_radius(field,data):
    return data['index','spherical_radius']


#pf.add_field(('gas', "radial_velocity2"), function=_ENZOTEST_gas_radial_velocity, particle_type=False, units="km/s", force_override=True)
def _ENZOTEST_gas_radial_velocity(field, data):
    normal = data.get_field_parameter('normal')
    center = data.get_field_parameter('center')
    bv = data.get_field_parameter("bulk_velocity")
    #pos = "%"
    pos = YTArray([data['index', ax].in_units('cm') for ax in "xyz"], 'cm')
    vel = "velocity_%s"
    vel = YTArray([data['gas', vel % ax].in_units('km/s') for ax in "xyz"], 'km/s')
    pos = np.reshape(pos, (3, pos.size/3))
    pos = pos - np.reshape(center, (3, 1))
    vel = np.reshape(vel, (3, vel.size/3))
    vel = vel - np.reshape(bv, (3, 1))
    theta = get_sph_theta(pos, normal)
    phi = get_sph_phi(pos, normal)
    sphr = get_sph_r_component(vel, theta, phi, normal)
    return sphr

        
#Define "normal", "north" and "east" Cutting Plane Velocities based on angular momentum vector of galaxy
def _gas_offaxis_velocity_north(field, data):
    bv   = data.get_field_parameter('bulk_velocity')
    Lgal = data.get_field_parameter("angular_momentum_vector").in_units('cm**2/s')
    face_on = Lgal/YTQuantity(np.linalg.norm(Lgal), 'cm**2/s')
    vec = np.array([0,0,1])
    north_vector = np.cross(vec, face_on)
    north_vector = north_vector / np.linalg.norm(north_vector)
    ptype='gas'
    v = north_vector[0]*(data[ptype,'velocity_x']-bv[0]) + north_vector[1]*(data[ptype,'velocity_y']-bv[1]) + north_vector[2]*(data[ptype,'velocity_z']-bv[2])
    return v

def _gas_offaxis_velocity_east(field, data):
    bv   = data.get_field_parameter('bulk_velocity')
    Lgal = data.get_field_parameter("angular_momentum_vector").in_units('cm**2/s')
    face_on = Lgal/YTQuantity(np.linalg.norm(Lgal), 'cm**2/s')
    vec = np.array([0,0,1])
    north_vector = np.cross(vec, face_on)
    north_vector = north_vector / np.linalg.norm(north_vector)
    edge_on = np.cross(north_vector,face_on)
    edge_on = edge_on / np.linalg.norm(edge_on) 
    ptype='gas'
    v =      edge_on[0]*(data[ptype,'velocity_x']-bv[0]) +      edge_on[1]*(data[ptype,'velocity_y']-bv[1]) +      edge_on[2]*(data[ptype,'velocity_z']-bv[2])
    return v

def _gas_offaxis_velocity_normal(field, data):
    bv   = data.get_field_parameter('bulk_velocity')
    Lgal = data.get_field_parameter("angular_momentum_vector").in_units('cm**2/s')
    face_on = Lgal/YTQuantity(np.linalg.norm(Lgal), 'cm**2/s')
    ptype='gas'
    v =      face_on[0]*(data[ptype,'velocity_x']-bv[0]) +      face_on[1]*(data[ptype,'velocity_y']-bv[1]) +      face_on[2]*(data[ptype,'velocity_z']-bv[2])
    return v

#pf.add_field(('gas','offaxis_velocity_north'),  function=_gas_offaxis_velocity_north, units='km/s', particle_type=False, take_log=False)
#pf.add_field(('gas','offaxis_velocity_east'),   function=_gas_offaxis_velocity_east,  units='km/s', particle_type=False, take_log=False)
#pf.add_field(('gas','offaxis_velocity_normal'), function=_gas_offaxis_velocity_normal,units='km/s', particle_type=False, take_log=False)
