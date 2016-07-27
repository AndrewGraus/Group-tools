#Calculates profiles for the largest high-res halo (unless another is specified).  If you use --pow (which requires visnap) radii below the calculated power radius will be saved as negative values as in AHF. -ODE

'''
Units:
densities in Msun/pc
vcirc in km/s
r in kpc
'''

#!/usr/bin/python

import sys
#from pylab import * # this import most of the stuff from the matplotlib library, see http://matplotlib.sourceforge.net/ for help using the plotting routines. It imports all numpy stuff as well.
import h5py
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
from numpy import *
import visnap

from optparse import OptionParser
usage = "usage: python %prog <IRATE_File> <Snapshot> [options]" #maks sure vizlength is in units you want (i.e. kpc if you set --usekpc, etc.)
parser = OptionParser(usage=usage)

parser.add_option("--halo",dest="hindex",type=int,help="The index of the halo you want to examine, if not the largest", default=-1)
parser.add_option("-r",dest="rockstar",action="store_true",help="Set if you want calculate based off Rockstar info")
parser.add_option("-a",dest="ahf",action="store_true",help="Set if you want calculate based off AHF info")
parser.add_option("--pow",dest="rpow",action="store_true",help="Set if you want calculate the Power radius using Visnap")
parser.add_option("--kpc",dest="kpc",action="store_true",help="Set if the positions are in kpc")


(ops,args) = parser.parse_args()
print len(args)
if len(args) < 2:
    parser.print_help()
    sys.exit(1337)

##### Declarations & Definitions #####

#the resolution we have for rotation curves in terms of the Power radius.  1.5 for CDM, .5 for SIDM.
Resolution = 0.5

rlist=[]
vrotlist=[]
sigvlist=[]
diffrholist=[]
rPowlist=[]
rhocrit = 277.545 #h^2 Msun/kpc^3
omega_m = 0.266
rho_b = rhocrit * omega_m
solmass=1.99 * 10**30 #convert solar masses to kg
kpc=3.24 * 10**-17 #convert km to kpc
G = 4.302e-6 #G in kpc * (km/s)^2 / solar masses
nbins = 70
h = 0.71

outbase = '_profile_information.txt'

#Miguel Rocha's code to calculate the Power radius of the sims:
if ops.rpow:
    def find_rpower(r, nenc, avg_dens):
        '''
        Find the radius at which the numerical two body scattering time-scale is
        comparable to the Hubble time.  Eq. 20 in power et al. 2003 is used.
    
        Input:
    
         r - radial bins
    
         nenc - number of enclosed particles at r
    
         avg_dens - average enclosed density at r
    
        Output:
    
         r_power - minimum resolved radius
    
         argRes - r[argRes] gives all r > r_power  
    
         argNotRes - r[argNotRes] gives all r < r_power
        '''
        
        rhoCrit = visnap.rho_crit
        argNotRes = argwhere(sqrt(200.0)*asfarray(nenc)*sqrt(rhoCrit/avg_dens)/(8.0*log(asfarray(nenc))) < 1.0)[:,0]
        r_power = r[argNotRes].max()
        argRes = argwhere(r > r_power)[:,0]
        return r_power, argRes, argNotRes


#Routine to calculate profiles for a snapshot and add them to the lists for saving
def savepoints(snap):

    global rlist
    global vrotlist
    global sigvlist
    global diffrholist
    #global rPowlist
               
    halodata = snap['HaloCatalog_Rockstar']
    hal=snap['ParticleData']['Dark_Halo']
    
    print index1
    halocen=halodata['Center'][index1]
    print halocen
    
    pos=hal['Position'][:]
    vel=hal['Velocity'][:]
    x,y,z=pos[:,0],pos[:,1],pos[:,2]
    Vx,Vy,Vz=vel[:,0],vel[:,1],vel[:,2]
    mass = hal['Mass'][:]*(10**10)/h
    
    x=x-halocen[0]
    y=y-halocen[1]
    z=z-halocen[2]
    r = (x**2 + y**2 + z**2)**0.5
    r = r * 1000/h
    
    rvir = halodata['Rvir'][index1]
    print rvir
    
    inhal = r < rvir/h
    
    rin = r[inhal]
    vx = Vx[inhal]
    vy = Vy[inhal]
    vz = Vz[inhal]
    m = mass[inhal]
    
    radbins = logspace(log10(rvir/200.),log10(rvir),num=nbins)
    rmid =  zeros(len(radbins))
    diffrho = zeros(len(radbins))     
    vdisp = zeros(len(radbins),dtype='float')
    DMMass = zeros(len(radbins),dtype='float')
    vrot = zeros(len(radbins),dtype='float') 
    sphereden = zeros(len(radbins), dtype='float')
    nenclosed = zeros(len(radbins), dtype='float')
    rMin=[]
    rMax=[]
    
    for i in range(len(radbins)):      
        if i == 0:     
            lowr = 0
        else:
            lowr = radbins[i-1]
        highr = radbins[i]
        largeenough = rin > lowr         
        smallenough = rin < highr        
        inshell = largeenough & smallenough
        rmid[i] = (lowr+highr)/2.0

        # Density 
        enclosedmasses = m[inshell]
        vol = (highr**3 - lowr**3)*4.0*pi/3.0       
        menc = sum(enclosedmasses)
        diffrho[i] = menc/(vol* 10**9)  
               
        # Total Velocity Dispersion
        meanVx = vx[inshell].mean()
        meanV2x = (vx[inshell]*vx[inshell]).mean()
        vxdisp = sqrt(meanV2x - meanVx*meanVx)

        meanVy = vy[inshell].mean()
        meanV2y = (vy[inshell]*vy[inshell]).mean()
        vydisp = sqrt(meanV2y - meanVy*meanVy)

        meanVz = vz[inshell].mean()
        meanV2z =  (vz[inshell]*vz[inshell]).mean()
        vzdisp = sqrt(meanV2z - meanVz*meanVz)

        vdisp[i] = sqrt(vxdisp**2 + vydisp**2 + vzdisp**2)

        # Rotational Velocity
        minsphere = sum(m[smallenough])
        nenclosed[i] = len(m[smallenough])
        vr2 = sqrt(G*minsphere/rmid[i])
        sphereden[i]=minsphere/((4.0*pi/3.0)*rmid[i]**3)
		
        vrot[i] = vr2
    if ops.rpow:
        rPow, rmin, rMax= find_rpower(rmid, nenclosed, sphereden)
        print "power radius found to be", rPow
	    #print "Diffrho", min(diffrho), max(diffrho)
    
        for i in range(len(rmid)):
            if rmid[i] < rPow:
                #print "found unconverged radius"
                rmid[i] = rmid[i]*-1

    #Save results
    diffrholist.append(diffrho)
    vrotlist.append(vrot)
    rlist.append(rmid)
    sigvlist.append(vdisp)
    if ops.rpow:
        rPowlist.append(rPow)
    print "finished a halo"

def find_biggest(snap):
    if ops.rockstar:
        halodata = snap['HaloCatalog_Rockstar']
        m=halodata['Mvir'][:]
        ids=halodata['ID'][:]
        return ids[m==max(m)][0]
    elif ops.ahf:
        halodata = snap['HaloCatalog_AHF']
        m1=halodata['Mvir'][:]
        ids=halodata['ID'][:]
        fmhr=halodata['fMhires'][:]
        hires=(fmhr==1.)
        mhr=m1[hires]
        hid2=ids[hires]
        return hid2[mhr==max(mhr)][0]
    else:
        print "please specify the halofinder you wish to use"
        sys.exit(10203)

#####################   Star Main Program ########################

try:
    f = h5py.File(sys.argv[1])
except IOError:
    print "{0} doesn't appear to be an HDF5 format file.  Please provide a valid IRATE file.".format(sys.argv[1])

snapnum = int(sys.argv[2]) #reads the number of the snapshot as an integer    
snap = f['Snapshot{0:05}'.format(snapnum)]

if ops.hindex==-1:
    index1=find_biggest(snap)
else:
    index1=ops.hindex
print 'Begining analysis of halo ', index1 , ' in snapshot ', snapnum

savepoints(snap)

fname = 'bin_halo_' + str(index1) + '_snap_' + sys.argv[2] + outbase

print "Saving data as ", fname

#File formatted as Radius, DiffDen, Vcirc, VDisp with units of kpc, solar masses per pc^3 and km/s 
savetxt(fname, zip(rlist[0], diffrholist[0], vrotlist[0], sigvlist[0]))

#and save a separate file with the power radius:
print "Saving power radius..."
f = open('Power_Radius.txt',"w")
f.write("%s" %rPowlist[0])
f.close()

print "done!"
