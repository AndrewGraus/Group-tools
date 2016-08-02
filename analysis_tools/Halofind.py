#assuming an IRATE file with AHF/Rockstar and shit in it
#Arguments are the file, then ahf or rockstar, then z if 

#!/usr/bin/python

import sys

import h5py
import numpy

from optparse import OptionParser
usage = "usage: python %prog <IRATE_File> <halo_catalog> [options]" #maks sure vizlength is in units you want (i.e. kpc if you set --usekpc, etc.)
parser = OptionParser(usage=usage)

parser.add_option("-z",dest="zoom",action="store_true",help="Set if you're using AHF on a zoom-in simulation")
parser.add_option("-n",dest="neighbor",action="store_true",help="Set if you want to find the nearest more massive halos as well")

(ops,args) = parser.parse_args()

if len(args) < 2:
    parser.print_help()
    sys.exit(1337)

h = 0.71

massmin = 9.8 * 10**9
massmax = 1.5 * 10**10
rcut=100.
Vcut=10.
#the mass range desired

isofact = 3 #multiple of rvir for large things we want halo outside of
isoself = 5 #multiple of rvir for halo that we want large things outside of

mconv=1e10 #mass conversion from gadget to m_sun

goodhalos = [] #list of halos that match the mass and isolation 1 criteria
goodhalos2 = [] #list of halos that match the mass and isolation 1 criteria
nearestneighbor = []

try:
    f = h5py.File(args[0])
except IOError:
    print "{0} doesn't appear to be an HDF5 format file.  Please provide a valid IRATE file.".format(sys.argv[1])

try:
    cat_type = args[1]

except IOError:
    print "please specify what halo catalog to use".format(sys.argv[1])

simprop=f['SimulationProperties']

snaps = []
for key in f.keys():      
    if key.startswith('Snapshot'):  snaps.append(key)

snap = f[str(snaps[-1])]#select the last snapshot
print snap.keys()
part=snap['ParticleData']
hal=part['Dark_Halo']
masshr=hal['Mass'][:]*mconv
print max(masshr)

if cat_type == "AHF":
    cat=snap['HaloCatalog_AHF']
    allc = cat['cNFW'][:]
    allspin = cat['lambda'][:]
    allhost = cat['hostHalo'][:]
if cat_type == "Rockstar":
    cat=snap['HaloCatalog_Rockstar']
    allrs = cat['Rs'][:]/h
    allspin = cat['spin_bullock'][:]

hrx = hal['Position'][:,0]*1000/h
hry = hal['Position'][:,1]*1000/h
hrz = hal['Position'][:,2]*1000/h

allmass = cat['Mvir'][:]/h
allcen = cat['Center'][:]*1000/h
allr = cat['Rvir'][:]/h
allid = cat['ID'][:]
alln = cat['npart'][:]
allvmax=cat['Vmax'][:]
allrmax=cat['Rmax'][:]

if cat_type == "Rockstar":
    allc = allr/allrs

#Now everything is in kpc and M_sun

print "There are {0} halos in the box".format(len(allmass))

#print min(allid)

if ops.zoom:#Rockstar only halo finds on HR particles so only need this for AHF
    boundary=part['Dark_Boundary']
    masslr=boundary['Mass'][:]*mconv
    lrx = boundary['Position'][:,0]*1000/h
    lry = boundary['Position'][:,1]*1000/h
    lrz = boundary['Position'][:,2]*1000/h
    #select only HR halos
    allFmass = cat['fMhires'][:]
    hr = allFmass >= .999
    hrmass = allmass[hr]
    hrcen = allcen[hr]
    hrr = allr[hr]
    hrid = allid[hr]
    hrn = alln[hr]
    hrhost = allhost[hr]
    hrrm = allrmax[hr]
    hrvmax=allvmax[hr]
    hrconc=allc[hr]
    hrspin=allspin[hr]
    #only necessary if halofinding on a zoom simulation.  Else use below:
    print "{0} hi-res halos found".format(len(hrmass))
else:
    hrmass = allmass
    hrcen = allcen
    hrr = allr
    hrid = allid
    hrn = alln
    hrrm = allrmax
    hrvmax=allvmax
    hrconc=allc
    hrspin=allspin
    if cat_type == "AHF":
        hrhost = allhost

if cat_type == "AHF":
    #next is to make sure we're not looking at subhalos:
    primary  = hrhost == -1.0
    Mass = hrmass[primary]
    Cen = hrcen[primary]
    Rad = hrr[primary]
    IDs = hrid[primary]
    npart = hrn[primary]
    concentration=hrconc[primary]
    spin=hrspin[primary]
    rmax=hrrm[primary]
    vmax=hrvmax[primary]
    print "{0} hr primary systems found".format(len(Mass))

if cat_type == "Rockstar":    
    #Rockstar implementation to come once I figure out how find_parent works.  For now we'll assume the isolation cut handles it
    Mass = hrmass
    Cen = hrcen
    Rad = hrr
    IDs = hrid
    npart = hrn
    concentration=hrconc
    spin=hrspin
    rmax=hrrm
    vmax=hrvmax
    print "{0} hr primary systems found".format(len(Mass))


#Then examine mass range
bigenough = Mass >= massmin
smallenough = Mass < massmax

Mint = Mass[bigenough & smallenough]
Cint = Cen[bigenough & smallenough]
Rint = Rad[bigenough & smallenough]
IDint = IDs[bigenough & smallenough]
npartint = npart[bigenough & smallenough]
concint = concentration[bigenough & smallenough]
spinint = spin[bigenough & smallenough]
rmint = rmax[bigenough & smallenough]
vmint = vmax[bigenough & smallenough]
print '{0} hires. primary halos found in mass range'.format(len(Mint))

allhalodata=numpy.column_stack((IDint, Mint, concint, spinint, rmint, vmint))
numpy.savetxt('rock_allhalos.txt',allhalodata)

#Then implement isolation criterion

xpos = allcen[:,0]
ypos = allcen[:,1]
zpos = allcen[:,2]

nIDint=[]

for i in range(len(IDint)):
    
    #print 'checking halo with Rvir {0}'.format(Rint[i])

    halid = IDint[i]
    notself = allid != halid
    otherhals = allid[notself]
    masscheck = allmass[notself]
    otherr = allr[notself]

    cen=Cint[i]

    distvects = allcen[notself]-cen
    dx,dy,dz=distvects[:,0],distvects[:,1],distvects[:,2]

    distances = (dx**2+dy**2+dz**2)**0.5

    massmask = (masscheck>Mint[i])&(distances >0)

    centers=allcen[notself][massmask]
    rmask=allr[notself][massmask]

    nearest = distances==min(distances[massmask])
    print 'clostest big halo is ', min(distances[massmask]), 'mpc away'
    print 'with distance cut = ', isofact*otherr[nearest]
    nn = otherhals[nearest][0]
    print 'closest halo is ', nn
    
    isolated=isofact*rmask < distances[massmask]
    
    if isolated.all():
        goodhalos.append(IDint[i])
        nearestneighbor.append(nn)
    #else:
    #    print 'Halo failed: it is too close to something big'
    
    
    xmin = Cint[i][0]-isoself*Rint[i]
    xmax = Cint[i][0]+isoself*Rint[i]
    ymin = Cint[i][1]-isoself*Rint[i]
    ymax = Cint[i][1]+isoself*Rint[i]
    zmin = Cint[i][2]-isoself*Rint[i]
    zmax = Cint[i][2]+isoself*Rint[i]
    
    xgood1 = xpos > xmin 
    xgood2 = xpos < xmax
    xgood = xgood1 & xgood2
    ygood1 = ypos > ymin 
    ygood2 = ypos < ymax
    ygood = ygood1 & ygood2
    zgood1 = zpos > zmin 
    zgood2 = zpos < zmax
    zgood = zgood1 & zgood2
    inrange = xgood & ygood & zgood

    #halid = IDint[i]
    
    notself2 = allid[inrange] != halid
    
    massrange =  allmass[inrange]
    massrange_ns = massrange[notself2]
    for halom in massrange_ns:
        if halom >= Mint[i] * 0.5:
            #print "halo {0} failed".format(IDint[i])
            #print "halo mass is", Mint[i], "and found a halo with mass ", halom
            break
    else:
        goodhalos2.append(IDint[i])
    
print nearestneighbor

halwant = sorted(set(goodhalos).intersection(goodhalos2))

print len(halwant), " halos found that match all criteria.  They are:"
print halwant

#print "halos found at ", allcen[halwant]

hwindex=[]
for hnum in halwant:
    hwindex.append(int(hnum))

Halocens=allcen[hwindex]

cenx=Halocens[:,0]
ceny=Halocens[:,1]
cenz=Halocens[:,2]

neardens=[]
lowR=allrmax[hwindex]<=rcut
goodV=allvmax[hwindex][lowR]>Vcut
highspin = allspin[hwindex][lowR][goodV]>(numpy.mean(allspin[hwindex])-numpy.std(allspin[hwindex]))
lowspin = allspin[hwindex][lowR][goodV]<(numpy.mean(allspin[hwindex])+numpy.std(allspin[hwindex]))
goodspin = lowspin & highspin

isofiledata=numpy.column_stack((allid[hwindex], allmass[hwindex], allc[hwindex], allspin[hwindex], allrmax[hwindex], allvmax[hwindex]))
numpy.savetxt('isohalos.txt',isofiledata)

isocutdata=numpy.column_stack((allid[hwindex][lowR][goodV][goodspin], allmass[hwindex][lowR][goodV][goodspin], allc[hwindex][lowR][goodV][goodspin], allspin[hwindex][lowR][goodV][goodspin], allrmax[hwindex][lowR][goodV][goodspin], allvmax[hwindex][lowR][goodV][goodspin]))
numpy.savetxt('isohalos2.txt',isocutdata)

print len(allid[hwindex][lowR][goodV][goodspin])
print allid[hwindex][lowR][goodV][goodspin]

gcenx=allcen[hwindex][lowR][goodV][goodspin][:,0]
gceny=allcen[hwindex][lowR][goodV][goodspin][:,1]
gcenz=allcen[hwindex][lowR][goodV][goodspin][:,2]

gid2=allid[hwindex][lowR][goodV][goodspin]



if ops.neighbor:
    for i in range(len(halwant)):
        index=halwant[i]
        print index
        print allmass[index]
        print allcen[index]
        
        bighal = allmass>allmass[index]
        bigcen=allcen[bighal]
        bigm=allmass[bighal]
        bigid=allid[bighal]
        cendiff = allcen[bighal]-allcen[index]
        dist=(sum(cendiff**2,axis=1))**0.5
        neighbor=dist==min(dist)
        print "found a larger halo this far: ", dist
        print bigm[neighbor]
        print bigid[neighbor]
        print bigcen[neighbor]
