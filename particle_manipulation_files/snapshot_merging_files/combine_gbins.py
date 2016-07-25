#!/bin/python

#from readgb import read_gbin
import sys
from glob import glob
from struct import pack,unpack
from shutil import copyfileobj
from os import remove
from numpy import array
import numpy

#I'm going to do a file for each block, then at the end concatenate all those files together
def readgb_as_binary(fname):
    """Reads the blocks of the gadget binary file <fname>, but only unpacks the
    header. Returns a list [header,gas blocks, halo blocks,...], where the blocks
    are in order coordinates, velocities, IDs, (mass, gas blocks)."""

    print "\nOpening "+fname    
    f = open(fname,'rb')

    #First read the header so I know how many particles of each type I have
    header_size = unpack('<I',f.read(4))[0]

    #number of particles of each type in this file
    nfile = unpack('<6I',f.read(24)) #Number of particles in this file

    masstable = unpack('<6d',f.read(48))  #masses of the particle groups
        
    a = unpack('<d',f.read(8))[0]        #expansion factor
    z = unpack('<d',f.read(8))[0]        #redshift

    flag_sfr = unpack('<i',f.read(4))[0] #star formation included?
    flag_feed = unpack('<i',f.read(4))[0] #feedback included?

    ntot = unpack('<6i',f.read(24))      #total number of particles in the simulation (= nfile if numfiles == 1)
        
    flag_cool = unpack('<i',f.read(4))[0]  #cooling included?
    numfiles = unpack('<i',f.read(4))[0]   #number of files in each snapshot
    boxsize = unpack('<d',f.read(8))[0] #Size of the box, if periodic
    omega0 = unpack('<d',f.read(8))[0]  #matter density at z = 0
    omegaL = unpack('<d',f.read(8))[0]  #vacuum energy density at z = 0
    h = unpack('<d',f.read(8))[0] #hubble parameter in units of 100 km/s/Mpc
    flag_age = unpack('<i',f.read(4))[0]  #stellar age included?
    flag_metals = unpack('<i',f.read(4))[0]  #use metals?
    nhighword = unpack('<6i',f.read(24))   #contains the most significant word of 64-bit particle numbers (if npart > 2^32)

    flag_entropy = unpack('<i',f.read(4))[0] #entropy instead of thermal energy in initial conditions?

    f.seek(264,0)   #Moves to the end of the header (and block that tells you size of header)

    header = [nfile,masstable,a,z,flag_sfr,flag_feed,ntot,flag_cool,numfiles,boxsize,omega0,omegaL,h,flag_age,flag_metals,nhighword,flag_entropy]
    
    gas,halo,disk,bulge,star,bndry = [],[],[],[],[],[]
    ngas = nfile[0]
    nhalo = nfile[1]
    ndisk = nfile[2]
    nbulge = nfile[3]
    nstar = nfile[4]
    nbndry = nfile[5]
    
    bsize = f.read(4)
    gas.append(f.read(12*ngas))
    halo.append(f.read(12*nhalo))
    disk.append(f.read(12*ndisk))
    bulge.append(f.read(12*nbulge))
    star.append(f.read(12*nstar))
    bndry.append(f.read(12*nbndry))
    if f.read(4) != bsize:
        raise StandardError("The block size at the end of the coordinate block doesn't match that at the beginning.  This is an issue.")
        
    bsize = f.read(4)
    gas.append(f.read(12*ngas))
    halo.append(f.read(12*nhalo))
    disk.append(f.read(12*ndisk))
    bulge.append(f.read(12*nbulge))
    star.append(f.read(12*nstar))
    bndry.append(f.read(12*nbndry))
    if f.read(4) != bsize:
        raise StandardError("The block size at the end of the velocity block doesn't match that at the beginning.  This is an issue.")
    
    bsize = f.read(4)
    gas.append(f.read(4*ngas))
    halo.append(f.read(4*nhalo))
    disk.append(f.read(4*ndisk))
    bulge.append(f.read(4*nbulge))
    star.append(f.read(4*nstar))
    bndry.append(f.read(4*nbndry))
    if f.read(4) != bsize:
        raise StandardError("The block size at the end of the ID block doesn't match that at the beginning.  This is an issue.")
        
    mtrue = array([False,False,False,False,False,False])
    if (ngas > 0 and masstable[0] == 0) or (nhalo > 0 and masstable[1] == 0) or (ndisk > 0 and masstable[2] == 0) or (nbulge > 0 and masstable[3] == 0) or (nstar > 0 and masstable[4] == 0) or (nbndry > 0 and masstable[5] == 0):
    #In other words, only read the size of the mass block if there is a mass block (for any of the groups)
        bsize = f.read(4)
        if ngas > 0 and masstable[0] == 0:    
            #There are particles in the group, but their masses aren't in the header (so they must be in the file)
            print "Grabbing variable masses for gas group"
            gas.append(f.read(4*ngas))
            mtrue[0] = True
        if nhalo > 0 and masstable[1] == 0:
            print "Grabbing variable masses for halo group"
            halo.append(f.read(4*nhalo))
            mtrue[1] = True
        if ndisk > 0 and masstable[2] == 0:
            mtrue[2] = True
            print "Grabbing variable masses for disk group"
            disk.append(f.read(4*ndisk))
        if nbulge > 0 and masstable[3] == 0:
            mtrue[3] = True
            print "Grabbing variable masses for bulge group"
            bulge.append(f.read(4*nbulge))
        if nstar > 0 and masstable[4] == 0:
            print "Grabbing variable masses for star group"
            star.append(f.read(4*nstar))
            mtrue[4] = True
        if nbndry > 0 and masstable[5] == 0:
            print "Grabbing variable masses for boundary group"
            bndry.append(f.read(4*nbndry))
            mtrue[5] = True
        if f.read(4) != bsize:   #And read the size of the block again.
            raise StandardError("The block size at the end of the mass block doesn't match that at the beginning.  This is an issue.")


    if ngas > 0:
        bsize = f.read(4)
        gas.append(f.read(4*ngas))
        if bsize != f.read(4):
            raise StandardError("The block size at the end of the internal energy block doesn't match that at the beginning.  This is an issue.")

        bsize = f.read(4)
        gas.append(f.read(4*ngas))
        if bsize != f.read(4):
            raise StandardError("The block size at the end of the density block doesn't match that at the beginning.  This is an issue.")        

        bsize = f.read(4)
        gas.append(f.read(4*ngas))
        if bsize != f.read(4):
            raise StandardError("The block size at the end of the smoothing length block doesn't match that at the beginning.  This is an issue.")
            
    current_pos = f.tell()
    f.seek(0,2) #Jump to the end of the file
    if current_pos == f.tell():
        print "Read "+fname
    else:
        print "Completed reading "+fname+" but there remain {0} bytes at the end of the file unread.".format(f.tell()-current_pos)
    f.close()
    
    return [header,gas,halo,disk,bulge,star,bndry,mtrue]

if len(sys.argv) < 3:
    print "Usage:  python combine_gbins.py <input file base> <output file>"
    print "Combines portions of snapshots (e.g. snapshot_001.0 and snapshot_001.1"
    print "into one file, which is specified.  Provide the filename up to the number"
    print "after the final dot."
    sys.exit(1337)
    
    
files = glob(sys.argv[1]+'*')
files.sort()

if len(files) == 0:
    print "No files found!  Exiting."
    sys.exit(1337)
else:
    print "Will combine {0} files into {1}".format(len(files),sys.argv[2])

headerfile = 'tmp_header.bin'
headerf = open(headerfile,'wb')

gascoordfile = 'tmp_gascoords.bin'
gascoords = open(gascoordfile,'wb')
gasvelfile = 'tmp_gasvels.bin'
gasvels = open(gasvelfile,'wb')
gasidfile = 'tmp_gasids.bin'
gasids = open(gasidfile,'wb')
gasmassfile = 'tmp_gasmass.bin'
gasmass = open(gasmassfile,'wb')
ufile = 'tmp_intenergy.bin'
u = open(ufile,'wb')
rhofile = 'tmp_rho.bin'
rho = open(rhofile,'wb')
hsmlfile = 'tmp_hsml.bin'
hsml = open(hsmlfile,'wb')

halocoordfile = 'tmp_halocoords.bin'
halocoords = open(halocoordfile,'wb')
halovelfile = 'tmp_halovels.bin'
halovels = open(halovelfile,'wb')
haloidfile = 'tmp_haloids.bin'
haloids = open(haloidfile,'wb')
halomassfile = 'tmp_halomass.bin'
halomass = open(halomassfile,'wb')

diskcoordfile = 'tmp_diskcoords.bin'
diskcoords = open(diskcoordfile,'wb')
diskvelfile = 'tmp_diskvels.bin'
diskvels = open(diskvelfile,'wb')
diskidfile = 'tmp_diskids.bin'
diskids = open(diskidfile,'wb')
diskmassfile = 'tmp_diskmass.bin'
diskmass = open(diskmassfile,'wb')

bulgecoordfile = 'tmp_bulgecoords.bin'
bulgecoords = open(bulgecoordfile,'wb')
bulgevelfile = 'tmp_bulgevels.bin'
bulgevels = open(bulgevelfile,'wb')
bulgeidfile = 'tmp_bulgeids.bin'
bulgeids = open(bulgeidfile,'wb')
bulgemassfile = 'tmp_bulgemass.bin'
bulgemass = open(bulgemassfile,'wb')

starcoordfile = 'tmp_starcoords.bin'
starcoords = open(starcoordfile,'wb')
starvelfile = 'tmp_starvels.bin'
starvels = open(starvelfile,'wb')
staridfile = 'tmp_starids.bin'
starids = open(staridfile,'wb')
starmassfile = 'tmp_starmass.bin'
starmass = open(starmassfile,'wb')

bndrycoordfile = 'tmp_bndrycoords.bin'
bndrycoords = open(bndrycoordfile,'wb')
bndryvelfile = 'tmp_bndryvels.bin'
bndryvels = open(bndryvelfile,'wb')
bndryidfile = 'tmp_bndryids.bin'
bndryids = open(bndryidfile,'wb')
bndrymassfile = 'tmp_bndrymass.bin'
bndrymass = open(bndrymassfile,'wb')



nmass = 0

for i in range(len(files)):
    fname = files[i]
    [header,gas,halo,disk,bulge,star,bndry,mtrue] = readgb_as_binary(fname)
    if i == 0:  ntot = sum(header[6])
    if ntot >= 4294967296:  #(2**32--the move files a gadget snapshot can hold)
        print "The snapshots contain more particles than can be represented in a single file.  Exiting."
        sys.exit(1337)
        
    nmass = nmass + sum(array(header[0])[mtrue])  #mtrue only tells me for a given file, so I can't base this all on the first file, so I'll have to jump back in the mass file to fix the block header
    
    if i == 0:
        #vecsize = pack('<I',12*ntot)
        #scalarsize = pack('<I',4*ntot)
        #gassize = pack('<I',4*header[6][0])
        vecsize = array([12*ntot],dtype='I').tostring()
        scalarsize = array([4*ntot],dtype='I').tostring()
        gassize = array([4*header[6][0]],dtype='I').tostring()
        print "Writing block sizes."
        gascoords.write(vecsize)    #Write to the gas files, since they'll always come first
        gasvels.write(vecsize)
        gasids.write(scalarsize)
        gasmass.write(pack('<I',0))   #have to go back and change this later, but that's easy enough
        if header[6][0] > 0:
            u.write(gassize)
            rho.write(gassize)
            hsml.write(gassize)
                
        #Now write the data to file
        #header only gets done by the first file
        
        print "Writing header."
        headerf.write(pack('<I',256))
        headerf.write(pack('<6I',header[6][0],header[6][1],header[6][2],header[6][3],header[6][4],header[6][5]))
        headerf.write(pack('<6d',header[1][0],header[1][1],header[1][2],header[1][3],header[1][4],header[1][5]))    #Hopefully this isn't wrong--that is, hopefully the masses will be written in the header of the first file if they exist anywhere
        headerf.write(pack('<d',header[2]))
        headerf.write(pack('<d',header[3]))
        headerf.write(pack('<i',header[4]))
        headerf.write(pack('<i',header[5]))
        headerf.write(pack('<6i',header[6][0],header[6][1],header[6][2],header[6][3],header[6][4],header[6][5]))
        headerf.write(pack('<i',header[7]))
        headerf.write(pack('<i',1))
        headerf.write(pack('<d',header[9]))
        headerf.write(pack('<d',header[10]))
        headerf.write(pack('<d',header[11]))
        headerf.write(pack('<d',header[12]))
        headerf.write(pack('<i',header[13]))
        headerf.write(pack('<i',header[14]))
        headerf.write(pack('<6i',header[15][0],header[15][1],header[15][2],header[15][3],header[15][4],header[15][5]))
        headerf.write(pack('<i',header[16]))
        header_bytes_left = 260 - headerf.tell()
        
        for i in range(header_bytes_left):
            headerf.write(pack('<x'))    #Fill out the rest of the header with pad bytes
        assert headerf.tell() == 260
        headerf.write(pack('<I',256))
        
    #Now I want to write the blocks to their files            
    gascoords.write(gas[0])
    gasvels.write(gas[1])
    gasids.write(gas[2])
    if header[0][0] > 0:
        if mtrue[0]:
            gasmass.write(gas[3])
            u.write(gas[4])
            rho.write(gas[5])
            hsml.write(gas[6])
        else:
            u.write(gas[3])
            rho.write(gas[4])
            hsml.write(gas[5])
    #Now I'm done with gas, so let's free up that memory
    del gas
    
    #Now repeat all that for the other particle groups
    halocoords.write(halo[0])
    halovels.write(halo[1])
    haloids.write(halo[2])
    if mtrue[1]:
        halomass.write(halo[3])       
    del halo
        
    diskcoords.write(disk[0])
    diskvels.write(disk[1])
    diskids.write(disk[2])
    if mtrue[2]:
        diskmass.write(disk[3])     
    del disk
        
    bulgecoords.write(bulge[0])
    bulgevels.write(bulge[1])
    bulgeids.write(bulge[2])
    if mtrue[3]:
        bulgemass.write(bulge[3])
    del bulge
        
    starcoords.write(star[0])
    starvels.write(star[1])
    starids.write(star[2])
    if mtrue[4]:
        starmass.write(star[3])
    del star
        
    bndrycoords.write(bndry[0])
    bndryvels.write(bndry[1])
    bndryids.write(bndry[2])
    if mtrue[5]:
        bndrymass.write(bndry[3])
    del bndry
    
    print "Finished moving blocks from "+fname
    
#Now let's close the blocks by writing their sizes again
bndrycoords.write(vecsize)      #This time, write to bndry files, since they'll always come last.
bndryvels.write(vecsize)
bndryids.write(scalarsize)
u.write(gassize)
rho.write(gassize)
hsml.write(gassize)

#For mass I have to go back and fill in the correct block size
masssize = pack('<I',4*nmass)
bndrymass.write(masssize)
gasmass.seek(0)
gasmass.write(masssize)

headerf.close()
gascoords.close()
gasvels.close()
gasids.close()
gasmass.close()
u.close()
rho.close()
hsml.close()

halocoords.close()
halovels.close()
haloids.close()
halomass.close()

diskcoords.close()
diskvels.close()
diskids.close()
diskmass.close()

bulgecoords.close()
bulgevels.close()
bulgeids.close()
bulgemass.close()

starcoords.close()
starvels.close()
starids.close()
starmass.close()

bndrycoords.close()
bndryvels.close()
bndryids.close()
bndrymass.close()

#Now to cat all these files together
print "\nCatting together the blocks for different datasets"
master = open(sys.argv[2],'wb')
copyfileobj(open(headerfile,'rb'),master)

#Now all the coordinates
copyfileobj(open(gascoordfile,'rb'),master)
copyfileobj(open(halocoordfile,'rb'),master)
copyfileobj(open(diskcoordfile,'rb'),master)
copyfileobj(open(bulgecoordfile,'rb'),master)
copyfileobj(open(starcoordfile,'rb'),master)
copyfileobj(open(bndrycoordfile,'rb'),master)

#And velocities
copyfileobj(open(gasvelfile,'rb'),master)
copyfileobj(open(halovelfile,'rb'),master)
copyfileobj(open(diskvelfile,'rb'),master)
copyfileobj(open(bulgevelfile,'rb'),master)
copyfileobj(open(starvelfile,'rb'),master)
copyfileobj(open(bndryvelfile,'rb'),master)

copyfileobj(open(gasidfile,'rb'),master)
copyfileobj(open(haloidfile,'rb'),master)
copyfileobj(open(diskidfile,'rb'),master)
copyfileobj(open(bulgeidfile,'rb'),master)
copyfileobj(open(staridfile,'rb'),master)
copyfileobj(open(bndryidfile,'rb'),master)

copyfileobj(open(gasmassfile,'rb'),master)
copyfileobj(open(halomassfile,'rb'),master)
copyfileobj(open(diskmassfile,'rb'),master)
copyfileobj(open(bulgemassfile,'rb'),master)
copyfileobj(open(starmassfile,'rb'),master)
copyfileobj(open(bndrymassfile,'rb'),master)

copyfileobj(open(ufile,'rb'),master)
copyfileobj(open(rhofile,'rb'),master)
copyfileobj(open(hsmlfile,'rb'),master)

print "Wrote "+sys.argv[2]+".  Cleaning up now."
remove(headerfile)
remove(gascoordfile)
remove(gasvelfile)
remove(gasidfile)
remove(gasmassfile)
remove(ufile)
remove(rhofile)
remove(hsmlfile)

remove(halocoordfile)
remove(halovelfile)
remove(haloidfile)
remove(halomassfile)

remove(diskcoordfile)
remove(diskvelfile)
remove(diskidfile)
remove(diskmassfile)

remove(bulgecoordfile)
remove(bulgevelfile)
remove(bulgeidfile)
remove(bulgemassfile)

remove(starcoordfile)
remove(starvelfile)
remove(staridfile)
remove(starmassfile)

remove(bndrycoordfile)
remove(bndryvelfile)
remove(bndryidfile)
remove(bndrymassfile)



