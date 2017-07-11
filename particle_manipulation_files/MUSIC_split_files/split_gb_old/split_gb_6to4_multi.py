#!/bin/python

import struct,sys
from numpy import fromstring,array,append,empty,unique

def read_gbin(fname,verbose=False):
    print "\nOpening "+fname    
    f = open(fname,'rb')

    #First read the header so I know how many particles of each type I have
    header_size = struct.unpack('<I',f.read(4))[0]

    print header_size

    #number of particles of each type in this file
    nfile = struct.unpack('<6I',f.read(24)) #Number of particles in this file

    print nfile

    masstable = struct.unpack('<6d',f.read(48))  #masses of the particle groups
        
    print masstable
    masstable = (masstable[0],0.0,0.0,0.0,0.0,0.0)

    a = struct.unpack('<d',f.read(8))[0]        #expansion factor
    z = struct.unpack('<d',f.read(8))[0]        #redshift

    flag_sfr = struct.unpack('<i',f.read(4))[0] #star formation included?
    flag_feed = struct.unpack('<i',f.read(4))[0] #feedback included?

    ntot = struct.unpack('<6i',f.read(24))      #total number of particles in the simulation (= nfile if numfiles == 1)

    print ntot

    flag_cool = struct.unpack('<i',f.read(4))[0]  #cooling included?
    numfiles = struct.unpack('<i',f.read(4))[0]   #number of files in each snapshot
    boxsize = struct.unpack('<d',f.read(8))[0] #Size of the box, if periodic
    omega0 = struct.unpack('<d',f.read(8))[0]  #matter density at z = 0
    omegaL = struct.unpack('<d',f.read(8))[0]  #vacuum energy density at z = 0
    h = struct.unpack('<d',f.read(8))[0] #hubble parameter in units of 100 km/s/Mpc
    flag_age = struct.unpack('<i',f.read(4))[0]  #stellar age included?
    flag_metals = struct.unpack('<i',f.read(4))[0]  #use metals?
    nhighword = struct.unpack('<6i',f.read(24))   #contains the most significant word of 64-bit particle numbers (if npart > 2^32)

    flag_entropy = struct.unpack('<i',f.read(4))[0] #entropy instead of thermal energy in initial conditions?

    f.seek(264,0)   #Moves to the end of the header (and block that tells you size of header)

    header = [nfile,masstable,a,z,flag_sfr,flag_feed,ntot,flag_cool,numfiles,boxsize,omega0,omegaL,h,flag_age,flag_metals,nhighword,flag_entropy]
    
    gas,halo,disk,bulge,star,bndry = [],[],[],[],[],[]
    ngas = nfile[0]
    nhalo = nfile[1]
    ndisk = nfile[2]
    nbulge = nfile[3]
    nstar = nfile[4]
    nbndry = nfile[5]
    
    print ngas, nhalo, ndisk, nbulge, nstar, nbndry

    print numfiles

    print "Reading coordinates"
    coord_size = struct.unpack('<I',f.read(4))[0]

    gas.append(fromstring(f.read(12*ngas),dtype='f').reshape((-1,3)))
    halo.append(fromstring(f.read(12*nhalo),dtype='f').reshape((-1,3)))
    disk.append(fromstring(f.read(12*ndisk),dtype='f').reshape((-1,3)))
    print "Read coordinates for {0} of {1} particles".format(ngas+nhalo+ndisk,sum(nfile))
    bulge.append(fromstring(f.read(12*nbulge),dtype='f').reshape((-1,3)))
    star.append(fromstring(f.read(12*nstar),dtype='f').reshape((-1,3)))
    bndry.append(fromstring(f.read(12*nbndry),dtype='f').reshape((-1,3)))
    
    if struct.unpack('<I',f.read(4))[0] != coord_size:
        raise StandardError("The block size at the end of the coordinate block doesn't match that at the beginning.  This is an issue.")

    #next up is velocities.  pretty identical to the coordinates.
    vel_size = struct.unpack('<I',f.read(4))[0]
    print "Reading velocities"

    gas.append(fromstring(f.read(12*ngas),dtype='f').reshape((-1,3)))
    halo.append(fromstring(f.read(12*nhalo),dtype='f').reshape((-1,3)))
    disk.append(fromstring(f.read(12*ndisk),dtype='f').reshape((-1,3)))
    print "Read velocities for {0} of {1} particles".format(ngas+nhalo+ndisk,sum(nfile))
    bulge.append(fromstring(f.read(12*nbulge),dtype='f').reshape((-1,3)))
    star.append(fromstring(f.read(12*nstar),dtype='f').reshape((-1,3)))
    bndry.append(fromstring(f.read(12*nbndry),dtype='f').reshape((-1,3)))

    if struct.unpack('<I',f.read(4))[0] != vel_size:   #And read the size of the block again.
        raise StandardError("The block size at the end of the velocity block doesn't match that at the beginning.  This is an issue.")

    #Next up are the particle IDs, which are unsigned integers.
    id_size = struct.unpack('<I',f.read(4))[0]
    print "Reading particle IDs"
    gas.append(fromstring(f.read(4*ngas),dtype='I'))
    halo.append(fromstring(f.read(4*nhalo),dtype='I'))
    disk.append(fromstring(f.read(4*ndisk),dtype='I'))
    print "Read IDs for {0} of {1} particles".format(ngas+nhalo+ndisk,sum(nfile))
    bulge.append(fromstring(f.read(4*nbulge),dtype='I'))
    star.append(fromstring(f.read(4*nstar),dtype='I'))
    bndry.append(fromstring(f.read(4*nbndry),dtype='I'))

    if struct.unpack('<I',f.read(4))[0] != id_size:   #And read the size of the block again.
        raise StandardError("The block size at the end of the IDs block doesn't match that at the beginning.  This is an issue.")

    gasmass = False
    #Now I have to do the mass block.  Do this by checking if there are particles in a block that have a zero in the mass table.    
    if (ngas > 0 and masstable[0] == 0) or (nhalo > 0 and masstable[1] == 0) or (ndisk > 0 and masstable[2] == 0) or (nbulge > 0 and masstable[3] == 0) or (nstar > 0 and masstable[4] == 0) or (nbndry > 0 and masstable[5] == 0):
    #In other words, only read the size of the mass block if there is a mass block (for any of the groups)
        mass_size = struct.unpack('<I',f.read(4))[0]
        
        if ngas > 0 and masstable[0] == 0:    #There are particles in the group, but their masses aren't in the header (so they must be in the file)
            print "Reading variable masses for gas group"
            gas.append(fromstring(f.read(4*ngas),dtype='f'))
            gasmass = True
            
        if nhalo > 0 and masstable[1] == 0:
            print "Reading variable masses for halo group"
            halo.append(fromstring(f.read(4*nhalo),dtype='f'))
        if ndisk > 0 and masstable[2] == 0:
            print "Reading variable masses for disk group"
            disk.append(fromstring(f.read(4*ndisk),dtype='f'))
        if nbulge > 0 and masstable[3] == 0:
            print "Reading variable masses for bulge group"
            bulge.append(fromstring(f.read(4*nbulge),dtype='f'))
        if nstar > 0 and masstable[4] == 0:
            print "Reading variable masses for star group"
            star.append(fromstring(f.read(4*nstar),dtype='f'))
        if nbndry > 0 and masstable[5] == 0:
            print "Reading variable masses for boundary group"
            bndry.append(fromstring(f.read(4*nbndry),dtype='f'))
        
        if struct.unpack('<I',f.read(4))[0] != mass_size:   #And read the size of the block again.
            raise StandardError("The block size at the end of the mass block doesn't match that at the beginning.  This is an issue.")


    if ngas > 0:
        print "Reading gas specific data."
        
        #Put all this inside the if statement because I don't want it reading block sizes for gas blocks if there is no gas data (cause then there will be no block size and I will accidentally read into other data or past the end of the file.)
        
        #Internal energy:
        
        u_size = struct.unpack('<I',f.read(4))[0]
        print ngas, 4*ngas, u_size
        gas.append(fromstring(f.read(4*ngas),dtype='f'))
        print struct.unpack('<I',f.read(4))[0], u_size
        if struct.unpack('<I',f.read(4))[0] != u_size:
            raise StandardError("The block size at the end of the internal energy block doesn't match that at the beginning.  This is an issue.")
        print "Read gas internal energy; not attempting to read gas density or smoothing length because this is assumed to be an IC file"
        
        
    current_pos = f.tell()
    f.seek(0,2) #Jump to the end of the file
    if current_pos == f.tell():
        print "Read "+fname
    else:
        print "Completed reading "+fname+" but there remain {0} bytes at the end of the file unread.".format(f.tell()-current_pos)
    f.close()

    return [header,gas,halo,disk,bulge,star,bndry,gasmass]


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print "Usage:  python split_gb.py <input file base> <output directory>"
        sys.exit(1337)

    inbase = sys.argv[1]
    outdir = sys.argv[2]
    if not outdir.endswith('/'):    outdir = outdir + '/'
    
    from shutil import copyfile
    from os import path,mkdir
    if not path.isdir(outdir):
        print "Making output directory {0}".format(outdir)
        mkdir(outdir)
        
    from glob import glob
    files = glob(inbase+'*')
    files.sort()
    outfiles = []
    
    bndrymasses = []
    bulgemasses = []
    diskmasses = []
        
    if len(files) == 0:
        print "No files found."
        sys.exit(1337)
    elif len(files) == 1:
        print "Only one file found; please use the non-multi version of this script for that."
        sys.exit(1337)
        
    for fname in files:
        name = fname.split('/')[-1]
        stripped = name.rstrip('0123456789')
        num = name[len(stripped):]
        outname = outdir + stripped.rstrip('.')+"_split."+num
        outfiles.append(outname)
        
        [header,gasdata,halodata,diskdata,bulgedata,stardata,bndrydata,gasmass] = read_gbin(fname,verbose=True)

        if len(diskdata[0]) > 0 or len(bulgedata[0]) > 0:
            print "There are already particles in either disk or bulge.  Exiting."
            sys.exit(1337)

        #My goal is to leave gasdata, halo, and stars unchanged, but split the bndry group into disk, bulge, and bndry based on masses
        #Scratch that--in this case, I want to leave gas, halo, and stars unchanged, but split the boundary group into disk, then the next 2 into bulge
        #This is gonna be hardcoded like a motherfucker, so be warned.
        
        if len(bndrydata[0]) == 0:
            print "There are no particles in the boundary group.  Copying {0} to {1}".format(fname,outname)
            copyfile(fname,outname)
            continue

        allbndrypos = bndrydata[0]
        allbndryvel = bndrydata[1]
        allbndryid = bndrydata[2]
        allbndrymass = bndrydata[3]
        
        print "Finding masses in file."
        masses = unique(allbndrymass)    
        masses.sort()

        if len(masses) != 5:
            print "There are {0} different masses in the boundary group; this script was no more than 5.  Do it by hand.".format(len(masses))
            sys.exit(1337)

        level11mass = masses[0]
        level10mass = masses[1]
        level9mass = masses[2]
        level8mass = masses[3]
        level7mass = masses[4]

        diskpos = allbndrypos[allbndrymass == level11mass]
        diskvel = allbndryvel[allbndrymass == level11mass]
        diskID = allbndryid[allbndrymass == level11mass]
        
        #OK, that's the disk group.  Now bulge is going to be 9 and 10.  Also, DO THEY SPEAK ENGLISH IN WHAT MOTHERFUCKER?
        #Also, that means that I have to do the mass block for that group.
        
        l10pos = allbndrypos[allbndrymass == level10mass]
        l10vel = allbndryvel[allbndrymass == level10mass]
        l10ID = allbndryid[allbndrymass == level10mass]
        l10mass = empty(len(l10pos),dtype='f')
        l10mass[:] = level10mass
        
        nl10 = len(l10pos)
        
        l9pos = allbndrypos[allbndrymass == level9mass]
        l9vel = allbndryvel[allbndrymass == level9mass]
        l9ID = allbndryid[allbndrymass == level9mass]
        l9mass = empty(len(l9pos),dtype='f')
        l9mass[:] = level9mass
        
        nl9 = len(l9pos)
        
        #Now to join them together 

        bulgepos = append(l10pos,l9pos,axis=0)
        bulgevel = append(l10vel,l9vel,axis=0)
        bulgeID = append(l10ID,l9ID)
        bulgemass = append(l10mass,l9mass)
        
        assert len(bulgemass) == nl10 + nl9
        assert len(bulgepos) == nl10 + nl9  #Make sure I didn't fuck a duck
        assert len(bulgepos) == len(bulgemass)
        
        #Free up some memory
        del l10pos,l10vel,l10ID,l10mass,l9pos,l9vel,l9ID,l9mass
        
        #Now do the same shit for levels 7 and 8
        l8pos = allbndrypos[allbndrymass == level8mass]
        l8vel = allbndryvel[allbndrymass == level8mass]
        l8ID = allbndryid[allbndrymass == level8mass]
        l8mass = empty(len(l8pos),dtype='f')
        l8mass[:] = level8mass
        
        nl8 = len(l8pos)
        
        l7pos = allbndrypos[allbndrymass == level7mass]
        l7vel = allbndryvel[allbndrymass == level7mass]
        l7ID = allbndryid[allbndrymass == level7mass]
        l7mass = empty(len(l7pos),dtype='f')
        l7mass[:] = level7mass

        nl7 = len(l7pos)

        bndrypos = append(l8pos,l7pos,axis=0)
        bndryvel = append(l8vel,l7vel,axis=0)
        bndryID = append(l8ID,l7ID)
        bndrymass = append(l8mass,l7mass)
        
        assert len(bndrymass) == nl8 + nl7
        assert len(bndrypos) ==  nl8 + nl7  #Make sure I'm not still fucking a chicken
        assert len(bndrypos) == len(bndrymass)
        
        #Free up some memory
        del l8pos,l8vel,l8ID,l8mass,l7pos,l7vel,l7ID,l7mass

        ndisk = len(diskpos)
        nbulge = len(bulgepos)
        nbndry = len(bndrypos)

        if gasmass:  print "There are {0} particles in the gas group with variable masses".format(header[0][0])
        else:        print "There are {0} particles in the gas group with a mass of {1:.4g}".format(header[0][0],float(header[1][0]))
        print "There are {0} particles in the halo group with a mass of {1:.4g}".format(header[0][1],float(header[1][1]))
        print "Placed {0} particles in the disk group with a mass of {1:.4g}".format(ndisk,float(level11mass))
        print "Placed {0} particles in the bulge group with a masss of {1:.4g} and {2} particles with a mass of {3:.4g}".format(nl10,float(level10mass),nl9,float(level9mass))
        print "There are {0} particles in the star group with a mass of {1:.4g}".format(header[0][4],float(header[1][4]))
        print "Placed {0} particles in the bndry group with a mass of {1:.4g} and {2} particles with a mass of {3:.4g}".format(nl8,float(level8mass),nl7,float(level7mass))

        ntot = ndisk+nbulge+nbndry+header[0][1]+header[0][0]+header[0][4]
        gas = header[0][0] > 0

        assert ndisk+nbulge+nbndry == header[0][5]      #Check that all the particles that were in the bndry group originally made it somewhere

        #pack the header and write it:
        
        f = open(outname,'wb')
        f.write(struct.pack('<I',256))
        f.write(struct.pack('<6I',header[0][0],header[0][1],ndisk,nbulge,header[0][4],nbndry))    #Npart_file
        f.write(struct.pack('<6d',header[1][0],header[1][1],level11mass,0,header[1][4],0)) #Mass table 
        f.write(struct.pack('<d',header[2]))   #time
        f.write(struct.pack('<d',header[3]))   #z
        f.write(struct.pack('<i',header[4]))   #FlagSfr
        f.write(struct.pack('<i',header[5]))   #FlagFeedback
        f.write(struct.pack('<6i',header[6][0],header[6][1],ndisk,nbulge,header[6][4],nbndry)) #Npart_total
        f.write(struct.pack('<i',header[7]))   #FlagCooling
        f.write(struct.pack('<i',header[8]))   #NumFiles
        f.write(struct.pack('<d',header[9]))   #BoxSize
        f.write(struct.pack('<d',header[10]))  #Omega0
        f.write(struct.pack('<d',header[11]))  #OmegaL
        f.write(struct.pack('<d',header[12]))  #h
        f.write(struct.pack('<i',header[13]))  #FlagAge
        f.write(struct.pack('<i',header[14]))  #FlagMetals
        f.write(struct.pack('<6i',header[15][0],header[15][1],header[15][2],header[15][3],header[15][4],header[15][5]))    #NallHW -- might have to modify this for large sims
        f.write(struct.pack('<i',header[16]))  #flag_entr_ics

        header_bytes_left = 260 - f.tell()
        for i in range(header_bytes_left):
            f.write(struct.pack('<x'))    #Fill out the rest of the header with pad bytes
        assert f.tell() == 260
        f.write(struct.pack('<I',256))
        
        #Now pack the data using tostring, but I have to check that the dtypes of all my arrays are right
        gasdata[0] = gasdata[0].astype('f')
        gasdata[1] = gasdata[1].astype('f')
        gasdata[2] = gasdata[2].astype('I')
        if gas:
            gasdata[3] = gasdata[3].astype('f')      #If gasmass, this is the mass block, otherwise, it's the internal energy
            if gasmass:
                gasdata[4] = gasdata[4].astype('f')  #This is the internal energy if it exists, but it only exists if gasmass is true
        
        halodata[0] = halodata[0].astype('f')
        halodata[1] = halodata[1].astype('f')
        halodata[2] = halodata[2].astype('I')
        try:
            halodata[3] = halodata[3].astype('f')
        except IndexError:
            print "There is no mass block for the halo group, as expected."

        bndrypos = bndrypos.astype('f')
        bndryvel = bndryvel.astype('f')
        bndryID = bndryID.astype('I')
        bndrymass = bndrymass.astype('f')
        
        bulgepos = bulgepos.astype('f')
        bulgevel = bulgevel.astype('f')
        bulgeID = bulgeID.astype('I')
        bulgemass = bulgemass.astype('f')
        
        diskpos = diskpos.astype('f')
        diskvel = diskvel.astype('f')
        diskID = diskID.astype('I')
        
        stardata[0] = stardata[0].astype('f')
        stardata[1] = stardata[1].astype('f')
        stardata[2] = stardata[2].astype('I')
        try:
            stardata[3] = stardata[3].astype('f')
        except IndexError:
            print "There is no mass block for star particles, as expected."
            
        #Now I have to write it all to the file
        print "Writing coordinates."
        f.write(struct.pack('<I',12*ntot))
        f.write(gasdata[0].tostring())
        f.write(halodata[0].tostring())
        f.write(diskpos.tostring())
        f.write(bulgepos.tostring())
        f.write(stardata[0].tostring())
        f.write(bndrypos.tostring())
        f.write(struct.pack('<I',12*ntot))
        
        
        print "Writing velocities."
        f.write(struct.pack('<I',12*ntot))
        f.write(gasdata[1].tostring())
        f.write(halodata[1].tostring())
        f.write(diskvel.tostring())
        f.write(bulgevel.tostring())
        f.write(stardata[1].tostring())
        f.write(bndryvel.tostring())
        f.write(struct.pack('<I',12*ntot))
        
        
        print "Writing particles IDs."
        f.write(struct.pack('<I',4*ntot))
        f.write(gasdata[2].tostring())
        f.write(halodata[2].tostring())
        f.write(diskID.tostring())
        f.write(bulgeID.tostring())
        f.write(stardata[2].tostring())
        f.write(bndryID.tostring())
        f.write(struct.pack('<I',4*ntot))
        
        nmass = len(bulgemass)+len(bndrymass)
        if len(stardata) > 3:
            nmass = nmass + len(stardata[3])
        if len(halodata) > 3:
            nmass = nmass + len(halodata[3])
        if gasmass:
            nmass = nmass + len(gasdata[3])
        
        print "Writing masses."
        
        f.write(struct.pack('<I',4*nmass))
        if gasmass:
            print "Writing variable masses for gas particles."
            f.write(gasdata[3].tostring())
        if len(halodata) > 3:
            print "Writing variable masses for halo particles."
            f.write(halodata[3].tostring())
        #Don't give a shit about disk group--it's definitely one species
        f.write(bulgemass.tostring())
        if len(stardata) > 3:
            print "Writing variable masses for star particles."
            f.write(stardata[3].tostring())
        f.write(bndrymass.tostring())
        f.write(struct.pack('<I',4*nmass))
        
        
        if gas:
            print "Writing gas internal energy."
            #Then I have to do gas specific stuff
            f.write(struct.pack('<I',4*len(gasdata[0])))
            if gasmass:   f.write(gasdata[4].tostring())
            else:         f.write(gasdata[3].tostring())
            f.write(struct.pack('<I',4*len(gasdata[0])))
        
        f.close()
        
        print "Wrote "+outname
        
    print "Correcting nfile and making sure the mass table is correct."
    ntot = array([0,0,0,0,0,0])
    for fname in outfiles:
        f = open(fname,'rb')
        f.seek(4)
        nfile = fromstring(f.read(24),dtype='i')
        ntot = ntot + nfile

        f.close()
        
    for fname in outfiles:
        f = open(fname,'rb+')
        f.seek(4+24+48+8+8+4+4)    #skip the block size, nfile table, mass table, time, redshift, flagSfr, FlagFeedback, to get to the Nall table
        f.write(struct.pack('<6i',ntot[0],ntot[1],ntot[2],ntot[3],ntot[4],ntot[5]))
        f.close()
    
    
    
    





