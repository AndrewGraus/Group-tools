#!/usr/bin/env python
import struct
import numpy as np
from numpy import fromstring

def read_gbin(fname,verbose=False):
    print "\nOpening "+fname
    f = open(fname,'rb')

    #First read the header so I know how many particles of each type I have
    header_size = struct.unpack('<I',f.read(4))[0]

    #number of particles of each type in this file
    nfile = struct.unpack('<6I',f.read(24)) #Number of particles in this file

    masstable = struct.unpack('<6d',f.read(48))  #masses of the particle groups

    a = struct.unpack('<d',f.read(8))[0]        #expansion factor
    z = struct.unpack('<d',f.read(8))[0]        #redshift

    flag_sfr = struct.unpack('<i',f.read(4))[0] #star formation included?
    flag_feed = struct.unpack('<i',f.read(4))[0] #feedback included?

    ntot = struct.unpack('<6i',f.read(24))      #total number of particles in the simulation (= nfile if numfiles == 1)

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

    header = {'NumPart_ThisFile':nfile,'MassTable':masstable,'Time':a,'Redshift':z,'Flag_Sfr':flag_sfr,'Flag_Feedback':flag_feed,
    'NumPart_Total':ntot,'Flag_Cooling':flag_cool,'NumFilesPerSnapshot':numfiles,'BoxSize':boxsize,'Omega0':omega0,'OmegaLambda':omegaL,
    'HubbleParam':h,'Flag_StellarAge':flag_age,'Flag_Metals':flag_metals,'NumPart_Total_HW':nhighword,'Flag_Entropy_ICs':flag_entropy}

    # header = [nfile,masstable,a,z,flag_sfr,flag_feed,ntot,flag_cool,numfiles,boxsize,omega0,omegaL,h,flag_age,flag_metals,nhighword,flag_entropy]

    gas,halo,disk,bulge,star,bndry = {},{},{},{},{},{}
    ngas = nfile[0]
    nhalo = nfile[1]
    ndisk = nfile[2]
    nbulge = nfile[3]
    nstar = nfile[4]
    nbndry = nfile[5]

    print "Reading coordinates"
    coord_size = struct.unpack('<I',f.read(4))[0]

    gas['Coordinates'] = fromstring(f.read(12*ngas),dtype='f').reshape((-1,3))
    halo['Coordinates'] = fromstring(f.read(12*nhalo),dtype='f').reshape((-1,3))
    disk['Coordinates'] = fromstring(f.read(12*ndisk),dtype='f').reshape((-1,3))
    print "Read coordinates for {0} of {1} particles".format(ngas+nhalo+ndisk,sum(nfile))
    bulge['Coordinates'] = fromstring(f.read(12*nbulge),dtype='f').reshape((-1,3))
    star['Coordinates'] = fromstring(f.read(12*nstar),dtype='f').reshape((-1,3))
    bndry['Coordinates'] = fromstring(f.read(12*nbndry),dtype='f').reshape((-1,3))

    if struct.unpack('<I',f.read(4))[0] != coord_size:
        raise StandardError("The block size at the end of the coordinate block doesn't match that at the beginning.  This is an issue.")

    #next up is velocities.  pretty identical to the coordinates.
    vel_size = struct.unpack('<I',f.read(4))[0]
    print "Reading velocities"

    gas['Velocities'] = fromstring(f.read(12*ngas),dtype='f').reshape((-1,3))
    halo['Velocities'] = fromstring(f.read(12*nhalo),dtype='f').reshape((-1,3))
    disk['Velocities'] = fromstring(f.read(12*ndisk),dtype='f').reshape((-1,3))
    print "Read coordinates for {0} of {1} particles".format(ngas+nhalo+ndisk,sum(nfile))
    bulge['Velocities'] = fromstring(f.read(12*nbulge),dtype='f').reshape((-1,3))
    star['Velocities'] = fromstring(f.read(12*nstar),dtype='f').reshape((-1,3))
    bndry['Velocities'] = fromstring(f.read(12*nbndry),dtype='f').reshape((-1,3))

    if struct.unpack('<I',f.read(4))[0] != vel_size:   #And read the size of the block again.
        raise StandardError("The block size at the end of the velocity block doesn't match that at the beginning.  This is an issue.")

    #Next up are the particle IDs, which are unsigned integers.
    id_size = struct.unpack('<I',f.read(4))[0]
    print "Reading particle IDs"
    gas['ParticleIDs'] = fromstring(f.read(4*ngas),dtype='I')
    halo['ParticleIDs'] = fromstring(f.read(4*nhalo),dtype='I')
    disk['ParticleIDs'] = fromstring(f.read(4*ndisk),dtype='I')
    print "Read IDs for {0} of {1} particles".format(ngas+nhalo+ndisk,sum(nfile))
    bulge['ParticleIDs'] = fromstring(f.read(4*nbulge),dtype='I')
    star['ParticleIDs'] = fromstring(f.read(4*nstar),dtype='I')
    bndry['ParticleIDs'] = fromstring(f.read(4*nbndry),dtype='I')

    if struct.unpack('<I',f.read(4))[0] != id_size:   #And read the size of the block again.
        raise StandardError("The block size at the end of the IDs block doesn't match that at the beginning.  This is an issue.")

    #Now I have to do the mass block.  Do this by checking if there are particles in a block that have a zero in the mass table.
    if True in [(nfile[ii] > 0 and masstable[ii] == 0) for ii in range(len(nfile))]:
    #In other words, only read the size of the mass block if there is a mass block (for any of the groups)
        mass_size = struct.unpack('<I',f.read(4))[0]

        if ngas > 0 and masstable[0] == 0:    #There are particles in the group, but their masses aren't in the header (so they must be in the file)
            print "Reading variable masses for gas group"
            gas['Masses'] = fromstring(f.read(4*ngas),dtype='f')
        if nhalo > 0 and masstable[1] == 0:
            print "Reading variable masses for halo group"
            halo['Masses'] = fromstring(f.read(4*nhalo),dtype='f')
        if ndisk > 0 and masstable[2] == 0:
            print "Reading variable masses for disk group"
            disk['Masses'] = fromstring(f.read(4*ndisk),dtype='f')
        if nbulge > 0 and masstable[3] == 0:
            print "Reading variable masses for bulge group"
            bulge['Masses'] = fromstring(f.read(4*nbulge),dtype='f')
        if nstar > 0 and masstable[4] == 0:
            print "Reading variable masses for star group"
            star['Masses'] = fromstring(f.read(4*nstar),dtype='f')
        if nbndry > 0 and masstable[5] == 0:
            print "Reading variable masses for boundary group"
            bndry['Masses'] = fromstring(f.read(4*nbndry),dtype='f')

        if struct.unpack('<I',f.read(4))[0] != mass_size:   #And read the size of the block again.
            raise StandardError("The block size at the end of the mass block doesn't match that at the beginning.  This is an issue.")

    #Put all this inside the if statement because I don't want it reading block sizes for gas blocks if there is no gas data (cause then there will be no block size and I will accidentally read into other data or past the end of the file.)
    if ngas > 0:
        print "Reading gas specific data."
        #Internal energy:
        u_size = struct.unpack('<I',f.read(4))[0]
        gas['InternalEnergy'] = fromstring(f.read(4*ngas),dtype='f')
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

    return {'Header':header,'PartType0':gas,'PartType1':halo,'PartType2':disk,'PartType3':bulge,'PartType4':star,'PartType5':bndry}




if __name__ == '__main__':
    import sys
    from optparse import OptionParser

    parser = OptionParser(usage="%prog <inname> <outname> [options]")

    parser.add_option('--usebh',dest='usebh',help="Put particles in BH (type 5).  Otherwise, only puts particles in Disk and Bulge",default=False,action="store_true")
    # parser.add_option('--usestars',dest='usestars',help="Put particles in Stars (type 4).  Otherwise, only puts particles in Disk and Bulge",default=False,action="store_true")

    ops,args = parser.parse_args()

    try:
        inname,outname = args
    except Exception:
        parser.print_help()
        sys.exit(1)

    #first check if it's a multi-part file:
    from glob import glob
    files = glob(inname+'*')
    if len(files) > 1:
        print "Sorry, haven't written this yet!"
        sys.exit(1)

    data = read_gbin(inname)
    header = data['Header']
    #figure out which particle group we're splitting
    #know we're not splitting gas, halo, or star

    if header['NumPart_Total'][2] > 0:
        assert header['NumPart_Total'][3] == 0
        splitto = ['PartType2','PartType3']
        if ops.usebh:
            assert header['NumPart_Total'][5] == 0
            splitto.append('PartType5')
        tosplit = 'PartType2'
    elif header['NumPart_Total'][3] > 0:
        splitto = ['PartType2','PartType3']
        assert header['NumPart_Total'][2] == 0
        if ops.usebh: 
            assert header['NumPart_Total'][5] == 0
            splitto.append('PartType5')
        tosplit = 'PartType3'
    elif header['NumPart_Total'][5] > 0:
        splitto = ['PartType2','PartType3']
        assert header['NumPart_Total'][2] == 0
        assert header['NumPart_Total'][3] == 0
        tosplit = 'PartType5'
        if ops.usebh:
            splitto.append('PartType5')
    else:
        print "No particles found in any of the low res groups!"
        sys.exit(1)

    splitto.sort()

    nfile = np.array(header['NumPart_ThisFile'],copy=True)
    ntot = np.array(header['NumPart_Total'],copy=True)
    mtable = np.array(header['MassTable'],copy=True)
    assert (nfile == ntot).all()

    splitidx = int(tosplit[-1])
    nsplit = ntot[splitidx]
    masses = data[tosplit]['Masses']
    umasses = np.unique(masses)

    #ok, now to figure out where to put the particles....
    def chunkIt(seq, num):
        avg = len(seq) / float(num)
        out = []
        last = 0.0

        while last < len(seq):
            out.append(seq[int(last):int(last + avg)])
            last += avg
        lengths = np.array([len(kk) for kk in out])
        if (np.diff(lengths) < 0).any():
            for ii in range(len(out)-1):
                if len(out[ii]) > len(out[ii+1]):
                    temp1 = np.array(out[ii][:-1],copy=True)
                    temp2 = np.hstack( ( out[ii][-1], out[ii+1] ) )
                    out[ii] = temp1
                    out[ii+1] = temp2
        return out

    if len(umasses) < len(splitto):
        outmasses = np.empty(len(splitto))
        outmasses[:] = -1
        outmasses[:len(umasses)] = umasses
    else:
        outmasses = chunkIt(np.sort(umasses),len(splitto))

    print "\n###########"
    print "Splitting {0} particles masses in {1} to ".format(umasses.size,tosplit)+', '.join(splitto)
    nfile[int(tosplit[-1])] = 0     #clear the data in that group for now; will get filled back in below if necessary

    splitdata = {}
    datatosplit = data[tosplit]
    loc = 0
    for key in splitto:
        mydict = {}
        allii = int(key[-1])
        msk = np.empty(nsplit,dtype=bool)
        msk[:] = False
        for m in outmasses[loc]:
            msk = (msk)|(masses==m)
        nfile[allii] = np.count_nonzero(msk)
        if len(outmasses[loc]) == 1:
            mtable[allii] = outmasses[loc]
        else:
            mydict['Masses'] = masses[msk]
        mydict['Coordinates'] = datatosplit['Coordinates'][msk]
        mydict['Velocities'] = datatosplit['Velocities'][msk]
        mydict['ParticleIDs'] = datatosplit['ParticleIDs'][msk]
        splitdata[key] = mydict
        loc += 1

    for ii in range(len(splitto)):
        allii = int(splitto[ii][-1])
        string = '{0} particles with masses '.format(nfile[allii])
        for jj in range(len(outmasses[ii])):            
            string += "{0:.3g}".format(outmasses[ii][jj])
            if jj != len(outmasses[ii]) - 1:
                string += ", "
            else:
                string += " "
        string += "are going in "+splitto[ii]
        print string

    data[tosplit] = {'Coordinates':np.array([],dtype='float32'),'Velocities':np.array([],dtype='float32'),'ParticleIDs':np.array([],dtype='uint32')}  #clear out the data in the original group
    assert nfile.sum() == ntot.sum()        #check that I still have all the particles -- nfile has been modified; ntot hasn't

    print "###########\n"

    #first figure out the new header
    print "Opening "+outname+" to write"
    f = open(outname,'wb')
    f.write(struct.pack('<I',256))
    f.write(nfile.astype('uint32').tostring())
    f.write(mtable.astype('float64').tostring())
    f.write(struct.pack('<d',header['Time']))   #time
    f.write(struct.pack('<d',header['Redshift']))   #z
    f.write(struct.pack('<i',header['Flag_Sfr']))   #FlagSfr
    f.write(struct.pack('<i',header['Flag_Feedback']))   #FlagFeedback
    f.write(nfile.astype('int32').tostring())   #ntotal
    f.write(struct.pack('<i',header['Flag_Cooling']))   #FlagCooling
    f.write(struct.pack('<i',header['NumFilesPerSnapshot']))   #NumFiles
    f.write(struct.pack('<d',header['BoxSize']))   #BoxSize
    f.write(struct.pack('<d',header['Omega0']))  #Omega0
    f.write(struct.pack('<d',header['OmegaLambda']))  #OmegaL
    f.write(struct.pack('<d',header['HubbleParam']))  #h
    f.write(struct.pack('<i',header['Flag_StellarAge']))  #FlagAge
    f.write(struct.pack('<i',header['Flag_Metals']))  #FlagMetals
    f.write(np.array(header['NumPart_Total_HW'],dtype='int32').tostring())    #NallHW -- might have to modify this for large sims
    f.write(struct.pack('<i',header['Flag_Entropy_ICs']))  #flag_entr_ics

    header_bytes_left = 260 - f.tell()
    for i in range(header_bytes_left):
        f.write(struct.pack('<x'))    #Fill out the rest of the header with pad bytes
    assert f.tell() == 260
    f.write(struct.pack('<I',256))

    #now pack and write the data
    vec_size = np.array([12*nfile.sum()],dtype='uint32').tostring()
    f.write(vec_size)
    f.write(data['PartType0']['Coordinates'].astype('float32').tostring())
    print "Wrote {0} coordinates for PartType0".format(data['PartType0']['Coordinates'].shape[0])
    f.write(data['PartType1']['Coordinates'].astype('float32').tostring())
    print "Wrote {0} coordinates for PartType1".format(data['PartType1']['Coordinates'].shape[0])
    for k in splitto:
        if k == 'PartType5':    continue    #do this after stars
        f.write(splitdata[k]['Coordinates'].astype('float32').tostring())
        print "Wrote {0} coordinates for {1}".format(splitdata[k]['Coordinates'].shape[0],k)
    f.write(data['PartType4']['Coordinates'].astype('float32').tostring())
    print "Wrote {0} coordinates for PartType4".format(data['PartType4']['Coordinates'].shape[0])
    if ops.usebh:
        f.write(splitdata['PartType5']['Coordinates'].astype('float32').tostring())
        print "Wrote {0} coordinates for PartType5".format(splitdata['PartType5']['Coordinates'].shape[0])
    else:
        f.write(data['PartType5']['Coordinates'].astype('float32').tostring())
        print "Wrote {0} coordinates for PartType5".format(data['PartType5']['Coordinates'].shape[0])
    f.write(vec_size)


    f.write(vec_size)
    f.write(data['PartType0']['Velocities'].astype('float32').tostring())
    f.write(data['PartType1']['Velocities'].astype('float32').tostring())
    for k in splitto:
        if k == 'PartType5':    continue    #do this after stars
        f.write(splitdata[k]['Velocities'].astype('float32').tostring())
    f.write(data['PartType4']['Velocities'].astype('float32').tostring())
    if ops.usebh:
        f.write(splitdata['PartType5']['Velocities'].astype('float32').tostring())
    else:
        f.write(data['PartType5']['Velocities'].astype('float32').tostring())
    f.write(vec_size)

    scalar_size = np.array([4*nfile.sum()],dtype='uint32').tostring()
    f.write(scalar_size)
    f.write(data['PartType0']['ParticleIDs'].astype('uint32').tostring())
    f.write(data['PartType1']['ParticleIDs'].astype('uint32').tostring())
    for k in splitto:
        if k == 'PartType5':    continue    #do this after stars
        f.write(splitdata[k]['ParticleIDs'].astype('uint32').tostring())
    f.write(data['PartType4']['ParticleIDs'].astype('uint32').tostring())
    if ops.usebh:
        f.write(splitdata['PartType5']['ParticleIDs'].astype('uint32').tostring())
    else:
        f.write(data['PartType5']['ParticleIDs'].astype('uint32').tostring())
    f.write(scalar_size)

    nmass = 0
    for ii in [0,1,4]:
        if 'Masses' in data['PartType'+str(ii)].keys():
            nmass += nfile[ii]
    for key in splitdata.keys():
        if 'Masses' in splitdata[key].keys():
            nmass += splitdata[key]['Masses'].size
    if not ops.usebh:
        if 'Masses' in data['PartType5'].keys():
            nmass += nfile[5]

    mass_size = np.array([4*nmass],dtype='uint32').tostring()
    f.write(mass_size)
    if 'Masses' in data['PartType0'].keys():
        print "Writing varaible masses for PartType0"
        f.write(data['PartType0']['Masses'].astype('float32').tostring())
    if 'Masses' in data['PartType1'].keys():
        print "Writing varaible masses for PartType1"
        f.write(data['PartType1']['Masses'].astype('float32').tostring())
    for key in splitto:
        if key == 'PartType5':  continue
        if 'Masses' in splitdata[key].keys():
            print "Writing varaible masses for "+key
            f.write(splitdata[key]['Masses'].astype('float32').tostring())
    if 'Masses' in data['PartType4'].keys():
        print "Writing varaible masses for PartType4"
        f.write(data['PartType4']['Masses'].astype('float32').tostring())
    if ops.usebh:
        if 'Masses' in splitdata['PartType5'].keys():
            print "Writing variable masses for PartType5"
            f.write(splitdata['PartType5']['Masses'].astype('float32').tostring())
    else:
        if 'Masses' in data['PartType5'].keys():
            print "Writing variable masses for PartType5"
            f.write(data['PartType5']['Masses'].astype('float32').tostring())
    f.write(mass_size)

    if nfile[0] > 0:
        print "Writing gas internal energy"
        gas_size = np.array([nfile[0]],dtype='uint32').tostring()
        f.write(gas_size)
        f.write(data['PartType0']['InternalEnergy'].astype('float32').tostring())
        f.write(gas_size)

    f.close()
    print "Wrote "+outname
