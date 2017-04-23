#This is a simple script that looks over all the gadget
#binaries in the target folder and if they don't have
#an associated split program, make it
#
#***************NOTE******************
#
#this is hard coded to only run the 5to4 split so that needs
#to be modified if I'm going to use something else

import os, sys, re, subprocess

if len(sys.argv) < 3:
    print 'usage: python run_split_prog.py <directory files are in> <number of species 4,5,6>'

total_part = sys.argv[2]

split_files=[]

if total_part==4:
    path_to_split = '/home/agraus/code/Group-tools/particle_manipulation_files/MUSIC_split_files/split_gb.py'
elif total_part==5:
    path_to_split = '/home/agraus/code/Group-tools/particle_manipulation_files/MUSIC_split_files/split_gb_5to4.py'
elif total_part==6:
    path_to_split = '/home/agraus/code/Group-tools/particle_manipulation_files/MUSIC_split_files/split_gb_6to4.py'
elif total_part==7:
    path_to_split = '/home/agraus/code/Group-tools/particle_manipulation_files/MUSIC_split_files/split_gb_7to4.py'

#loop over files in directory and check if they are
#gadget binaries
for import_file in os.listdir(sys.argv[1]):
    im_file_split = re.split('[.]',import_file)
    if str(im_file_split[-1])=='gadget':
        split_files.append(import_file)

#do the same thing as above but check to make sure that the
#file isn't a split file nor does it have a corresponding
#split file
for import_file in split_files:
    im_file_split = re.split('[.]',import_file)
    #if im_file_split[-1]=='gadget' and import_file not in split_files:
    if im_file_split[-2] != 'split' and im_file_split[-2]+'.split.gadget' not in os.listdir(sys.argv[1]):
        subprocess.call(['python',path_to_split,str(sys.argv[1])+import_file,str(sys.argv[1])+im_file_split[0]+'.split.gadget'])


        
