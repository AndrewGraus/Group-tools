#This is a simple script that looks over all the gadget
#binaries in the target folder and if they don't have
#an associated split program, make it
#
#***************NOTE******************
#
#this is hard coded to only run the 5to4 split so that needs
#to be modified if I'm going to use something else

import os, sys, re, subprocess

split_files=[]

#loop over files in directory and if they have a .split.gadget
#end then add them to a list that says they are already split
for import_file in os.listdir(sys.argv[1]):
    im_file_split = re.split('[.]',import_file)
    if str(im_file_split[-1])=='gadget' and str(im_file_split[-2])=='split':
        split_files.append(import_file)

#do the same thing as above but check to see if it has a .gadget
#but NOT a .split then run the split program
for import_file in os.listdir(sys.argv[1]):
    im_file_split = re.split('[.]',import_file)
    if im_file_split[-1]=='gadget' and import_file not in split_files:
        subprocess.call(['python','/home/agraus/code/python_generalscripts/split_gb_5to4.py',str(sys.argv[1])+import_file,str(sys.argv[1])+im_file_split[0]+'.split.gadget'])
        #print str(sys.argv[1])+im_file_split[0]+'.split.gadget', str(sys.argv[1])+import_file


        
