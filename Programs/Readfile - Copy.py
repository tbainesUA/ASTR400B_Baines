
# coding: utf-8

# In[27]:


"""Import Packages"""
import numpy as np
import astropy.units as u


# In[37]:


"""Loadfile Function"""
# Define a function that will read in a file
# Inputs: fname(without extension i.e. .txt, .dat, etc)
# Returns: Time (Myr), total number of particles and array with data
# CALLING FUNCTION: time, particles_tot (total), data = Read(filename) 

def Read(filename):
    
    #Set path to find where data files are stored
    path = "C:\\Users\Tyler\Documents\ASTR 400B\ASTR400B_Baines\Data\\"+filename+".txt"
    #print (path) #uncomment for print check
   
    #open file
    file = open(path, 'r')
    
    #read header info line by line (line will be a string)
    # read first two lines FIRST and store as variable
    
    # read and store time 
    line1 = file.readline()
    label, value = line1.split()
    time = float(value)*10.0*u.Myr
    
    # read and store total number of particles
    line2 = file.readline()
    label, value = line2.split()
    particles_tot = float(value)
    
    #close the file
    file.close()
   

    # read the remainder of the file, 
    # "dtype=None" means line is split using white spaces
    # "skip_header=3"  skipping the first 3 lines 
    # the flag "names=True" creates arrays to store the date
    #       with the column headers given in line 4 like "m", "x"
    data = np.genfromtxt(path, dtype = None, names = True, skip_header = 3)
    
    # this will return the time of the snapshot, 
    #total number of particles 
    #and an array that stores the remainder of the data
    return time, particles_tot, data

