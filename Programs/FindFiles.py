
# coding: utf-8

# In[1]:


# This programs find the specfic files based off the time stamped one wishes to investigate
# For Astr400 Research Project MW/M31 Galaxy Merger Simulation


# In[6]:


# Import packages
import numpy as np 
import astropy.units as u

from Readfile import Read


# In[45]:


def FileFinder(Galaxy, Time_array):
    # By Tyler Baines and Drew Weldon
    # For Dr. Gurtina Besla ASTR 400B class Research Project
    # Milky Way and M31 Merger Simulation
    
    # This Funcition takes in the arguements of: Galaxy Name (string)
    #                                            Array of times values you want
    #
    # This function requires known time values you wish to investigate
    #
    # This program does not save values to a file, but returns the list of files
    # you wish to examine without .txt extension
    
    
    # Define an empty list to store file names
    Target_list = []
    
    # Define file time incriments for each file.
    # Total of 801 snap shots from 0 to 800
    # Total time for simulation 11.4286 Gyr
    file_dt = 11.4286/801 
    
    # Array of Snapshot number (note these will be float need ints)
    Snaps = Time_array/file_dt
    
    # loop through Snaps and Store file name
    for ii in range(len(Snaps)):
        
        # First check: For snaps less than 100: want to return Galaxy_0(01-99)
        # Similar to OribitCOM assignment
        if int(Snaps[ii]) < 100:
            
            # String of number values with snap number wanted
            ilbl = "000" + str(int(Snaps[ii]))
            
            # Keep last 3 characters of string
            ilbl = ilbl[-3:]
            
            # Compose file name
            filename = Galaxy + "_" + ilbl
            
        else:
            # Compose file name
            filename = Galaxy + "_" + str(int(Snaps[ii]))
        # Append file to file name
        Target_list += [filename]
        
        
    
    return Target_list

    




