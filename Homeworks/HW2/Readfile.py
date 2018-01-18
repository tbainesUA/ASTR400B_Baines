
# coding: utf-8

# In[1]:


"""Import Packages"""
import numpy as np
import astropy.units as u


# In[3]:


"""Loadfile Function"""
def Read(filename):
    file = open(filename, 'r')
    line1 = file.readline()
    label, value = line1.split
    time = float(value)*10.0*u.Myr
    
    line2 = file.readline()
    label, value = line2.split
    particles = float(values)
    
    file.close()
    
    data = np.genfromtxt(filename, dtype = None, names = True, skip_header = 3)
    return time, particles, data

