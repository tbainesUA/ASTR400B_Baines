
# coding: utf-8

# In[4]:


#scrpit that finds and located specified particle info and prints physical quantities.
"""Import Packages"""
import numpy as np
import astropy.units as u
from Readfile import Read


# In[9]:


"""Particle Info Function"""
def ParticleInfo(particle_type, particle_num):
    
    #import data
    time, P_total, data = Read("MW_000")
    
    #located specific particle type
    index = np.where(data["type"]==particle_type)
    
    
    #specify quantities 
    """Mass [Particle Type] [ith particle select]"""
    m = data['m'][index][particle_num-1] # array containing all mass values of particle types and ith particle
    #m = m*1e10*u.M_sun #Calculate quantity with appropiate units
    
    """Positions [Particle Type] [ith particle select]"""
    x = data['x'][index][particle_num-1] # array containing all x coordinates [Particle Type]
    y = data['y'][index][particle_num-1] # array containing all y coordinates
    z = data['z'][index][particle_num-1] # array containing all z coordinates
    
    """Velocities [Particle Type] [ith particle select]"""
    vx = data['vx'][index][particle_num-1] # array containing all vx coordinates
    vy = data['vy'][index][particle_num-1] # array containing all vy coordinates
    vz = data['vz'][index][particle_num-1] # array containing all vz coordinates
    
    """3D Magnitude Calculations"""
    R = x*x+y*y+z*z #x^2+y^2+z^2
    V_R = vx*vx+vy*vy+vz*vz #vx^2+vy^2+vz^2
    #calculate magnitues
    r = np.round(np.sqrt(R),3)*u.kpc #round 3D position magnitude to 3 decimal places
    v_r = np.round(np.sqrt(V_R),3)*u.km/u.s #round 3D velocity magnitude
    
    """Print Results"""
    #print ("Position: %s" % r)
    #print ("Position: %s" % np.round(r.to(u.lyr),3))
    #print ("Velocity: %s" % v_r)
    #print ("Mass: %s" % (m*1e10*u.M_sun))


# In[10]:


ParticleInfo(2,100)

