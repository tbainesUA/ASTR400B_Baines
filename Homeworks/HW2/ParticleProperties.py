"""Import Packages"""
import numpy as np
import astropy.units as u
from ReadFile import Read

"""Particle Info Function"""
def ParticleInfo(particle_type, particle_num):
    
    time, P_total, data = Read("MW_000.txt")

    index = np.where(data['x']>1)
"""Positions"""
    x = data['x'][index]
    y = data['y'][index]
    z = data['z'][index]
"""Velocities"""
    vx = data['vx'][index]
    vy = data['vy'][index]
    vz = data['vz'][index]
"""3D translations"""
    R = x*x+y*y+z*z
    VR = x*x+y*y+z*z
    r = np.sqrt(R)
    v_r = np.sqrt(V_R)
