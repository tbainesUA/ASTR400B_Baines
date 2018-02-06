
# coding: utf-8

# In[3]:


# Program that will Calculate the Center of Mass Position and Velcoity


# In[271]:


"""Import python Modules"""
import numpy as np
import astropy.units as u
import astropy.table as tbl
"""Import my Modules"""
from Readfile import Read



# In[182]:


# Define Class Function
class CenterOfMass:

    def __init__(self, filename, ptype):
        # Function that takes inputs: object(self), the filename, and particle type
         
        # read in the file and particle type                                                                           
        self.time, self.total, self.data = Read(filename)

        #create an array to store indexes of particles of desired Ptype                                                
        self.index = np.where(self.data['type'] == ptype)

        # store the mass of only the particles of the given type                                
        self.m = self.data['m'][self.index]
        
        # store the positions of only the particles of the given type  
        self.x = self.data['x'][self.index]
        self.y = self.data['y'][self.index]
        self.z = self.data['z'][self.index]
        
        # store the velocity of only the particles of the given type  
        self.vx = self.data['vx'][self.index]
        self.vy = self.data['vy'][self.index]
        self.vz = self.data['vz'][self.index]
        ##### PLACE other particle properties here: x,y,z,vx,vy,vz ##### 
        """Completed this section"""

    
    def total_mass(self):
        #Function that takes in object and Returns the total mass
        #Note: you can add other keyword arguments into the function, but 'self' must be first                         
        return np.sum(self.m)*u.Msun*1e10
    """Completed this section"""

    ##### PLACE OTHER FUNCTIONS BELOW #####                                                                            

    def COMdefine(self, m, x, y, z):
        # Function that takes in inputs of  3D Position and Velocity and Mass of that "thing"
        # and Outputs the 3D Center of Mass Position and Velocity
        
        # Calculate X,Y,Z Center of Mass
        XCOM = np.sum(x*m)/np.sum(m)
        YCOM = np.sum(y*m)/np.sum(m)
        ZCOM = np.sum(z*m)/np.sum(m)
        
        # return calculated values
        return XCOM,YCOM,ZCOM
    """Completed this section"""
        
    def COM_P(self, delta):
        # Function that in object and delta (small change) and 
        # Returns 3D position converted values
        
        # Compute Center of Mass (COM) of 3D vector space
        XCOM, YCOM, ZCOM = self.COMdefine(self.m, self.x, self.y, self.z)
        
        # Magnitude of Position/velocity of particles
        RCOM = np.sqrt(XCOM*XCOM + YCOM*YCOM + ZCOM*ZCOM)
        
        # Change Particle positions to COM reference frame
        XREF_COM = self.x - XCOM
        YREF_COM = self.y - YCOM
        ZREF_COM = self.z - ZCOM
        
        # Magnitude of position/velocity to COM reference frame
        RREF_COM = np.sqrt(XREF_COM*XREF_COM + YREF_COM*YREF_COM + ZREF_COM*ZREF_COM)
        
        # Find the maximum value in the RREF_COM array
        RMAX = np.max(RREF_COM)

        #Difference
        difference = 1000
        
        """Start of While Loop"""
        while np.abs(difference) > delta:
            # Loop that runs over shrinking volumes in R direction to verify COM convergence.
            
            # Index for particles less than half RMAX
            index = np.where(RREF_COM <= RMAX/2.)[0]
            
           
            # Compute Center of Mass (COM) of 3D vector space at index values
            XCOM2, YCOM2, ZCOM2 = self.COMdefine(self.m[index], self.x[index], self.y[index], self.z[index])
            
            
            # Magnitude of Position/velocity of particles of smaller volume selection
            RCOM2 = np.sqrt(XCOM2*XCOM2 + YCOM2*YCOM2 + ZCOM2*ZCOM2)
            #print (RCOM2)
            
            # Change Particle positions to COM rotating reference frame coordinates in new volume
            XREF_COM2 = self.x[ugh] - XCOM2
            YREF_COM2 = self.y[ugh] - YCOM2
            ZREF_COM2 = self.z[ugh] - ZCOM2
            
            # Magnitude of position/velocity to COM rotating reference frame of new volume
            RREF_COM2 = np.sqrt(XREF_COM2*XREF_COM2 + YREF_COM2*YREF_COM2 + ZREF_COM2*ZREF_COM2)

            # Find the maximum value in the RREF_COM array of new volume
            RMAX2 = np.max(RREF_COM2)
            
            
            # Update guess
            difference = RCOM - RCOM2 
            RCOM = RCOM2
            RREF_COM = RREF_COM2
            RMAX = RMAX2
            
            # Return 3D position values
            """End of While Loop"""
        return XCOM2, YCOM2, ZCOM2
    """Completed this section"""

    
    def COM_V(self, XCOM, YCOM, ZCOM):
        # Function that takes in object(self), and 3D position arguments
        
        #Determine Position Coordinates in Rotating Frame
        XREF_COM = self.x - XCOM
        YREF_COM = self.y - YCOM
        ZREF_COM = self.z - ZCOM
        
        # Determine magnitude position of COM
        RREF_COM = np.sqrt(XREF_COM*XREF_COM + YREF_COM*YREF_COM + ZREF_COM*ZREF_COM)
        
        # Find indices for all particles within 15 Kpc from COM position
        mark = np.where(RREF_COM <= 15)[0]
        
        # Identify and store all velocity values that meet mark condition
        VX = self.vx[mark]
        VY = self.vy[mark]
        VZ = self.vz[mark]
        M = self.m[mark]
        
        # Compute 3D COM Velcocity components
        VXCOM, VYCOM, VZCOM = self.COMdefine(M, VX, VY, VZ)
        
        #return 3D Velocity components
        return VXCOM, VYCOM, VZCOM
    """Completed this section"""
    


# In[279]:


# Initial Class funcitons for each file of 3 Galaxies, MW, M31, M33
MWCOM = CenterOfMass("MW_000",2)
M33COM = CenterOfMass("M33_000",2)
M31COM = CenterOfMass("M31_000",2)

#Compute Center of Mass Position coordinates
MWP = np.round(MWCOM.COM_P(10),2)
M33P = np.round(M33COM.COM_P(10),2)
M31P = np.round(M31COM.COM_P(10),2)

#Compute Center of Mass Velocity Coordinates 
MWV = np.round(MWCOM.COM_V(MWP[0], MWP[1], MWP[2]),2)
M31V = np.round(M31COM.COM_V(M31P[0], M31P[1], M31P[2]),2)
M33V = np.round(M33COM.COM_V(M33P[0], M33P[1], M33P[2]),2)



"""Making Table of Results"""
# Define array of zero values to star computed position and velocity values
tab_results = np.zeros(0)

# List to be iterated through to attach to tab_results
Results = ['MW COM', MWP, MWV, 'M31 COM', M31P, M31V, 'M33 COM', M33P, M33V ]

# Iterate through list and append values to tab_results
for i in Results:
    tab_results = np.append(tab_results, i)

# Restructure array into 3 by 7
tab_results = np.reshape(tab_results, (3,7))

#Make and plot
t = tbl.Table(tab_results, names = ['Galaxy', 'X-component (kpc)', 'Y-Component (kpc)', 'Z-Component (kpc)','VX-component (km/s)', 'VY-Component (km/s)', 'VZ-Component (km/s)'])
t.show_in_notebook()



# In[306]:


"""Quesition 2"""
print ("Question 3 Results")
# Calculate the seperation between MW and  M31 position
SepP = np.round(np.sqrt(np.sum((MWP - M31P)**2)),3)*u.kpc
# Calculate the seperation between MW and  M31 velocity
SepV = np.round(np.sqrt(np.sum((MWV - M31V)**2)),3)*u.km/u.s
# Print results
print ("The magnitude of current seperation of position between MW and M31 is: %s\nThe magnitude of current seperation of velocity between MW and M31 is: %s" % (SepP,SepV))


"""Question 3"""
# Calculate the seperation between M33 and  M31 position
SepP = np.round(np.sqrt(np.sum((M33P - M31P)**2)),3)*u.kpc
# Calculate the seperation between M33 and  M31 velocity
SepV = np.round(np.sqrt(np.sum((M33V - M31V)**2)),3)*u.km/u.s
# Print results
print ("The magnitude of current seperation of position between M33 and M31 is: %s\nThe magnitude of current seperation of velocity between M33 and M31 is: %s" % (SepP,SepV))


# Question 4
# 
# The iterative process it important be would like to examine how the center of mass of the two galaxies interact with one another.
# 
# 
