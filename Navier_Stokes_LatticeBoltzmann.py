#!/usr/bin/env python
# coding: utf-8

# In[2]:


#%matplotlib notebook
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

def main():
    """ Lattice Boltzmann Simulation """
    
    # Simulation parameters
    Nx                     = 400    # resolution x-dir
    Ny                     = 100    # resolution y-dir
    rho0                   = 90    # average density
    tau                    = 0.6    # collision timescale
    Nt                     = 700   # number of timesteps

    # Lattice speeds / weights
    NL = 9
    idxs = np.arange(NL)
    cxs = np.array([0, 0, 1, 1, 1, 0,-1,-1,-1])
    cys = np.array([0, 1, 1, 0,-1,-1,-1, 0, 1])
    weights = np.array([4/9,1/9,1/36,1/9,1/36,1/9,1/36,1/9,1/36]) 
    """
    8,1,2
    7,0,3
    6,5,4
    """
    
    # Initial Conditions
    F = np.ones((Ny,Nx,NL)) 
    F[1:-2,int(Nx/4),3] += 1
    
    rho = np.sum(F,2)
    ux  = np.sum(F*cxs,2) / rho
    uy  = np.sum(F*cys,2) / rho
    
    for i in idxs:
        F[:,:,i] *= rho0 / rho

    # boundary conditions
    X, Y = np.meshgrid(range(Nx), range(Ny))
    #Circular obstacle
    #BC = (X - 3*Nx/4)**2 + (Y - Ny/2)**2 < (Ny/4)**2 
    #Slit
    BC = np.full((Ny, Nx), False, dtype=bool)
    #BC[::2,3*Nx//4] = True
    BC[0:40,3*Nx//4] = True
    BC[60:-1,3*Nx//4] = True
    #Walls
    BC[:,[0,-1]] = True
    BC[[0,-1],:] = True
    
        
    # Simulation Main Loop
    for it in range(Nt):    
        # Calculate fluid variables
        rho = np.sum(F,2)
        ux  = np.sum(F*cxs,2) / rho
        uy  = np.sum(F*cys,2) / rho
        
        # Relax towards Maxwell-Boltzmann probability distribution function
        Feq = np.zeros(F.shape)
        for i, cx, cy, w in zip(idxs, cxs, cys, weights):
            Feq[:,:,i] = rho * w * ( 1 + 3*(cx*ux+cy*uy)  + (9/2)*(cx*ux+cy*uy)**2 - 3*(ux**2+uy**2)/2 )
        
        F += -(1.0/tau) * (F - Feq)
        
        # Calculate reflective boundaries
        bndryF = F[BC,:]
        bndryF = bndryF[:,[0,5,6,7,8,1,2,3,4]] 
        
        # Apply boundary 
        F[BC,:] = bndryF 

        # Drift
        for i, cx, cy in zip(idxs, cxs, cys):
            F[:,:,i] = np.roll(F[:,:,i], cx, axis=1)
            F[:,:,i] = np.roll(F[:,:,i], cy, axis=0)
           
        # plot
        
        if ((it % 100) == 0): 
            Density = np.ma.array(rho, mask=BC)
            
            fig = plt.figure(figsize=(8,4))
            plt.clf()
            plt.imshow(Density, cmap='cool')
            plt.imshow(~BC, cmap='gray', alpha=0.3)
            plt.show()
        
         
    #Plot 3D
    Density = np.ma.array(rho, mask=BC)
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"},figsize=(8,4))       
    surface = ax.plot_surface(X, Y, Density, cmap=cm.coolwarm,linewidth=0, antialiased=False)
    #boundary = ax.plot_surface(X, Y,~BC, cmap='gray', alpha=0.3)
    plt.show()
    
    return 0

if __name__== "__main__":
  main()


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




