import numpy as np
import matplotlib.pyplot as plt

import mol_dyn_vis as mol_dyn_vis

T = np.linspace(-2,100,100)      # Temperature [degree celsius]
p = np.linspace(0,60000,101)     # Pressure [dbar]
S = np.linspace(0,45,46)         # Salinity [psu]

mu = np.zeros((len(T),len(p),len(S)))

# Test values provided by 
# Dietrich, Kalle, Kraus & Siedler 1975: Allgemeine Meereskunde. 
#                                        Eine Einfuehrung in die Ozeanographie.
#                                        Gebrueder Borntraeger Berlin.
mu_check = mol_dyn_vis.mol_dyn_vis(3907,2.42,34.917) # returns units: kg/(m s)
# should result in 1.6934e-2 g/cm/s:
print('Check value 0.016934 = '+str(mu_check*10.))


# Since the formulation for dynamic molecular viscosity goes negative for high temperatures, 
# we introduce a new minimum value 
minval = mol_dyn_vis.mol_dyn_vis(10,45,0)
print('New minimum value: '+str(minval))

for i,temp in enumerate(T):
    for j,press in enumerate(p):
        for k,sal in enumerate(S): 
            mu[i,j,k] = mol_dyn_vis.mol_dyn_vis(press,temp,sal)

fig=plt.figure(figsize=(16,8))
fig.subplots_adjust(hspace=0.4, wspace=0.3)
plt.subplot(231)
plt.pcolormesh(T,S,mu[:,0,:].squeeze().transpose(),rasterized=True)
plt.xlabel('T')
plt.ylabel('S')
cb=plt.colorbar()
cb.ax.set_title(r'kg$\,$m$^{-1}\,$s$^{-1}$')
cb.ax.axhline(y=0., color='red')
cb.ax.axhline(y=minval, color='orchid')
plt.contour(T,S,mu[:,0,:].squeeze().transpose(), levels=[0,minval],colors=['red','orchid'])
plt.grid(ls='--',color='lightgray',alpha=0.5)
plt.title(r'$\mu$ at '+str(p[0])+' dbar')

idxp=50
plt.subplot(232)
plt.pcolormesh(T,S,mu[:,idxp,:].squeeze().transpose(),rasterized=True)
plt.xlabel('T')
plt.ylabel('S')
cb=plt.colorbar()
cb.ax.axhline(y=0., color='red')
cb.ax.axhline(y=minval, color='orchid')
cb.ax.set_title(r'kg$\,$m$^{-1}\,$s$^{-1}$')
plt.contour(T,S,mu[:,idxp,:].squeeze().transpose(), levels=[0,minval],colors=['red','orchid'])
plt.grid(ls='--',color='lightgray',alpha=0.5)
plt.title(r'$\mu$ at '+str(p[idxp])+' dbar')

idxp=100
plt.subplot(233)
plt.pcolormesh(T,S,mu[:,idxp,:].squeeze().transpose(),rasterized=True)
plt.xlabel('T')
plt.ylabel('S')
cb=plt.colorbar()
cb.ax.axhline(y=0., color='red')
cb.ax.axhline(y=minval, color='orchid')
cb.ax.set_title(r'kg$\,$m$^{-1}\,$s$^{-1}$')
plt.contour(T,S,mu[:,idxp,:].squeeze().transpose(), levels=[0,minval],colors=['red','orchid'])
plt.grid(ls='--',color='lightgray',alpha=0.5)
plt.title(r'$\mu$ at '+str(p[idxp])+' dbar')


idxS=0
plt.subplot(234)
plt.pcolormesh(T,p,mu[:,:,idxS].squeeze().transpose(),rasterized=True)
plt.gca().invert_yaxis()
plt.xlabel('T')
plt.ylabel('p')
cb=plt.colorbar()
cb.ax.axhline(y=0., color='red')
cb.ax.axhline(y=minval, color='orchid')
cb.ax.set_title(r'kg$\,$m$^{-1}\,$s$^{-1}$')
plt.contour(T,p,mu[:,:,idxS].squeeze().transpose(), levels=[0,minval],colors=['red','orchid'])
plt.grid(ls='--',color='lightgray',alpha=0.5)
plt.title(r'$\mu$ at salinity '+str(S[idxS])+' psu')

idxS=34
plt.subplot(235)
plt.pcolormesh(T,p,mu[:,:,idxS].squeeze().transpose(),rasterized=True)
plt.gca().invert_yaxis()
plt.xlabel('T')
plt.ylabel('p')
cb=plt.colorbar()
cb.ax.axhline(y=0., color='red')
cb.ax.axhline(y=minval, color='orchid')
cb.ax.set_title(r'kg$\,$m$^{-1}\,$s$^{-1}$')
plt.contour(T,p,mu[:,:,idxS].squeeze().transpose(), levels=[0,minval],colors=['red','orchid'])
plt.grid(ls='--',color='lightgray',alpha=0.5)
plt.title(r'$\mu$ at salinity '+str(S[idxS])+' psu')

idxS=45
plt.subplot(236)
plt.pcolormesh(T,p,mu[:,:,idxS].squeeze().transpose(),rasterized=True)
plt.gca().invert_yaxis()
plt.xlabel('T')
plt.ylabel('p')
cb=plt.colorbar()
cb.ax.axhline(y=0., color='red')
cb.ax.axhline(y=minval, color='orchid')
cb.ax.set_title(r'kg$\,$m$^{-1}\,$s$^{-1}$')
plt.contour(T,p,mu[:,:,idxS].squeeze().transpose(), levels=[0,minval],colors=['red','orchid'])
plt.grid(ls='--',color='lightgray',alpha=0.5)
plt.title(r'$\mu$ at salinity '+str(S[idxS])+' psu')

plt.suptitle('Testing dynamic molecular viscosity '+r'$\mu$:' +'   red - 0 line, orchid - new minimum value',fontweight='bold')
fig.savefig('test_mol_dyn_vis.pdf')
