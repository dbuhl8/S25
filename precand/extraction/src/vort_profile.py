import dbuhlMod as db
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors
import matplotlib.animation as animation
from matplotlib import gridspec
from matplotlib.offsetbox import AnchoredText
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.widgets import Slider
from netCDF4 import MFDataset

#plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["mathtext.fontset"] = 'stix'
plt.rcParams["font.size"] = 16

cm_data = [[0.2081, 0.1663, 0.5292], [0.2116238095, 0.1897809524, 0.5776761905],
     [0.212252381, 0.2137714286, 0.6269714286], [0.2081, 0.2386, 0.6770857143],
     [0.1959047619, 0.2644571429, 0.7279], [0.1707285714, 0.2919380952,
      0.779247619], [0.1252714286, 0.3242428571, 0.8302714286],
     [0.0591333333, 0.3598333333, 0.8683333333], [0.0116952381, 0.3875095238,
      0.8819571429], [0.0059571429, 0.4086142857, 0.8828428571],
     [0.0165142857, 0.4266, 0.8786333333], [0.032852381, 0.4430428571,
      0.8719571429], [0.0498142857, 0.4585714286, 0.8640571429],
     [0.0629333333, 0.4736904762, 0.8554380952], [0.0722666667, 0.4886666667,
      0.8467], [0.0779428571, 0.5039857143, 0.8383714286],
     [0.079347619, 0.5200238095, 0.8311809524], [0.0749428571, 0.5375428571,
      0.8262714286], [0.0640571429, 0.5569857143, 0.8239571429],
     [0.0487714286, 0.5772238095, 0.8228285714], [0.0343428571, 0.5965809524,
      0.819852381], [0.0265, 0.6137, 0.8135], [0.0238904762, 0.6286619048,
      0.8037619048], [0.0230904762, 0.6417857143, 0.7912666667],
     [0.0227714286, 0.6534857143, 0.7767571429], [0.0266619048, 0.6641952381,
      0.7607190476], [0.0383714286, 0.6742714286, 0.743552381],
     [0.0589714286, 0.6837571429, 0.7253857143],
     [0.0843, 0.6928333333, 0.7061666667], [0.1132952381, 0.7015, 0.6858571429],
     [0.1452714286, 0.7097571429, 0.6646285714], [0.1801333333, 0.7176571429,
      0.6424333333], [0.2178285714, 0.7250428571, 0.6192619048],
     [0.2586428571, 0.7317142857, 0.5954285714], [0.3021714286, 0.7376047619,
      0.5711857143], [0.3481666667, 0.7424333333, 0.5472666667],
     [0.3952571429, 0.7459, 0.5244428571], [0.4420095238, 0.7480809524,
      0.5033142857], [0.4871238095, 0.7490619048, 0.4839761905],
     [0.5300285714, 0.7491142857, 0.4661142857], [0.5708571429, 0.7485190476,
      0.4493904762], [0.609852381, 0.7473142857, 0.4336857143],
     [0.6473, 0.7456, 0.4188], [0.6834190476, 0.7434761905, 0.4044333333],
     [0.7184095238, 0.7411333333, 0.3904761905],
     [0.7524857143, 0.7384, 0.3768142857], [0.7858428571, 0.7355666667,
      0.3632714286], [0.8185047619, 0.7327333333, 0.3497904762],
     [0.8506571429, 0.7299, 0.3360285714], [0.8824333333, 0.7274333333, 0.3217],
     [0.9139333333, 0.7257857143, 0.3062761905], [0.9449571429, 0.7261142857,
      0.2886428571], [0.9738952381, 0.7313952381, 0.266647619],
     [0.9937714286, 0.7454571429, 0.240347619], [0.9990428571, 0.7653142857,
      0.2164142857], [0.9955333333, 0.7860571429, 0.196652381],
     [0.988, 0.8066, 0.1793666667], [0.9788571429, 0.8271428571, 0.1633142857],
     [0.9697, 0.8481380952, 0.147452381], [0.9625857143, 0.8705142857, 0.1309],
     [0.9588714286, 0.8949, 0.1132428571], [0.9598238095, 0.9218333333,
      0.0948380952], [0.9661, 0.9514428571, 0.0755333333],
     [0.9763, 0.9831, 0.0538]]
 
parula_map = LinearSegmentedColormap.from_list('parula', cm_data)

def make_norm(maxnum):
    return colors.Normalize(vmin=-maxnum,vmax=maxnum)

def sclmap(cmap, maxnum):
    return mpl.cm.ScalarMappable(cmap=cmap,norm=make_norm(maxnum))

fs = 16
dfs = 4
ptstp = -1

font_fam = 'TimesRoman'

# opening slice data files
fnxy = 'XYSLICE*.cdf'
#fnxz = 'XZSLICE3.cdf'
#fnyz = 'YZSLICE3.cdf'
cdf_filexy = MFDataset(fnxy)
#cdf_filexz = MFDataset(fnxz)
#cdf_fileyz = MFDataset(fnyz)
dtype = np.float32
cmap = 'RdYlBu_r'
num_contours = 100 # enforces smoothness

#print(cdf_filexy.variables)

#slice files need to manually assign Ro
Ro = 0.2

# obtaining discretization data
x = np.array(cdf_filexy.variables["x"])
y = np.array(cdf_filexy.variables["y"])
#z = np.array(cdf_filexz.variables["z"])
t = np.array(cdf_filexy.variables["t"][:])
Nx = len(x) # using MFDataset (these might be extra long)
Ny = len(y)
#Nz = len(z)
Nt = len(t)
print('Number of Timesteps: ', Nt)
print('Number of Gridpoints in x: ', Nx)
print('Number of Gridpoints in y: ', Ny)
gx = 4*np.pi
gy = 4*np.pi
#gz = np.pi
dx = gx/Nx
dy = gy/Ny
#dz = gz/Nz
dt = t[1]-t[0]



#these arrays are indexed by [t,z,y,x]
ux_xy =  np.array(cdf_filexy.variables["ux"][:])
uy_xy =  np.array(cdf_filexy.variables["uy"][:])
wz = db.FD6X_xyslice(uy_xy,Nx,dx) - db.FD6Y_xyslice(ux_xy,Ny,dy)
wzmax = np.max([np.abs(ux_xy).max(), np.abs(uy_xy).max()])
r2 = np.zeros((Ny,Nx))
uz =  np.array(cdf_filexy.variables["uz"][:])

del ux_xy
del uy_xy

max_dist = Nx/2

# find vortex core
#for k in range(1):
    # finds the max of the array
k = -1
wzmax = np.max(np.abs(wz[k,:,:]))
ind = np.unravel_index(np.argmax(wz[k,:,:],axis=None), wz[k,:,:].shape)
print(np.max(wz[k,:,:]))

fig, ax = plt.subplots(2,2,figsize=(16,9))
wz_plot = ax[0,0].imshow(wz[k,:,:],norm=colors.Normalize(vmin=-wzmax,vmax=wzmax),
    cmap = 'RdYlBu_r', origin='lower')
cb1 = fig.colorbar(wz_plot, ax=ax[0,0],fraction=0.02, pad=0.1)
cb1.ax.tick_params(labelsize=fs-dfs)
ax[0,0].set_xticks([])
ax[0,0].set_yticks([])


#print(ind)
indx = ind[1]
indy = ind[0]
for j in range(Ny):
    for i in range(Nx):
        xd =  abs(i-indx)
        yd  =  abs(j-indy)
        if (xd == max_dist):
            xd = max_dist
        elif (xd > max_dist):
            xd = 2*max_dist - xd
        if (yd == max_dist):
            yd = max_dist
        elif (yd > max_dist):
            yd = 2*max_dist - yd
        #if (abs(i-indy) == max_dist):
            #yd = max_dist
        r2[j,i] = np.sqrt((xd*dx)**2 + (yd*dy)**2)
r2[indy,indx] = 0

# this creates a new 1d array for r which is sorted from 0 to r_max and has no
# repeat values of r
r = np.unique(r2.flatten())
# defining new arrays for desired quantities
ombar = np.zeros_like(r)
ang_mom = np.zeros_like(r)
nr = len(r)

r_plot = ax[0,1].imshow(r2,norm=colors.Normalize(vmin=r[0],vmax=r[-1]),
    cmap = 'hsv', origin='lower')
cb2 = fig.colorbar(r_plot, ax=ax[0,1],fraction=0.02, pad=0.1)
cb2.ax.tick_params(labelsize=fs-dfs)
ax[0,1].set_xticks([])
ax[0,1].set_yticks([])

print('Number of radii: ', nr)
for i in range(nr):
    ind = np.where(r2 == r[i])
    num_pt = len(ind[0])
    #if (i%10 == 0) and (i  < 101):
        #print('Number of points in avg: ',num_pt)
        #print('Indices Array : ',ind)
    for j in range(num_pt):
        ombar[i] += wz[k,ind[0][j],ind[1][j]]
    ombar[i] = ombar[i]/num_pt
for i in range(nr-4):
    ang_mom[i+2] = (-(r[i+4]**2)*ombar[i+4] + 8*(r[i+3]**2)*ombar[i+3] \
        -8*(r[i+1]**2)*ombar[i+1]+(r[i]**2)*ombar[i])/\
        (12*(r[i+4]-r[i])*r[i+1]) + 1./Ro

ang_mom[0] = ang_mom[1]
ang_mom[nr-1] = ang_mom[nr-2]
ax[1,0].plot(r,r**2*ombar, 'rx')
ax[1,1].plot(r,ombar,'bo')
fig.tight_layout()
plt.savefig('vortex_profile.pdf')
#plt.show()
# Finding max values for the colorbar
#uzmax = (np.abs(uz).max())/2




# plotting with matplotlib

# makes a box plot for each subplot
#fig, ax = plt.subplots(2, 1,subplot_kw=dict(projection='3d'))
#fig = plt.figure(figsize=(6,3))

#gs = gridspec.GridSpec(1, 2, width_ratios=[1,1]) 
#ax1 = fig.add_subplot(gs[0])
#ax2 = fig.add_subplot(gs[1])
#ax2 = fig.add_subplot(gs[1],projection='3d')
#ax2.computed_zorder = False


#wz plot
#wz_top = ax1.imshow(big_wz[:,:], norm=colors.Normalize(vmin=-wzmax,vmax=wzmax),
    #cmap = 'RdYlBu_r', origin='lower')
#cb1 = fig.colorbar(wz_top, ax=ax1,fraction=0.02, pad=0.1, )
#cb1.ax.tick_params(labelsize=fs-dfs)
#ax1.add_artist(AnchoredText(r'$\omega_z$',loc='upper center',prop=dict(size=fs,
    #weight='bold'),frameon=False))
#ax1.set_title(r'$\omega_z$')
#ax1.set_ylabel(r'$y$',rotation=0,labelpad=10)
#ax1.set_xlabel(r'$x$')

#uz_top = ax2.imshow(big_uz[:,:], norm=colors.Normalize(vmin=-uzmax,vmax=uzmax),
    #cmap = 'seismic', origin='lower')
#cb2 = fig.colorbar(uz_top, ax=ax2,fraction=0.02, pad=0.1, )
#ax2.set_title(r'$u_z$')
#ax2.set_ylabel(r'$y$',rotation=0,labelpad=10)
#ax2.set_xlabel(r'$x$')


# Additional Plot Formating
# makes room for the slider
#fig.subplots_adjust(0, .05, .90, .90, .05,.05)

#ax1.set_xticks([])
#ax1.set_yticks([])
#ax2.set_xticks([])
#ax2.set_yticks([])

# horizaontally oriented slider
#taxis = plt.axes([0.15, 0.02, 0.7, 0.03], facecolor='blue')
#staxis = Slider(taxis, 'Time', t[0], t[-1], valinit=t[0], valstep=dt)

#fig.tight_layout()
# movie down vertical extent of domain
def update_frame(frame):
    big_wz[0:Ny,0:Nx] = wz[frame,:,:]
    big_wz[Ny:2*Ny,0:Nx] = wz[frame,:,:]
    big_wz[0:Ny,Nx:2*Nx] = wz[frame,:,:]
    big_wz[Ny:2*Ny,Nx:2*Nx] = wz[frame,:,:]
    big_uz[0:Ny,0:Nx] = uz[frame,:,:]
    big_uz[Ny:2*Ny,0:Nx] = uz[frame,:,:]
    big_uz[0:Ny,Nx:2*Nx] = uz[frame,:,:]
    big_uz[Ny:2*Ny,Nx:2*Nx] = uz[frame,:,:]
    wz_top.set_array(big_wz)
    uz_top.set_array(big_uz)
    staxis.set_val(t[frame]) 
    plt.savefig('xyslice_frame{:0{}d}.png'.format(frame,5),dpi=800)
    print('Done with frame: ', frame)
    return (uz_top, wz_top)

#for i in range(Nt):
    #plots = update_frame(i)

#ani = animation.FuncAnimation(fig=fig,
    #func=update_frame,frames=Nt,interval=100,blit=True)
#ani.save('Om1B100Re1000Pe100_xyslice.gif')
#plt.show()

