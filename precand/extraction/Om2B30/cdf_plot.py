import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

# this file is meant to work for a simdat.cdf file. 
# different configurations exsit for slice files which are 2D in nature. 

def FDX(field,dx):
    # takes in a field dataset with 4 dimensions (3 spatial, 1 time)
    # returns dataset (numpy) of same shape but computing the 1st derivative
    # with respect to x_comp
    dims = field.dims
    dx_field = field
    # foward / backwards finite difference formula (boundary)
    dx_field[:,0,:,:] = (field[:,1,:,:]-field[:,0,:,:])/dx
    dx_field[:,dims,:,:] = (field[:,1,:,:]-field[:,0,:,:])/dx
    # centered finite difference formula (inside)
    for i in range(dims[1]-2):
        dx_field[:,i+1,:,:] = (field[:,i+2,:,:]-field[:,i,:,:])/(2*dx)
    return dx_field

def FDY(field,dy):
    # takes in a field dataset with 4 dimensions (3 spatial, 1 time)
    # returns dataset (numpy) of same shape but computing the 1st derivative
    # with respect to x_comp
    dims = field.dims
    dy_field = field
    # foward / backwards finite difference formula (boundary)
    dy_field[:,0,:,:] = (field[:,:,1,:]-field[:,0,0,:])/dy
    dy_field[:,dims,:,:] = (field[:,:,1,:]-field[:,:,0,:])/dy
    # centered finite difference formula (inside)
    for i in range(dims[2]-2):
        dy_field[:,:,i+1,:] = (field[:,:,i+2,:]-field[:,:,i,:])/(2*dy)
    return dy_field

cdf_file = xr.open_dataset('simdat4.cdf')

print(cdf_file.variables)

x = cdf_file.x.as_numpy()
y = cdf_file.y.as_numpy()
z = cdf_file.z.as_numpy()

gx = cdf_file.Gammax
gy = cdf_file.Gammay
gz = cdf_file.Gammaz

dX = np.zeros(3)
dX[0] = len(x)/gx
dX[1] = len(y)/gy
dX[2] = len(z)/gz

# these matrices are oriented with index [t, z, y, x] (dont ask me why)
ux = cdf_file.ux.as_numpy()
uy = cdf_file.uy.as_numpy()
uz = cdf_file.uz.as_numpy()

wz = FDX(uy, dX[0]) - FDY(ux, dX[1])

fig, ax = plt.subplots()

# make a movie over timesteps
ax.pcolor(XX, ZZ, wz[-1,:,0,:],
