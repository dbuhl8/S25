import numpy as np
import dbuhlMod as db
from netCDF4 import MFDataset

# commented out because this breaks expanse plotting
# plt.rcParams['text.usetex'] = True

dtype = np.float32
ptstp = -1  # this choses which timestep to plot (-1 is the last one)
fs = 18
dfs = 6
font = {'family': 'Times New Roman',
                        'size'   : fs}

files_30 = ['Om0.5B30Re600Pe60/OUT.dat',\
    'Om1B30Re600Pe60/OUT.dat','Om2B30Re600Pe60/OUT.dat', \
    'Om3B30Re600Pe60/OUT.dat','Om4B30Re600Pe60/OUT.dat',\
    'Om5B30Re600Pe60/OUT.dat','Om8B30Re600Pe60/OUT.dat','Om10B30Re600Pe60/OUT.dat']

files_100 = ['Om0.5B100Re600Pe60/OUT.dat',\
    'Om1B100Re600Pe60/OUT.dat','Om2B100Re600Pe60/OUT.dat', \
    'Om3B100Re600Pe60/OUT.dat','Om4B100Re600Pe60/OUT.dat',\
    'Om8B100Re600Pe60/OUT.dat']
cor_idx = [0, 1, 2, 3, 4, 6]

simdat_files_30 = ['Om0.5B30Re600Pe60/simdat*.cdf',\
    'Om1B30Re600Pe60/simdat*.cdf','Om2B30Re600Pe60/simdat*.cdf',\
    'Om3B30Re600Pe60/simdat*.cdf','Om4B30Re600Pe60/simdat*.cdf',\
    'Om5B30Re600Pe60/simdat*.cdf', 'Om8B30Re600Pe60/simdat*.cdf',\
    'Om10B30Re600Pe60/simdat*.cdf']
sim_30_idx = [0, 1, 2, 3, 4, 5, 6, 7]

simdat_files_100 = ['Om0.5B100Re600Pe60/simdat*.cdf',\
    'Om1B100Re600Pe60/simdat*.cdf','Om2B100Re600Pe60/simdat*.cdf', \
    'Om3B100Re600Pe60/simdat*.cdf','Om4B100Re600Pe60/simdat*.cdf',\
    'Om8B100Re600Pe60/simdat*.cdf']
sim_100_idx = [0,1,2,3,4,5]

invRo_30 = [0.5, 1, 2, 3, 4, 5, 8, 10]
invRo_100 = [0.5, 1, 2, 3, 4, 8]
invFr = np.array([np.sqrt(30), 10.])
Fr = np.array([1./f for f in invFr])

labels_100 = [r'$\Omega=0.5$',r'$\Omega=1$',r'$\Omega=2$',r'$\Omega=3$',r'$\Omega=4$',r'$\Omega=8$']
labels_30 = [r'$\Omega=0.5$',r'$\Omega=1$',r'$\Omega=2$',r'$\Omega=3$',r'$\Omega=4$',r'$\Omega=5$',r'$\Omega=8$',r'$\Omega=10$']
cor_label_30 = [1,2,3,4,5,7]

num_files = [len(files_30), len(files_100)]
sim_num_files = [len(simdat_files_30), len(simdat_files_100)]
# this is redundant now, algorithm has changed so there is not a fixed number of
# samples. Keep this set to one so that the numpy array are initialized to be 2D
# (critical for the hstack function called later in simdat loop)
num_samples = 1
cvar = 0.2

def add_column(mat,l):
    return np.hstack((mat, np.zeros((l,1))))

# these files will have an arbitrary number of timesteps for different files,
# arrays will be made 2D such that a computation can be done for each timestep
# for initial testing this will be keep as a 1d array to tst MFDataset
vol_frac_30 = np.zeros([sim_num_files[0],num_samples])
horiz_urms_30 = np.zeros_like(vol_frac_30)
uhrms_in_turb_30 = np.zeros_like(vol_frac_30)
t_30 = np.zeros_like(vol_frac_30)
wt_30 = np.zeros_like(vol_frac_30)
wt_in_turb_30 = np.zeros_like(vol_frac_30)
wrms_30 = np.zeros_like(vol_frac_30)
wrms_in_turb_30 = np.zeros_like(vol_frac_30)
tdisp_30 = np.zeros_like(vol_frac_30)
tdisp_in_turb_30 = np.zeros_like(vol_frac_30)
tdisp_in_lam_30 = np.zeros_like(vol_frac_30)
enstr_30 = np.zeros_like(vol_frac_30)
enstr_in_turb_30 = np.zeros_like(vol_frac_30)
enstr_in_lam_30 = np.zeros_like(vol_frac_30)

vol_frac_100 = np.zeros([sim_num_files[1], num_samples])
t_100 = np.zeros_like(vol_frac_100)
wt_100 = np.zeros_like(vol_frac_100)
wrms_100 = np.zeros_like(vol_frac_100)
enstr_100 = np.zeros_like(vol_frac_100)
tdisp_100 = np.zeros_like(vol_frac_100)
horiz_urms_100 = np.zeros_like(vol_frac_100)
wt_in_turb_100 = np.zeros_like(vol_frac_100)
wrms_in_turb_100 = np.zeros_like(vol_frac_100)
tdisp_in_turb_100 = np.zeros_like(vol_frac_100)
tdisp_in_lam_100 = np.zeros_like(vol_frac_100)
uhrms_in_turb_100 = np.zeros_like(vol_frac_100)
enstr_in_turb_100 = np.zeros_like(vol_frac_100)
enstr_in_lam_100 = np.zeros_like(vol_frac_100)

# can change these definitions to have simple 2d structure and then append
eInvRo_30 = np.zeros([num_files[0], num_samples])
eInvRo_err_30 = np.zeros_like(eInvRo_30)
avg_wT_30 = np.zeros_like(eInvRo_30)
avg_wT_err_30 = np.zeros_like(eInvRo_30)

eInvRo_100 = np.zeros([num_files[1], num_samples])
eInvRo_err_100 = np.zeros_like(eInvRo_100)
avg_wT_100 = np.zeros_like(eInvRo_100)
avg_wT_err_100 = np.zeros_like(eInvRo_100)

for i, fn in enumerate(simdat_files_30):
    cdf_file = MFDataset(fn)

    #obtaining discretization data
    x = np.array(cdf_file.variables["x"])
    y = np.array(cdf_file.variables["y"])
    z = np.array(cdf_file.variables["z"])
    t = np.array(cdf_file.variables["t"][:])
    Nx = len(x) # using MFDataset (these might be extra long)
    Ny = len(y)
    Nz = len(z)
    Nt = len(t)
    gx = cdf_file.variables["Gammax"][0]
    gy = cdf_file.variables["Gammay"][0]
    gz = cdf_file.variables["Gammaz"][0]
    dx = gx/Nx
    dy = gy/Ny
    dz = gz/Nz

    #these arrays are indexed by [t,z,y,x]
    ux = np.array(cdf_file.variables["ux"][:])
    uy = np.array(cdf_file.variables["uy"][:])
    uz = np.array(cdf_file.variables["uz"][:])
    temp = np.array(cdf_file.variables["Temp"][:])
    wz =  db.FD6X(uy, Nx, dx) - db.FD6Y(ux, Ny, dy)
    enstr = (db.FD6Y(uz, Ny, dy) - db.FD4Z(uy, Nz, dz))**2 +\
            (db.FD4Z(ux, Nz, dz) - db.FD6X(uz, Nx, dx))**2 + wz**2
    tdisp = np.sqrt(db.FD6X(temp, Nx, dx)**2 + db.FD6Y(temp, Ny, dy)**2 +\
            db.FD4Z(temp, Nz, dz)**2)

    #compute volume fraction where wz >= .1 wzmax
    for j, tval in enumerate(t):
        # if array has j-1 columns, create new column
        if len(horiz_urms_30[0,:])-1 < j:
            # adds an extra column of zeros to the arrays
            l = sim_num_files[0]
            vol_frac_30 = add_column(vol_frac_30,l)
            horiz_urms_30 = add_column(horiz_urms_30,l)
            uhrms_in_turb_30 = add_column(uhrms_in_turb_30,l)
            wrms_30 = add_column(wrms_30,l)
            wrms_in_turb_30 = add_column(wrms_in_turb_30,l)
            wt_30 = add_column(wt_30,l)
            t_30 = add_column(t_30,l)
            wt_in_turb_30 = add_column(wt_in_turb_30,l)
            enstr_30 = add_column(enstr_30,l)
            enstr_in_turb_30 = add_column(enstr_in_turb_30,l)
            enstr_in_lam_30 = add_column(enstr_in_lam_30,l)
            tdisp_30 = add_column(tdisp_30,l)
            tdisp_in_turb_30 = add_column(tdisp_in_turb_30,l)
            tdisp_in_lam_30 = add_column(tdisp_in_lam_30,l)

        # finding vertical pencils for turbulent and laminar regions
        idx = np.where(Fr[0]*np.abs(np.sum(wz[j,:,:,:],axis=1)/Nz+invRo_30[i]) <= 1)
        lam = np.where(Fr[0]*np.abs(np.sum(wz[j,:,:,:],axis=1)/Nz+invRo_30[i]) > 1)

        # volume average over whole domain
        t_30[i,j] = tval
        horiz_urms_30[i,j] = np.sqrt(np.sum(ux[j,:,:,:]**2 +uy[j,:,:,:]**2)/(Nx*Ny*Nz))
        vol_frac_30[i,j] = len(idx[0])/(Nx*Ny)
        wrms_30[i,j] = np.sqrt(np.sum(uz[j,:,:,:]**2)/(Nz*Ny*Nx))
        wt_30[i,j] = np.sum(uz[j,:,:,:]*temp[j,:,:,:])/(Nz*Ny*Nx)
        enstr_30[i,j] = np.sum(enstr[j,:,:,:])/(Nz*Ny*Nx)
        tdisp_30[i,j] = (np.sum(tdisp[j,:,:,:])/(Nz*Ny*Nx))**2
        
        # compute rms values in turbulent and laminar regions
        wt_in_turb_30[i,j] = np.sum(np.sum(uz[j,:,idx[0],idx[1]]*\
            temp[j,:,idx[0],idx[1]],axis=1)/Nz)/len(idx[0])
        uhrms_in_turb_30[i,j] = np.sqrt(np.sum(\
            np.sum(ux[j,:,idx[0],idx[1]]**2,axis=1)/Nz +
            np.sum(uy[j,:,idx[0],idx[1]]**2,axis=1)/Nz)/len(idx[0]))
        wrms_in_turb_30[i,j] = np.sqrt(np.sum(\
            np.sum(uz[j,:,idx[0],idx[1]]**2,axis=1)/Nz)/len(idx[0]))
        tdisp_in_turb_30[i,j] = (np.sum(\
            np.sum(tdisp[j,:,idx[0],idx[1]],axis=1)/Nz)/len(idx[0]))**2
        tdisp_in_lam_30[i,j] = (np.sum(\
            np.sum(tdisp[j,:,lam[0],lam[1]],axis=1)/Nz)/len(lam[0]))**2
        enstr_in_turb_30[i,j] = np.sum(np.sum(enstr[j,:,idx[0],idx[1]],\
            axis=1)/Nz)/len(idx[0])
        enstr_in_lam_30[i,j] = np.sum(np.sum(enstr[j,:,lam[0],lam[1]],\
            axis=1)/Nz)/len(lam[0])

        print("Finished with file: ", fn, ".")
    cdf_file.close()

for i, fn in enumerate(simdat_files_100):
    cdf_file = MFDataset(fn)
    x = np.array(cdf_file.variables["x"]) 
    y = np.array(cdf_file.variables["y"]) 
    z = np.array(cdf_file.variables["z"])
    t = np.array(cdf_file.variables["t"][:])
    Nx = len(x)
    Ny = len(y)
    Nz = len(z)
    Nt = len(t)
    gx = cdf_file.variables["Gammax"][0]
    gy = cdf_file.variables["Gammay"][0]
    gz = cdf_file.variables["Gammaz"][0]
    dx = gx/Nx
    dy = gy/Ny
    dz = gz/Nz

    #these arrays are indexed by [t,z,y,x]
    ux = np.array(cdf_file.variables["ux"][:])
    uy = np.array(cdf_file.variables["uy"][:])
    uz = np.array(cdf_file.variables["uz"][:])
    temp = np.array(cdf_file.variables["Temp"][:])
    wz =  db.FD6X(uy, Nx, dx) - db.FD6Y(ux, Ny, dy)
    enstr = (db.FD6Y(uz, Ny, dy) - db.FD4Z(uy, Nz, dz))**2 +\
            (db.FD4Z(ux, Nz, dz) - db.FD6X(uz, Nx, dx))**2 + wz**2
    tdisp = np.sqrt(db.FD6X(temp, Nx, dx)**2 + db.FD6Y(temp, Ny, dy)**2 +\
            db.FD4Z(temp, Nz, dz)**2)


    #compute volume fraction where wz >= .1 wzmax
    for j, tval in enumerate(t):
        if len(horiz_urms_100[0,:])-1 < j:
            # ensures arrays have sufficient length
            l = sim_num_files[1]
            vol_frac_100 = add_column(vol_frac_100,l)
            horiz_urms_100 = add_column(horiz_urms_100,l)
            uhrms_in_turb_100 = add_column(uhrms_in_turb_100,l)
            wrms_100 = add_column(wrms_100,l)
            wrms_in_turb_100 = add_column(wrms_in_turb_100,l)
            t_100 = add_column(t_100,l)
            wt_100 = add_column(wt_100,l)
            wt_in_turb_100 = add_column(wt_in_turb_100,l)
            enstr_100 = add_column(enstr_100,l)
            enstr_in_turb_100 = add_column(enstr_in_turb_100,l)
            enstr_in_lam_100 = add_column(enstr_in_lam_100,l)
            tdisp_100 = add_column(tdisp_100,l)
            tdisp_in_turb_100 = add_column(tdisp_in_turb_100,l)
            tdisp_in_lam_100 = add_column(tdisp_in_lam_100,l)

        # finding vertical pencils for turbulent and laminar regions
        
        idx = np.where(Fr[1]*np.abs(np.sum(wz[j,:,:,:],axis=1)/Nz+invRo_100[i]) <= 1)
        lam = np.where(Fr[1]*np.abs(np.sum(wz[j,:,:,:],axis=1)/Nz+invRo_100[i]) > 1)

        # volume average over whole domain
        t_100[i,j] = tval
        horiz_urms_100[i,j] = np.sqrt(np.sum(ux[j,:,:,:]**2 +uy[j,:,:,:]**2)/(Nx*Ny*Nz))
        vol_frac_100[i,j] = len(idx[0])/(Nx*Ny)
        wrms_100[i,j] = np.sqrt(np.sum(uz[j,:,:,:]**2)/(Nz*Ny*Nx))
        wt_100[i,j] = np.sum(uz[j,:,:,:]*temp[j,:,:,:])/(Nz*Ny*Nx)
        enstr_100[i,j] = np.sum(enstr[j,:,:,:])/(Nz*Ny*Nx)
        tdisp_100[i,j] = (np.sum(tdisp[j,:,:,:])/(Nz*Ny*Nx))**2
        
        # compute rms values in turbulent and laminar regions
        wt_in_turb_100[i,j] = np.sum(np.sum(uz[j,:,idx[0],idx[1]]*\
            temp[j,:,idx[0],idx[1]],axis=1)/Nz)/len(idx[0])
        uhrms_in_turb_100[i,j] = np.sqrt(np.sum(\
            np.sum(ux[j,:,idx[0],idx[1]]**2,axis=1)/Nz +
            np.sum(uy[j,:,idx[0],idx[1]]**2,axis=1)/Nz)/len(idx[0]))
        wrms_in_turb_100[i,j] = np.sqrt(np.sum(\
            np.sum(uz[j,:,idx[0],idx[1]]**2,axis=1)/Nz)/len(idx[0]))
        tdisp_in_turb_100[i,j] = (np.sum(\
            np.sum(tdisp[j,:,idx[0],idx[1]],axis=1)/Nz)/len(idx[0]))**2
        tdisp_in_lam_100[i,j] = (np.sum(\
            np.sum(tdisp[j,:,lam[0],lam[1]],axis=1)/Nz)/len(lam[0]))**2
        enstr_in_turb_100[i,j] = np.sum(np.sum(enstr[j,:,idx[0],idx[1]],\
            axis=1)/Nz)/len(idx[0])
        enstr_in_lam_100[i,j] = np.sum(np.sum(enstr[j,:,lam[0],lam[1]],\
            axis=1)/Nz)/len(lam[0])

    print("Finished with file: ", fn, ".")
    cdf_file.close()

# note that some of the elements in these arrays may be zero
# not all files will have the same number of eddy turnover periods in the
# timeseries
for i, fn in enumerate(files_30):
    file = open(fn, 'r')
    urms = np.array([])
    t = np.array([])
    wT = np.array([])
    # load data from text file to np arrays
    for line in file:
        if not line.startswith('#') and line.strip():  # skips headers
            # save data from this line (can choose which data to save)
            line_elements = [float(x) for x in line.split()]
            t = np.append(t,line_elements[1])
            urms = np.append(urms,np.sqrt(line_elements[35]**2 +
            line_elements[36]**2))
            wT = np.append(wT,line_elements[8])

    # perform averaging over eddie turnover times (roughly)
    dt = 2*(np.pi**2)/urms[0]
    idx = np.nonzero(t == t[0])[0]
    # need to change this out for a while loop 
    # needs to iterate over timeseries until it reachs the end
    #for j in range(num_samples):
    j = 0
    while idx[-1] < len(t)-1:
        cent = t[idx[-1]] + dt/2
        idx = np.nonzero(np.abs(t-cent) <= dt/2)[0]
        if len(eInvRo_30[0,:])-1 < j:
            # adds an extra column of zeros to the arrays
            eInvRo_30 = np.hstack((eInvRo_30,np.zeros((num_files[0],1))))
            eInvRo_err_30 = np.hstack((eInvRo_err_30,np.zeros((num_files[0],1))))
            avg_wT_30 = np.hstack((avg_wT_30,np.zeros((num_files[0],1))))
            avg_wT_err_30 = np.hstack((avg_wT_err_30,np.zeros((num_files[0],1))))   
        eInvRo_30[i,j] = np.sum(1/urms[idx])*invRo_30[i]/len(idx)
        eInvRo_err_30[i,j] = np.std(urms[idx])
        avg_wT_30[i,j] = np.sum(wT[idx])/len(idx)
        avg_wT_err_30[i,j] = np.std(wT[idx])
        dt = 2*(np.pi**2)/urms[idx[-1]]
        j += 1
    print("Finished with file: ", fn, ".")
    file.close()

for i, fn in enumerate(files_100):
    file = open(fn, 'r')
    urms = np.array([])
    t = np.array([])
    wT = np.array([])
    # load data from text file to np arrays
    for line in file:
        if not line.startswith('#') and line.strip():  # skips headers
            # save data from this line (can choose which data to save)
            line_elements = [float(x) for x in line.split()]
            t = np.append(t,line_elements[1])
            urms = np.append(urms,np.sqrt(line_elements[35]**2 +
            line_elements[36]**2))

            wT = np.append(wT,line_elements[8])

    # perform averaging over eddie turnover times (roughly)
    dt = 2*(np.pi**2)/urms[0]
    idx = np.nonzero(t == t[0])[0]
    #for j in range(num_samples):
    j = 0
    while idx[-1] < len(t)-1:
        cent = t[idx[-1]] + dt/2
        idx = np.nonzero(np.abs(t-cent) <= dt/2)[0]
        if len(eInvRo_100[0,:])-1 < j:
            # adds an extra column of zeros to the arrays
            eInvRo_100 = np.hstack((eInvRo_100,np.zeros((num_files[1],1))))
            eInvRo_err_100 = np.hstack((eInvRo_err_100,np.zeros((num_files[1],1))))
            avg_wT_100 = np.hstack((avg_wT_100,np.zeros((num_files[1],1))))
            avg_wT_err_100 = np.hstack((avg_wT_err_100,np.zeros((num_files[1],1))))   
        eInvRo_100[i,j] = np.sum(1/urms[idx])*invRo_100[i]/len(idx)
        eInvRo_err_100[i,j] = np.std(urms[idx])
        avg_wT_100[i,j] = np.sum(wT[idx])/len(idx)
        avg_wT_err_100[i,j] = np.std(wT[idx])
        dt = 2*(np.pi**2)/urms[idx[-1]]
        j += 1
    print("Finished with file: ", fn, ".")
    file.close()

# save computed data as txt files
np.savez("vort_frac", t_30=t_30, t_100=t_100,eInvRo_30=eInvRo_30,\
    avg_wT_30=avg_wT_30, eInvRo_err_30=eInvRo_err_30,\
    avg_wT_err_30=avg_wT_err_30,invRo_30=invRo_30, horiz_urms_30=horiz_urms_30,\
    vol_frac_30=vol_frac_30, eInvRo_100=eInvRo_100, avg_wT_100=avg_wT_100,\
    eInvRo_err_100=eInvRo_err_100, avg_wT_err_100=avg_wT_err_100,\
    invRo_100=invRo_100, horiz_urms_100=horiz_urms_100,\
    vol_frac_100=vol_frac_100,\
    wt_in_turb_30=wt_in_turb_30,wt_in_turb_100=wt_in_turb_100,\
    uhrms_in_turb_30=uhrms_in_turb_30, uhrms_in_turb_100=uhrms_in_turb_100,\
    enstr_100=enstr_100,enstr_30=enstr_30,tdisp_30=tdisp_30,tdisp_100=tdisp_100,\
    enstr_in_lam_30=enstr_in_lam_30,enstr_in_lam_100=enstr_in_lam_100,\
    enstr_in_turb_30=enstr_in_turb_30,enstr_in_turb_100=enstr_in_turb_100,\
    tdisp_in_lam_30=tdisp_in_lam_30,tdisp_in_lam_100=tdisp_in_lam_100,\
    tdisp_in_turb_30=tdisp_in_turb_30,\
    tdisp_in_turb_100=tdisp_in_turb_100,wrms_30=wrms_30,wrms_100=wrms_100,\
    wrms_in_turb_30=wrms_in_turb_30,wrms_in_turb_100=wrms_in_turb_100,\
    wt_30=wt_30,wt_100=wt_100)




