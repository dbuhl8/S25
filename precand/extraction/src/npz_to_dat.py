import numpy as np


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
invRo = [[0.5, 0.5], [1,1], [2,2], [3,3], [4,4], [5,8], [8, 0], [10,0]]
invFr = [np.sqrt(30), 10]

labels_100 = [r'$\Omega=0.5$',r'$\Omega=1$',r'$\Omega=2$',r'$\Omega=3$',r'$\Omega=4$',r'$\Omega=8$']
labels_30 = [r'$\Omega=0.5$',r'$\Omega=1$',r'$\Omega=2$',r'$\Omega=3$',r'$\Omega=4$',r'$\Omega=5$',r'$\Omega=8$',r'$\Omega=10$']
cor_label_30 = [1,2,3,4,5,7]

num_files = [len(files_30), len(files_100)]
sim_num_files = [len(simdat_files_30), len(simdat_files_100)]

fn = "vort_frac.npz"
npzfile = np.load(fn)

eInvRo_30=npzfile["eInvRo_30"]
eInvRo_100=npzfile["eInvRo_100"]

eInvRo_err_30=npzfile["eInvRo_err_30"]
eInvRo_err_100=npzfile["eInvRo_err_100"]

avg_wT_30=npzfile["avg_wT_30"]
avg_wT_100=npzfile["avg_wT_100"]

avg_wT_err_30=npzfile["avg_wT_err_30"]
avg_wT_err_100=npzfile["avg_wT_err_100"]

invRo_30=npzfile["invRo_30"]
invRo_100=npzfile["invRo_100"]

horiz_urms_30=npzfile["horiz_urms_30"]
horiz_urms_100=npzfile["horiz_urms_100"]
uhrms_in_turb_30=npzfile["uhrms_in_turb_30"]
uhrms_in_turb_100=npzfile["uhrms_in_turb_100"]

vol_frac_30=npzfile["vol_frac_30"]
vol_frac_100=npzfile["vol_frac_100"]

wt_30=npzfile["wt_30"]
wt_100=npzfile["wt_100"]
wt_in_turb_100=npzfile["wt_in_turb_100"]
wt_in_turb_30=npzfile["wt_in_turb_30"]

wrms_30=npzfile["wrms_30"]
wrms_100=npzfile["wrms_100"]
wrms_in_turb_30=npzfile["wrms_in_turb_30"]
wrms_in_turb_100=npzfile["wrms_in_turb_100"]

enstr_30=npzfile["enstr_30"]
enstr_100=npzfile["enstr_100"]
enstr_in_turb_30=npzfile["enstr_in_turb_30"]
enstr_in_turb_100=npzfile["enstr_in_turb_100"]
enstr_in_lam_30=npzfile["enstr_in_lam_30"]
enstr_in_lam_100=npzfile["enstr_in_lam_100"]

tdisp_30=npzfile["tdisp_30"]
tdisp_100=npzfile["tdisp_100"]
tdisp_in_turb_30=npzfile["tdisp_in_turb_30"]
tdisp_in_turb_100=npzfile["tdisp_in_turb_100"]
tdisp_in_lam_30=npzfile["tdisp_in_lam_30"]
tdisp_in_lam_100=npzfile["tdisp_in_lam_100"]

# write all of the data to a text file for gnuplot

fn = "rotating_sim.dat"
file = open(fn, 'a')

# write file header
file.write("# Data file for Rotating Flows\n\n")
file.write("# Data is arranged in the following indices:\n")
file.write("# Indices 0-"+str(sim_num_files[0]-1)+" correspond to Fr = 0.18\n")
file.write("# Indices "+str(sim_num_files[0]-1)+"-"+\
    str(sim_num_files[0]+sim_num_files[1]-1)+\
    " correspond to Fr = 0.1\n\n")
file.write("# Columns correspond to the following values\n")
file.write("# Column 1: Inverse Rossby Number\n")
file.write("# Column 2: Inverse Froude Number\n")
file.write("# Column 3: Horizontal U rms\n")
file.write("# Column 4: w rms\n")
file.write("# Column 5: wT\n")
file.write("# Column 6: Enstrophy\n")
file.write("# Column 7: TDisp\n")
file.write("# Column 8: Volume Frac (Fr|w_z + invRo| <= 1|)\n")
file.write("# Column 9: Horizontal U rms (in turbulence)\n")
file.write("# Column 10: w rms (in turbulence)\n")
file.write("# Column 11: wT (in turbulence)\n")
file.write("# Column 12: Enstrophy (in turbulence)\n")
file.write("# Column 13: TDisp (in turbulence)\n")
file.write("# Column 14: Enstrophy (in laminar)\n")
file.write("# Column 15: TDisp (in laminar)\n")




file.write("\n\n")

# write F = .18 data
for i, iRo in enumerate(invRo_30):
    # write index header 
    file.write("# Index "+str(i)+"\n") 
    file.write("# File: "+simdat_files_30[i]+"\n") 
    file.write("\n")

    # iRo, iFr, uhrms, wrms, wT, enstr, tdisp, vol_frac, uhrms_in_turb,
    # wrms_in_turb, wT_in_turb, enstr_in_turb, tdisp_in_turb, enstr_in_lam,
    # tdisp_in_lam
    idx_list = np.where(horiz_urms_30[i,:] > 0)
    for j, idx in enumerate(idx_list[0]):
        data_list = [iRo, invFr[0], horiz_urms_30[i,idx], wrms_30[i,idx],\
        enstr_30[i,idx], tdisp_30[i,idx], vol_frac_30[i,idx],\
        uhrms_in_turb_30[i,j],wrms_in_turb_30[i,j], wt_in_turb_30[i,j],\
        enstr_in_turb_30[i,j], tdisp_in_turb_30[i,j], enstr_in_lam_30[i,j],\
        tdisp_in_lam_30[i,j]]
        data_list = [str(x) for x in data_list]
        file.write("        ".join(data_list)+"\n")

    file.write("\n\n") 


for i, iRo in enumerate(invRo_100):
    # write index header 
    file.write("# Index "+str(i+sim_num_files[0])+"\n") 
    file.write("# File: "+simdat_files_100[i]+"\n") 
    file.write("\n")

    # iRo, iFr, uhrms, wrms, wT, enstr, tdisp, vol_frac, uhrms_in_turb,
    # wrms_in_turb, wT_in_turb, enstr_in_turb, tdisp_in_turb, enstr_in_lam,
    # tdisp_in_lam
    idx_list = np.where(horiz_urms_100[i,:] > 0)
    for j, idx in enumerate(idx_list[0]):
        data_list = [iRo, invFr[1], horiz_urms_100[i,idx], wrms_100[i,idx],\
        enstr_100[i,idx], tdisp_100[i,idx], vol_frac_100[i,idx],\
        uhrms_in_turb_100[i,j],wrms_in_turb_100[i,j], wt_in_turb_100[i,j],\
        enstr_in_turb_100[i,j], tdisp_in_turb_100[i,j], enstr_in_lam_100[i,j],\
        tdisp_in_lam_100[i,j]]
        data_list = [str(x) for x in data_list]
        file.write("        ".join(data_list)+"\n")

    file.write("\n\n") 


file.close()

