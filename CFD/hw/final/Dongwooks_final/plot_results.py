import matplotlib.pyplot as plt
import numpy as np

fn_grid = ['grid_sod_fog.dat', 'grid_sod_plm_minmod.dat', 'grid_sod_ppm_minmod.dat']
fn_slug = ['slug_sod_fog_tot.dat', 'slug_sod_plm_minmod_tot.dat',
    'slug_sod_ppm_minmod_tot.dat']
sbplt_title = ['FOG + HLL', 'PLM + Minmod + Hll', 'PPM + Minmod + Hll']

fig, ax = plt.subplots(3,1)

for i in range(3):
    file = open(fn_grid[i],'r')
    nt, nx, xstart, xstop = [float(x) for x in file.readline().split()]
    nt = int(nt)
    nx = int(nx)
    tstop = 0.2

    t = np.zeros(nt)
    x = np.linspace(xstart, xstop, nx+1, True) + (xstop-xstart)/(2.*nx)
    x = x[:nx]

    dens = np.zeros((nt,nx))
    vel = np.zeros_like(dens)
    pres = np.zeros_like(dens)
    gamc = np.zeros_like(dens)
    game = np.zeros_like(dens)
    eint = np.zeros_like(dens)

    file = open(fn_slug[i],'r')

    tidx = 0
    xidx = 0

    for line in file:
        tval, xval, d, v, p, c, e, ei = [float(x) for x in line.split()]
        if not tval in t: 
            if tstop in t:
                break
            tidx += 1
            t[tidx] = tval
        xidx = np.argmax(x>=xval)
        dens[tidx,xidx] = d
        vel[tidx,xidx] = v
        pres[tidx,xidx] = p
        gamc[tidx,xidx] = c
        game[tidx,xidx] = e
        eint[tidx,xidx] = ei

    print('Finished reading file: ',fn_slug[i])

    tidx = np.argmax(t == tstop)
    ax[i].plot(x, dens[tidx,:], 'ro-', x, vel[tidx,:], 'go-', x, pres[tidx,:], 'bo-')
    ax[i].set_title(sbplt_title[i])

plt.show()

    
