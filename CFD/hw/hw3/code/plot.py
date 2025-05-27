import numpy as np
import matplotlib.pyplot as plt

fn = ['prob6.dat','prob7.dat','prob8.dat','prob9.dat']
xstart = [0, -1, 0, -1]
figtitle = ['Problem 6 Solutions','Problem 7 Solutions',\
    'Problem 8 Solutions','Problem 9 Solutions']
plotfile = ['prob6.pdf','prob7.pdf','prob8.pdf','prob9.pdf']

ptstp = 0

for i in range(4):
    file = open(fn[i],'r')
    ln = 1

    Nx, Nt = (int(x) for x in file.readline().split())
    u1 = np.zeros((Nx,Nt))
    # every time it catches a blank line it skips to next array to take input
    u = np.array([float(x) for x in file.readline().split()]).T
    fr = 0
    while len(u) > 0:
        u1[:,fr] = u
        fr += 1
        u = [float(x) for x in file.readline().split()]
    print('Done with first block')

    Nx, Nt = (int(x) for x in file.readline().split())
    u2 = np.zeros((Nx,Nt))
    u = np.array([float(x) for x in file.readline().split()]).T
    fr = 0
    while len(u) > 0:
        u2[:,fr] = u
        fr += 1
        u = [float(x) for x in file.readline().split()]
    print('Done with second block')
            
    Nx, Nt = (int(x) for x in file.readline().split())
    u3 = np.zeros((Nx,Nt))
    u = np.array([float(x) for x in file.readline().split()]).T
    fr = 0
    while len(u) > 0:
        u3[:,fr] = u
        fr += 1
        u = [float(x) for x in file.readline().split()]

    print('Done with third block')
    Nx, Nt = (int(x) for x in file.readline().split())
    u4 = np.zeros((Nx,Nt))
    u = np.array([float(x) for x in file.readline().split()]).T
    fr = 0
    while len(u) > 0:
        u4[:,fr] = u
        fr += 1
        u = [float(x) for x in file.readline().split()]

    print('Done with fourth block')
    Nx, Nt = (int(x) for x in file.readline().split())
    u5 = np.zeros((Nx,Nt))
    u = np.array([float(x) for x in file.readline().split()]).T
    fr = 0
    while len(u) > 0:
        u5[:,fr] = u
        fr += 1
        u = [float(x) for x in file.readline().split()]

    print('Done with fifth block')
    Nx, Nt = (int(x) for x in file.readline().split())
    u6 = np.zeros((Nx,Nt))
    u = np.array([float(x) for x in file.readline().split()]).T
    fr = 0
    while len(u) > 0:
        u6[:,fr] = u
        fr += 1
        try: 
            u = [float(x) for x in file.readline().split()]
        except: 
            break
    print('Done with 6th block')

    # plot all 6 solutions for this problem 
    fig, ax = plt.subplots(2,3,sharex=True,sharey=True,figsize=(16,9))
    plt.ylim(-1,1)

    #plot 1
    x32 = np.linspace(xstart[i], 1, 32)
    x128 = np.linspace(xstart[i], 1, 128)

    ax[0,0].plot(x32, u1[:,ptstp])
    ax[0,0].set_title(r'Nx = 32, CFL = 0.8')
    ax[0,1].plot(x32, u3[:,ptstp])
    ax[0,1].set_title(r'Nx = 32, CFL = 1.0')
    ax[0,2].plot(x32, u5[:,ptstp])
    ax[0,2].set_title(r'Nx = 32, CFL = 1.2')
    ax[1,0].plot(x128, u2[:,ptstp])
    ax[1,0].set_title(r'Nx = 128, CFL = 0.8')
    ax[1,1].plot(x128, u4[:,ptstp])
    ax[1,1].set_title(r'Nx = 128, CFL = 1.0')
    ax[1,2].plot(x128, u6[:,ptstp])
    ax[1,2].set_title(r'Nx = 128, CFL = 1.2')

    fig.suptitle(figtitle[i])

    for axisset in ax:
        for axis in axisset:
            axis.set_xlabel(r'x') 
            axis.set_ylabel(r'U') 

    plt.savefig(plotfile[i])

    print('\n Finished with ',i+1,'th file\n\n')
    
