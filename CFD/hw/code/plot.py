import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

fn = ['prob6.dat','prob7.dat','prob8.dat','prob9.dat','prob67check.dat','prob89check.dat']
xstart = [0, -1, 0, -1, 0, -1]
figtitle = ['Problem 6 Solutions','Problem 7 Solutions',\
    'Problem 8 Solutions','Problem 9 Solutions','Upwind Solutions (Sin)',
    'Upwind Solutions (Shock)']
plotfile = ['prob6','prob7','prob8','prob9','upwind_sin','upwind_shock']
ylow = [-1, -0.5, -1, -0.5, -1, -0.5]
ytop = [1, 1.5, 1, 1.5, 1, 1.5]

ptstp = -1

for i in range(6):
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

    minNt = min([len(u1[0,:]),len(u2[0,:]),len(u3[0,:]),len(u4[0,:]),\
        len(u5[0,:]),len(u6[0,:])])
    
    # plot all 6 solutions for this problem 
    fig, ax = plt.subplots(2,3,sharex=True,sharey=True,figsize=(8,4))
    plt.ylim(ylow[i],ytop[i])

    #plot 1
    x32 = np.linspace(xstart[i], 1, 32)
    x128 = np.linspace(xstart[i], 1, 128)

    p1 = ax[0,0].plot(x32, u1[:,ptstp])[0]
    ax[0,0].set_title(r'Nx = 32, CFL = 0.8')
    p2 = ax[0,1].plot(x32, u3[:,ptstp])[0]
    ax[0,1].set_title(r'Nx = 32, CFL = 1.0')
    p3 = ax[0,2].plot(x32, u5[:,ptstp])[0]
    ax[0,2].set_title(r'Nx = 32, CFL = 1.2')
    p4 = ax[1,0].plot(x128, u2[:,ptstp])[0]
    ax[1,0].set_title(r'Nx = 128, CFL = 0.8')
    p5 = ax[1,1].plot(x128, u4[:,ptstp])[0]
    ax[1,1].set_title(r'Nx = 128, CFL = 1.0')
    p6 = ax[1,2].plot(x128, u6[:,ptstp])[0]
    ax[1,2].set_title(r'Nx = 128, CFL = 1.2')

    fig.suptitle(figtitle[i])
    fig.tight_layout()
    plt.savefig(plotfile[i]+'_tstop.png',dpi=800)

    for axisset in ax:
        for axis in axisset:
            axis.set_xlabel(r'x') 
            axis.set_ylabel(r'U') 

    def update_frame(frame):
        fr = frame
        p1.set_ydata(u1[:,fr])
        p2.set_ydata(u3[:,fr])
        p3.set_ydata(u5[:,fr])
        p4.set_ydata(u2[:,fr])
        p5.set_ydata(u4[:,fr])
        p6.set_ydata(u6[:,fr])
        #plt.savefig(plotfile[i])
        plt.savefig(plotfile[i]+'_frame{:0{}d}.png'.format(frame,5),dpi=800)
        print('Done with frame:', frame)
        return (p1, p2, p3, p4, p5, p6)

    #ani = animation.FuncAnimation(fig=fig,
        #func=update_frame,frames=minNt,interval=66,blit=True)
    #ani.save(plotfile[i]+'.gif')
    #plt.savefig(plotfile[i])

    print('\n Finished with ',i+1,'th file\n\n')
    
