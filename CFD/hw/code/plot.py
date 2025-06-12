import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

IC = [1,2,3,4,5,6,7,8,9,10,11]
method = ['fu','pu','pd','pc','pm','po','pv']
xstart = 0
xstop = 1
figtitle = ['FOG','PLM + Upwind SL',\
    'PLM + Downwind SL','PLM + Centered SL','PLM + Minmod SL',\
    'PLM + MC SL','PLM + VanLeers SL']


ptstp = -1

for i in range(11):
    for j in range(7):
        fn = "IC"+str(IC[i])+"Method"+method[j]
        file = open("U_"+fn+".dat",'r')
        print('Reading file: '+fn)
        ln = 1

        #Nx, Nt = (int(x) for x in file.readline().split())
        Nx, Nt = (36, 1000)
        u1 = np.zeros((Nx,Nt))
        t = np.zeros(Nt)
        # every time it catches a blank line it skips to next array to take input
        u = np.array([float(x) for x in file.readline().split()]).T
        fr = 0
        while len(u) > 0:
            u1[:,fr] = u
            fr += 1
            try:
                u = [float(x) for x in file.readline().split()]
            except: 
                u = []

        file.close()

        file = open("t_"+fn+".dat",'r')
        for k in range(Nt):
            t[k] = float(file.readline())

        idx = np.where(t >= 0.3)
        print(idx)
        if (len(idx[0]) == 0):
            nidx = np.where(u1 > 0)
            pNt = max(nidx[1])
        else:
            pNt = idx[0][-1]-1
        dx = 1./32.
        x = np.linspace(0-2*dx, 1+2*dx, Nx)
        # plot all 6 solutions for this problem 
        fig = plt.figure()
        ax1 = fig.add_subplot()
        #plt.ylim(ylow[i],ytop[i])

        ax1.plot(x,u1[:,pNt])

        fig.suptitle('IC'+str(IC[i])+' '+figtitle[j])
        fig.tight_layout()
        plt.savefig(fn+'_plot.png',dpi=800)
        plt.close()
        plt.clf()
        plt.cla()


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

#print('\n Finished with ',i+1,'th file\n\n')
    
