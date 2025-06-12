import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

IC = [1,2,3,4,5,6,7,8,9,10,11]
method = ['fu','pu','pd','pc','pm','po','pv']
xstart = 0
xstop = 1
ylow = -6
yhigh = 6
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
            try:
                t[k] = float(file.readline())
            except:
                t[k] = 0


        idx = np.where(t >= 0.3)
        #print(idx)
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
        plt.ylim(ylow,yhigh)
        p1 = ax1.plot(x,u1[:,pNt],'o-')[0]
        fig.suptitle('IC'+str(IC[i])+' '+figtitle[j])
        fig.tight_layout()
        plt.savefig(fn+'_plot.png',dpi=800)

        def update_frame(frame):
            p1.set_ydata(u1[:,frame])
            #print(p1)
            #plt.savefig(fn+'_frame{:0{}d}.png'.format(frame,5),dpi=800)
            print('Done with frame:', frame)
            return (p1,p1)

        ani = animation.FuncAnimation(fig=fig,
            func=update_frame,frames=pNt,interval=66,blit=True)
        ani.save(fn+'_animation.gif')

        plt.close()
        plt.clf()
        plt.cla()


#plt.savefig(plotfile[i])

#print('\n Finished with ',i+1,'th file\n\n')
    
