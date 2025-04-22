import numpy as np
import matplotlib.pyplot as plt

num_points = 40

x = np.linspace(-2,2,num_points)

fig, ax = plt.subplots()

for i, xi in enumerate(x):
    if (np.abs(xi) < .5): 
        # plot in blue
        dx = .5-xi
        px = np.linspace(xi, .5, 100)
        pt = np.linspace(0, 2*dx, 100)
        ax.plot(px, pt, 'red')
    elif (xi > 0):
        dx = xi-.5
        px = np.linspace(xi, .5, 100)
        pt = np.linspace(0, dx, 100)
        ax.plot(px, pt, 'blue')
    else: 
        dx = -2 -xi
        px = np.linspace(xi, -2, 100)
        pt = np.linspace(0, -dx, 100)
        ax.plot(px, pt, 'black')

fig.savefig('prob2_plot.pdf')

