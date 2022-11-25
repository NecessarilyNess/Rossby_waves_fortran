import matplotlib.pyplot as plt
#import csv
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter
time = 21

with open("streamfunctionn.csv") as file_name:
    array = np.loadtxt(file_name, delimiter=",")
    #print(array)
array1 = np.zeros((64,64,time))
for ii in range(time):
    array1[:,:,ii] = array[4096*ii:4096*(ii+1)].reshape(64,64)


x = np.linspace(-5,5,64)
y = np.linspace(-5,5,64)
xx, yy = np.meshgrid(x, y)
Y, X = np.meshgrid(y, x)
fig = plt.figure(figsize = (8,8))
ax = fig.subplots()

def init_func():
    plt.cla()
    plot = plt.contourf(X, Y, array1[:,:,ii], 50, cmap="coolwarm")
    cbar = plt.colorbar(plot, pad=0.05)
    cbar.ax.set_ylabel('streamfunction')

def update_plot(ii):
    plt.cla()
    ax.clear()
    plt.xlabel('X')
    plt.ylabel('Y')
    #plt.title(f'Streamfunction of a single plane wave, k=-1, l=-0.7, phase=0, t={ii}')
    plt.title(f'Streamfunction of 10 plane waves chosen randomly from the unit circle, t={ii}')
    plt.contourf(X, Y, array1[:,:,ii], 50, cmap="coolwarm")

anim = FuncAnimation(fig,
                        update_plot,
                        frames=np.arange(0, time),
                        init_func=init_func)

writergif = PillowWriter(fps=3)
anim.save('streamfunction10R.gif', writer=writergif)
