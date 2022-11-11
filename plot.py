import matplotlib.pyplot as plt
#import csv
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter
time = 81

with open("data.csv") as file_name:
    array = np.loadtxt(file_name, delimiter=",")
    #print(array)
array1 = np.zeros((64,64,2*time))
for ii in range(time):
    array1[:,:,2*ii] = array[4096*ii:4096*(ii+1),0].reshape(64,64)
    array1[:,:,2*ii+1] = array[4096*ii:4096*(ii+1),1].reshape(64,64)


x = np.linspace(-5,5,64)
y = np.linspace(-5,5,64)
xx, yy = np.meshgrid(x, y)
Y, X = np.meshgrid(y, x)
fig = plt.figure(figsize = (12,12))

def init_func():
    plt.cla()

def update_plot(ii):
    plt.cla()
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title(f'Velocity of an isotropic spectrum of 25 random Rossby waves, t={ii}')
    plt.quiver(xx, yy, array1[:,:,2*ii], array1[:,:,2*ii+1], np.arctan2(array1[:,:,2*ii],array1[:,:,2*ii+1]))

anim = FuncAnimation(fig,
                        update_plot,
                        frames=np.arange(0, time),
                        init_func=init_func)

writergif = PillowWriter(fps=3)
anim.save('velocityanim25R.gif', writer=writergif)

#x = np.linspace(-5,5,64)
#y = np.linspace(-5,5,64)
#X, Y = np.meshgrid(x, y)
#plt.figure(figsize=(12, 12))
#plt.quiver(X, Y, array1[:,:,0], array1[:,:,1], np.arctan2(array1[:,:,0],array1[:,:,1]))
#plt.xlabel('X')
#plt.ylabel('Y')
#plt.savefig('Velocity0.png')
