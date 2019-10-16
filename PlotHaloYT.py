import h5py as h5
import numpy as np
from matplotlib.legend_handler import HandlerLine2D
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm



f=h5.File("catalog.0.h5","r")

datasetNames = [n for n in f.keys()]
for n in datasetNames:
	print(n)

particle_identifier=f['particle_identifier']
particle_position_x=f['particle_position_x']
particle_position_y=f['particle_position_y']
particle_position_z=f['particle_position_z']
particle_mass=f['particle_mass']


ID0=np.array(particle_identifier)
x0=np.array(particle_position_x)
y0=np.array(particle_position_y)
z0=np.array(particle_position_z)
mass0=np.array(particle_mass)
mass0/=1.989e33# 5.02785e34

print(mass0)

x=x0[mass0>1.32e12]/(3.08568e24)
y=y0[mass0>1.32e12]/(3.08568e24)
z=z0[mass0>1.32e12]/(3.08568e24)
ID=ID0[mass0>1.32e12]

print(len(ID))

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
label =str(particle_identifier)

for i in range(len(ID)):
	#ax.scatter(particle_position_x[i],particle_position_y[i],particle_position_z[i],c='lime', alpha=0.6, marker='+',s=1)
	ax.scatter(x[i],y[i],z[i],c='black', alpha=0.6, marker='+',s=7)
	ax.text(x[i],y[i],z[i],'%s' %(str(ID[i])), size=7,c='grey')

ax.set_xlabel('X (Mpc)')
ax.set_ylabel('Y (Mpc)')
ax.set_zlabel('Z (Mpc)')

plt.show()