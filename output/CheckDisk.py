from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
mydt = np.dtype([
                   ('col1', np.double),
                   ('col2', np.double),
                   ('col3', np.double)
                   ])


coord = np.fromfile('Disk.out', dtype=mydt)

#N=10000
N=len(coord)
print N
R=range(0,N)
Fi=range(0,N)
Z0=range(0,N)

for i in range(0, N):
	R[i]=coord[i][0]
	Fi[i]=coord[i][1]
	Z0[i]=coord[i][2]

Ra=np.array(R)
Fia=np.array(Fi)


X=Ra*np.cos(Fia)
Y=Ra*np.sin(Fia)
Z=np.array(Z0)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.scatter(X,Y,Z,c='b', marker='.')

fig.set_size_inches(14,8)
ax.set_xlabel('X (kpc)')
ax.set_ylabel('Y (kpc)')
ax.set_zlabel('Z (kpc)')

plt.show()