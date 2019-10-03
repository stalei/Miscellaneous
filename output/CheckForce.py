
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import *
np.set_printoptions(threshold=np.inf)
mydt = np.dtype([
                   ('col1', np.double),
                   ('col2', np.double),
                   ('col3', np.double),
				   ('col4', np.double)
                   ])

points = np.fromfile('PF.bin', dtype=mydt)
#####################################
R0=1.9

######################################
G=6.67e-8
N=len(points)
print N
R=range(0,N)
Fi=range(0,N)
Z=range(0,N)
FZ=range(0,N)
Ffit=range(0,N)

for i in range(0, N):
	R[i]=points[i][0]
	Fi[i]=points[i][1]
	Z[i]=points[i][2]
	FZ[i]=points[i][3]

r=np.array(R)
fi=np.array(Fi)
z=np.array(Z)
Fz=np.array(FZ)
X=np.cos(fi)*R
Y=np.sin(fi)*R
print Fz

fig = plt.figure(1)
ax = fig.add_subplot(111, projection='3d')

ax.quiver(X, Y, Z, 0, 0, -Fz, length=1.0e46)

fig.set_size_inches(12,9)
ax.set_xlabel('X (kpc)')
ax.set_ylabel('Y (kpc)')
ax.set_zlabel('Z (kpc)')

print min(Fz),max(Fz)
plt.show()