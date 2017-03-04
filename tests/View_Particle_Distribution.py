import matplotlib.pyplot as plt
import numpy as np

time,x,y,z,px,py,pz=np.loadtxt("BeamProfile.txt",delimiter='\t',unpack='True')

fig=plt.figure()

ax1=fig.add_subplot(221)
ax1.scatter(z,x)

ax2=fig.add_subplot(222)
ax2.scatter(z,y)

ax3=fig.add_subplot(223)
ax3.scatter(x,y)

plt.show()
