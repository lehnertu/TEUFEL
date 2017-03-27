import matplotlib.pyplot as plt
import numpy as np

t,x,y,z,px,py,pz=np.loadtxt("InitialDistribution.txt",delimiter='\t', unpack='True')
plt.figure(0)
plt.plot(x,y,'.')
t1,x1,y1,z1,px1,py1,pz1=np.loadtxt("trajectory.txt",delimiter='\t', unpack='True')

plt.plot(x1,y1,'*')

plt.show()
