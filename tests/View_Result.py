import matplotlib.pyplot as plt
import numpy as np

x,y,z=np.loadtxt("trajectory.txt",delimiter='\t',unpack='True')
plt.figure(0)
plt.scatter(z,x)
plt.xlabel("Displacement in z(m)")
plt.ylabel("Displacement in y(m)")
plt.grid()
plt.figure(1)
plt.scatter(z,y)
plt.xlabel("Displacement in z(m)")
plt.ylabel("Displacement in x(m)")
plt.grid()

plt.show()
