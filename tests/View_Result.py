import matplotlib.pyplot as plt
import numpy as np

x,y,z=np.loadtxt("test.txt",delimiter='\t',unpack='True')
plt.figure(0)
plt.plot(z,x)
plt.grid()
plt.figure(1)
plt.plot(z,y)
plt.grid()

plt.show()
