import matplotlib.pyplot as plt
import numpy as np
t,Efield,Bfield=np.loadtxt("test.txt",delimiter='\t',unpack='True')

plt.figure(0)
plt.subplot(221)
plt.plot(t,Efield)
plt.grid()
plt.xlabel("time (seconds)")
plt.ylabel("Electric Field Magnitude (V/m)")
plt.subplot(222)
plt.plot(t,Bfield)
plt.grid()
plt.xlabel("time (seconds)")
plt.ylabel("Poynting Vector Magnitude(dP/dA)")
plt.show()
