import matplotlib.pyplot as plt
import numpy as np

t,Efield,Bfield=np.loadtxt("test.txt",delimiter='\t',unpack='True')


plt.figure(0)
plt.plot(t,Efield)
plt.grid()
plt.xlabel("time (seconds)")
plt.ylabel("Electric Field Magnitude (V/m)")

plt.figure(1)
plt.plot(t,Bfield)
plt.grid()
plt.xlabel("time (seconds)")
plt.ylabel("Magnetic Field Magnitude (T)")

plt.show()
