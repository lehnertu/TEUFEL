import matplotlib.pyplot as plt
import numpy as np

t,x,y,z,px,py,pz,ax,ay,az=np.loadtxt("StepVay_Trial.txt",delimiter='\t',unpack='True')

plt.figure(0)
plt.subplot(223)
plt.plot(t,z)
plt.xlabel("time(seconds)")
plt.ylabel("Undulator Length(m)")
plt.grid()


plt.subplot(221)
plt.plot(z,x)
plt.xlabel("Undulator Length(m)")
plt.ylabel("x (m)")
plt.grid()

plt.subplot(222)
plt.plot(z,y)
plt.xlabel("Undulator Length(m)")
plt.ylabel("y(m)")
plt.grid()

plt.figure(1)
plt.subplot(221)
plt.plot(z,px)
plt.xlabel("Undulator Length(m)")
plt.ylabel("px (m)")
plt.grid()

plt.subplot(222)
plt.plot(z,py)
plt.xlabel("Undulator Length(m)")
plt.ylabel("py (m)")
plt.grid()

plt.subplot(223)
plt.plot(z,pz)
plt.xlabel("Undulator Length(m)")
plt.ylabel("pz (m)")
plt.grid()

plt.figure(2)
plt.subplot(221)
plt.plot(z,ax)
plt.xlabel("Undulator Length(m)")
plt.ylabel("ax (m)")
plt.grid()

plt.subplot(222)
plt.plot(z,ay)
plt.xlabel("Undulator Length(m)")
plt.ylabel("ay (m)")
plt.grid()

plt.subplot(223)
plt.plot(z,az)
plt.xlabel("Undulator Length(m)")
plt.ylabel("az (m)")
plt.grid()

plt.show()

