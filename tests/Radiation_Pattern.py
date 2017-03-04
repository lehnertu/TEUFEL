from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
from scipy.integrate import simps
from scipy.interpolate import interp1d


x,y,z,field=np.loadtxt("test1.txt",delimiter='\t',unpack=True)
fig=plt.figure(0)
ax=fig.gca(projection='3d')
surf=ax.scatter(x,y,field,cmap=cm.CMRmap_r,c=y,linewidth=0,antialiased=False)
ax.set_xlabel("x(m)")
ax.set_ylabel("y(m)")
ax.set_zlabel("Power(Watts)")
plt.show()
