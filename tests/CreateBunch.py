import particletoolbox as pb
import numpy as np
import matplotlib.pyplot as plt

"""time is scaled from the longitudinal position as time=z/(betaz*c)"""
def GenerateBunch(nop1,mean_x1,mean_y1,mean_z1,sigma_x1,sigma_y1,sigma_z1,sigma_px1,sigma_py1,mean_E1,sigma_E1):

  x1,y1,z1,px1,py1,pz1,ekin1,Gamma1=pb.GenerateGaussianBunch(nop1,sigma_x1,sigma_y1,sigma_z1,sigma_px1,sigma_py1,mean_E1,sigma_E1)
  gamma=np.mean(Gamma1)
  beta=np.sqrt(1.0-1.0/pow(gamma,2.0))
  c=3.0e8
  if(abs(beta)>0.0):
    time=z1/(beta*c)
  else:
    time = z1*0.0;
  """Shift the distribution that the mean position is (mean_x,mean_y,mean_z)"""
  x1=np.add(x1,mean_x1)
  y1=np.add(y1,mean_y1)
  z1=np.add(z1,mean_z1)
  pz1=np.divide(np.subtract(pz1,0),0.511e6)*beta
  return (time,x1,y1,z1,px1,py1,pz1)

nop=500
mean_x=0.0
mean_y=0.0
mean_z=-0.1
sigma_x=0.001/2.355
sigma_y=0.001/2.355
sigma_z=0.00001
sigma_px=0.001
sigma_py=0.001
mean_E=8e6
sigma_E=10e3

time,x,y,z,px,py,pz=GenerateBunch(nop,mean_x,mean_y,mean_z,sigma_x,sigma_y,sigma_z,sigma_px,sigma_py,mean_E,sigma_E)
f=open("BeamProfile.txt",'w+')
for i in range(0,len(x)):
	f.write(str(time[i])+"\t"+str(x[i])+"\t"+str(y[i])+"\t"+str(z[i])+"\t"+str(px[i])+"\t"+str(py[i])+"\t"+str(pz[i])+"\n")
nop=500
mean_x=0.0
mean_y=0.0
mean_z=-0.1-0.0001
sigma_x=0.001/2.355
sigma_y=0.001/2.355
sigma_z=0.00001
sigma_px=0.001
sigma_py=0.001
mean_E=8e6
sigma_E=10e3
del time, x,y,z,px,py,pz
time,x,y,z,px,py,pz=GenerateBunch(nop,mean_x,mean_y,mean_z,sigma_x,sigma_y,sigma_z,sigma_px,sigma_py,mean_E,sigma_E)
for i in range(0,len(x)):
	f.write(str(time[i])+"\t"+str(x[i])+"\t"+str(y[i])+"\t"+str(z[i])+"\t"+str(px[i])+"\t"+str(py[i])+"\t"+str(pz[i])+"\n")
f.close()
pb.PlotPS(x, y, xlabel='x', ylabel='y')
pb.AnalyzeParticleDistribution(x,y,z,px,py,pz)
plt.show()
