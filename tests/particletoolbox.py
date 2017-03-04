from math import *
import numpy as np
from scipy.special import ndtri
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
import time

# rest mass of the electron in eV
mec2 = 510975.0
# speed of light in m/s
SpeedOfLight = 299793000.0

def GenerateRandomCube(nop):
  """
  A generator for coordinate points that are randomly distributed
  within a cube spanning [-1,1] in all axis.
  The first coordinate point always is [0,0,0].
  """
  start=time.time()
  pos = 2.0*np.random.rand(nop,3)-1.0
  pos[0,:]=0.0
  stop=time.time()
  print "used ", str(stop-start)," seconds"
  return pos

def GenerateUniformCube(
  nop,
  mindist = 0.0,
  binnumber = 0):
  """
  A generator for coordinate points that are uniformly distributed
  within a cube spanning [-1,1] in all axis.
  The first coordinate point always is [0,0,0].

  All points are generated randomly
  but are forced to have a minimum distance (mindist) from each other.
  The minimum distance is checked to be smaller than
  the cubic root of 4/nop. In that case 25% of the volume would be filled
  with balls of mindist diameter.

  To speed up the checking, the cube is divided in a number of
  bins (binnumber is the grid size per axis).
  No bin can contain more than one particle and only the neighbouring
  26 bins are checked for too-close particles.

  The binnumber is checked for there being at least three times as many
  bins as points requested. If the binnumber is not given, this is choosen as a default.
  """
  # check mindist
  check = exp(log(4.0/nop)/3.0)
  if mindist == 0.0 :
    mindist = exp(log(1.0/nop)/3.0)
    print "mindist set to ", mindist
  if mindist > check :
    print "mindist limited to ", check
    mindist = check
  # check : nop*2 bins at minimum
  cbroot = int(ceil(exp(log(3.0*nop)/3.0)))
  if binnumber < cbroot :
    binnumber = cbroot
    print "bin number set to ", binnumber
  # all grid cells are empty (-1)
  grid = -np.ones((binnumber, binnumber, binnumber), dtype=int)
  pos = np.zeros([nop, 3], dtype=float)
  # start with the reference particle only
  grid[floor(0.5*binnumber),floor(0.5*binnumber),floor(0.5*binnumber)] = 0
  start=time.time()
  count_rejected_bin = 0
  count_rejected_neighbour = 0
  for i in range(1,int(nop)) :
    accept = False
    while accept==False :
      ptest = 2.0*np.random.rand(3)-1.0	# uniform in [-1,1]
      xbin = int(floor(0.5*(ptest[0]+1.0)*binnumber))
      ybin = int(floor(0.5*(ptest[1]+1.0)*binnumber))
      zbin = int(floor(0.5*(ptest[2]+1.0)*binnumber))
      # print ptest, [xbin, ybin, zbin]
      if grid[xbin, ybin, zbin] == -1 :
	# we have hit an empty bin, so we could accept the point
	accept = True
	# check distance to neighbouring points
	x1=xbin-1
	if x1 < 0 : x1 = 0
	x2=xbin+2
	if x2 > binnumber : x2 = binnumber
	y1=ybin-1
	if y1 < 0 : y1 = 0
	y2=ybin+2
	if y2 > binnumber : y2 = binnumber
	z1=zbin-1
	if z1 < 0 : z1 = 0
	z2=zbin+2
	if z2 > binnumber : z2 = binnumber
	for ix in np.arange(x1,x2) :
	  for iy in np.arange(y1,y2) :
	    for iz in np.arange(z1,z2) :
	      igrid = grid[ix,iy,iz]
	      if  igrid != -1 :
		# we have a particle in a neighbour cell
		# compute distance
		delta = ptest - pos[igrid]
		dist = sqrt(np.dot(delta,delta))
		if dist < mindist :
		  accept = False
		  count_rejected_neighbour += 1
	if accept==True :
	  grid[xbin, ybin, zbin] = i
	  pos[i] = ptest
      else :
	count_rejected_bin += 1
  stop=time.time()
  print "used ", str(stop-start)," seconds"
  print np.count_nonzero(grid+1), " points generated"
  print count_rejected_bin, " points rejected due to hitting an occupied bin"
  print count_rejected_neighbour, " points rejected due to close neighbours"
  return pos, grid

def GenerateDistributionCube6D(nop):
  """
  Generate particles evenly distributed in a 6D cube
  spanning (-1,1) in each direction
  The first 3 dimensions are generated with GenerateUniformCube()
  to obtain a "quiet" distribution
  The last 3 dimenions are simply random numbers.
  The first coordinate point always is [0,0,0, 0,0,0] (reference particle).

  nop ... number of particles requested
  return ... 6 arrays of coordinates
  """
  print "generating spatial distribution ..."
  pos, grid = GenerateUniformCube(nop)
  print "generating momentum distribution ..."
  momentum = GenerateRandomCube(nop)
  x = pos[:,0]
  y = pos[:,1]
  z = pos[:,2]
  px = momentum[:,0]
  py = momentum[:,1]
  pz = momentum[:,2]
  return(x,y,z,px,py,pz)

def MakeGaussianDist(x,sigma):
  """
  redistribute (shift) positions as to make a
  gaussian distribution from a uniform ones
  x ... array of input coordinates range(-1,1)
  sigma ... sigma of output distribution
  return ... array of output coordinates
  """
  return sigma * ndtri(0.5*(x+1.0))

def MakeGaussianDist2D(x,px,sigma_x,sigma_px,dpdx):
  """
  generate a gaussian distribution in a 2D phase-space
  from coordinates distributed evenly on a (-1,1) square
  x, px ... arrays of input coordinates range(-1,1)
  sigma_x = sqrt(beta epsilon) ... sigma [mm] of output distribution
  sigma_px = sqrt(epsilon/beta) ... uncorrelated sigma [mrad] of px
  dpdx = -alpha/beta ... correlation [mrad/mm] to be imposed (positive=divergent)
  return ... array of output coordinates
  """
  rx = sigma_x * ndtri(0.5*(x+1.0))
  rpx = rx*dpdx + sigma_px*ndtri(0.5*(px+1.0))
  # one should limit the output to something like 6 sigma
  return (rx,rpx)

def GenerateGaussianBunch(
  nop,		# number of particles
  sigma_x,	# position spread in mm
  sigma_y,
  sigma_z,
  sigma_px,	# angular spread in mrad
  sigma_py,
  mean_E,	# reference kinetic energy in eV
  sigma_E	# energy spread in eV
):
  """
  Generate a bunch of electrons with a 6D-gaussian distribution
  about a reference particle with a momentum (mean_pz) in z direction.
  First particle is the reference particle which is included in the particle number.

  The particle positions are created with a "quiet start" condition.
  All particles are forced to have a certain minimum distance from each other.

  return values are :
  x,y,z : coordinates in mm
  px, py : angle with respect to the z-axis in mrad
  pz : longitudinal momentum p*c in MeV
  ekin : perticle kinetic energy in MeV
  Gamma : relativistic mass of th reference particle
  """
  (x,y,z,px,py,ekin) = GenerateDistributionCube6D(nop)
  print "assembling phase space ..."
  (x, px) = MakeGaussianDist2D(x,px,sigma_x,sigma_px,0.0)
  (y, py) = MakeGaussianDist2D(y,py,sigma_y,sigma_py,0.0)
  z = MakeGaussianDist(z,sigma_z)
  ekin = mean_E + MakeGaussianDist(ekin,sigma_E)
  # first particle is reference
  g = ekin/mec2 +1.0
  pz = mec2*np.sqrt(g*g-1.0)
  Gamma = mean_E/mec2 + 1.0
  return (x,y,z,px,py,pz,ekin,Gamma)

def WriteParmelaParticles(
  filename,
  nop,
  x, y, z, px, py, ekin,
  freq	# accelerator frequency in Hz
  ) :
  """
  Write nop particle coordinates to a file Particles.txt
  for usage with an INPUT 40 line in PARMELA.

  The accelerator frequency f ist used to translate the z
  coordinate of the particles into a starting phase value
  effectively shifting the particle to z=0.
  Transverse coordinates are trajectory-corrected,
  the particles are translated along x' and y' direction.

  Input coordinates are assumed
  in mm and mrad, energy in eV
  accelerator frequency in Hz

  Output coordinates are :
  x, x', y, y', d(phi), d(W)

  x and y are positions in cm
  x' and y' are px/p in radian
  d(phi) is in degree
  d(W) is in MeV
  """
  # first particle is reference
  eref=ekin[0]
  Gamma = eref/mec2 +1.0
  Beta = sqrt(1.0-1.0/pow(Gamma,2))
  dt = -0.001*z/(Beta*SpeedOfLight) # in seconds
  dphi = dt * freq * 360.0 # in degree phase
  dW = 1e-6*(ekin - eref) # in MeV
  # write to file
  f = open(filename,'w')
  # PARMELA adds its own reference particle, so we drop the first one
  for i in range(1,int(nop)) :
    f.write('{0:10.6f}  {1:10.8f}  {2:10.6f}  {3:10.8f}  {4:10.6f}  {5:10.8f}\n'.format(
      0.1*x[i], 0.001*px[i], 0.1*y[i], 0.001*py[i], dphi[i], dW[i]))
  f.close()

def WriteAstraParticles(
  filename,
  x, y, z, px, py, ekin,
  charge
  ) :
  """
  Write particle coordinates to a file
  for usage with ASTRA.

  Input coordinates are assumed
  in mm and mrad, energy in eV
  charge is the total charge of all perticles in nC

  ASTRA coordinates are :
  -----------------------
  x,y in m
  z in m (relative to first "reference" particle)
  px,py in eV/c (absolute value)
  pz in eV/c (relative to first "reference" particle)
  clock in ns
  charge in nC
  particle type
  status flag
  -----------------------
  """
  nop=len(x)
  # charge per particle
  cpp=charge/float(nop)
  # first particle is reference
  eref=ekin[0]
  Gamma = eref/mec2 +1.0
  Beta = sqrt(1.0-1.0/pow(Gamma,2))
  pzref = Beta*Gamma*mec2
  pz = [sqrt(pow(e+mec2,2)-pow(mec2,2)) for e in ekin]
  pz[0] = pzref
  # write to file
  f = open(filename,'w')
  # write the reference particle
  f.write('{0:14.9f}  {1:14.9f}  {2:14.9f}  {3:14.6f}  {4:14.6f}  {5:14.6f}  {6:6.3f}  {7:14.6g}  {8:2d}  {9:2d}\n'.format(
      0.0, 0.0, 0.0, 0.0, 0.0, pzref, 0.0, cpp, 1, 5))
  # write all other particles
  for i in range(1,int(nop)) :
    f.write('{0:14.9f}  {1:14.9f}  {2:14.9f}  {3:14.6f}  {4:14.6f}  {5:14.6f}  {6:6.3f}  {7:14.6g}  {8:2d}  {9:2d}\n'.format(
	0.001*x[i], 0.001*y[i], 0.001*z[i], 0.001*px[i]*pzref,  0.001*py[i]*pzref,
	pz[i]-pzref, 0.0, cpp, 1, 5))
  f.close()


def ReadAstraParticles(partfile):
  """
  read a particle distribution from an ASTRA run
  return particle coordinates (x,y,z) in mm
  particle momenta (px,py) in mrad pz in MeV
  particle kintic energy ekin in MeV
  relativistic beta and gamma factors of reference particle
  """
  X=[]
  Y=[]
  Z=[]
  PX=[]
  PY=[]
  PZ=[]
  dropped = 0
  f = open(partfile,'r')
  for line in f:
    (x, y, z, px, py, pz, a, b, c, d) = line.split() # to deal with blank
    if int(float(d)) == 5:
      X.append(float(x))
      Y.append(float(y))
      Z.append(float(z))
      PX.append(float(px))
      PY.append(float(py))
      PZ.append(float(pz))
    else:
      dropped += 1
  f.close()
  x=1000.0*np.array(X)
  px=np.array(PX)
  y=1000.0*np.array(Y)
  py=np.array(PY)
  z=1000.0*np.array(Z)
  print len(z)," particles read"
  print dropped," particles dropped"
  zref = z[0]
  print "z(ref.) = ",zref," mm"
  z[0] = 0	# das ist das Referenzteilchen, darauf sind alle anderen bezogen
  pz=np.array(PZ)
  pzref = pz[0]
  print "pz(ref.) = ",1e-6*pzref,"MeV"
  Gamma = sqrt(pow(pzref/mec2,2)+1.0)
  Beta = sqrt(1.0-1.0/pow(Gamma,2))
  pz[0] = 0	# das ist das Referenzteilchen, darauf sind alle anderen bezogen
  pz = pz + pzref
  px = 1e3*px/pz	# change to mrad
  py = 1e3*py/pz	# change to mrad
  ekin = [sqrt(pow(p,2)+pow(mec2,2))-mec2 for p in pz]
  return (x,y,z,px,py,pz,zref,ekin,Beta,Gamma)


def AnalyzeParticleDistribution(x,y,z,px,py,pz):
  """
  analyze the particle distribution
  and print momenta and emittances
  """
  pzref = pz[0]
  bg = pz / mec2
  g = np.sqrt(bg*bg+1.0)
  ekin = (g-1.0)*mec2
  ekin_mean = sum(ekin)/float(len(ekin))
  Gamma = g[0]
  Beta = sqrt(1.0-1.0/pow(Gamma,2))
  print
  print "ref. energy = ", 1e-6*(Gamma-1.0)*mec2, " MeV"
  print "ref. momentum = ", 1e-6*pz[0], " MeV"
  print "mean energy = ", 1e-6*ekin_mean, " MeV"
  print "gamma = ", Gamma
  print

  x_rms = sqrt(np.dot(x,x)/float(len(x)))
  print "x(rms) = ",x_rms," mm"
  px_rms = sqrt(np.dot(px,px)/float(len(px)))
  print "px(rms) = ",px_rms," mrad"
  ex_rms = sqrt((np.dot(x,x)*np.dot(px,px)-pow(np.dot(x,px),2))/pow(float(len(x)),2))
  print "ex(rms) = ",Beta*Gamma*ex_rms," mm mrad"
  print

  y_rms = sqrt(np.dot(y,y)/float(len(y)))
  print "y(rms) = ",y_rms," mm"
  py_rms = sqrt(np.dot(py,py)/float(len(py)))
  print "py(rms) = ",py_rms," mrad"
  ey_rms = sqrt((np.dot(y,y)*np.dot(py,py)-pow(np.dot(y,py),2))/pow(float(len(y)),2))
  print "ey(rms) = ",Beta*Gamma*ey_rms," mm mrad"
  print

  dt = 1.0e9*z/(Beta*SpeedOfLight)
  t_rms = sqrt(np.dot(dt,dt)/float(len(dt)))
  print "t(rms) = ",t_rms," ps"
  dE = ekin-ekin_mean
  e_rms = sqrt(np.dot(dE,dE)/float(len(dt)))
  print "E(rms) = ",1e-3*e_rms," keV"
  ez_rms = sqrt((np.dot(dt,dt)*np.dot(dE,dE)-pow(np.dot(dt,dE),2))/pow(float(len(dt)),2))
  print "ez(rms) = ",1e-3*ez_rms," keV ps"
  print
  return


def AnalyzeParticles(x,y,z,px,py,pz):
  """
  analyze the particle distribution
  return :
  particle number
  mean beam energy in MeV
  RMS beamsize in mm
  bunchlength in ps
  energy spread in keV
  """
  pzref = pz[0]
  bg = pz / mec2
  g = np.sqrt(bg*bg+1.0)
  ekin = (g-1.0)*mec2
  ekin_mean = sum(ekin)/float(len(ekin))
  Gamma = g[0]
  Beta = sqrt(1.0-1.0/pow(Gamma,2))
  x_rms = sqrt(np.dot(x,x)/float(len(x)))
  y_rms = sqrt(np.dot(y,y)/float(len(y)))
  dt = 1.0e9*z/(Beta*SpeedOfLight)
  t_rms = sqrt(np.dot(dt,dt)/float(len(dt)))
  dE = ekin-ekin_mean
  e_rms = sqrt(np.dot(dE,dE)/float(len(dt)))
  return (len(ekin), 1e-6*ekin_mean, x_rms, y_rms, t_rms, 1e-3*e_rms)


def AnalyzeParticlesRMS(x,y,z,px,py,pz):
  """
  analyze the particle distribution
  return RMS beamsize in mm
  normalized transverse emittance in mm.mrad
  bunchlength in ps
  energy spread in keV
  longitudinal emittance in keV.ps
  """
  pzref = pz[0]
  bg = pz / mec2
  g = np.sqrt(bg*bg+1.0)
  ekin = (g-1.0)*mec2
  ekin_mean = sum(ekin)/float(len(ekin))
  Gamma = g[0]
  Beta = sqrt(1.0-1.0/pow(Gamma,2))
  x_rms = sqrt(np.dot(x,x)/float(len(x)))
  ex_rms = sqrt((np.dot(x,x)*np.dot(px,px)-pow(np.dot(x,px),2))/pow(float(len(x)),2))
  y_rms = sqrt(np.dot(y,y)/float(len(y)))
  ey_rms = sqrt((np.dot(y,y)*np.dot(py,py)-pow(np.dot(y,py),2))/pow(float(len(y)),2))
  z_mean = sum(z)/float(len(z))
  dt = 1.0e9*(z-z_mean)/(Beta*SpeedOfLight)
  t_rms = sqrt(np.dot(dt,dt)/float(len(dt)))
  dE = ekin-ekin_mean
  e_rms = sqrt(np.dot(dE,dE)/float(len(dt)))
  ez_rms = sqrt((np.dot(dt,dt)*np.dot(dE,dE)-pow(np.dot(dt,dE),2))/pow(float(len(dt)),2))
  return (x_rms, Beta*Gamma*ex_rms, y_rms, Beta*Gamma*ey_rms, t_rms, 1e-3*e_rms, 1e-3*ez_rms)


def PlotPS(x, y, xlabel='x', ylabel='y'):
  fig = plt.figure(figsize=(6,6),dpi=150)
  # positioning of the plots
  left, width = 0.1, 0.65
  bottom, height = 0.1, 0.65
  rect_dens = [left, bottom, width, height]
  rect_histx = [left, bottom+height+0.02, width, 0.2]
  rect_histy = [left+width+0.02, bottom, 0.2, height]
  # determine the scale range
  xmin=min(x)
  xmax=max(x)
  ymin=min(y)
  ymax=max(y)
  # Histogram in 2D
  H, xticks, yticks = np.histogram2d(x,y,bins=200,
    range=[[xmin,xmax],[ymin,ymax]])
  # wegen Matrickonvention muss H transponiert werden
  # da von links oben beginnend gezeichnet wird, muss man vertikal flippen
  M = np.flipud(H.T)
  # erster Index ist jetzt Energie, wird von unten nach oben dargestellt
  # zweiter Index ist Zeit, wird von links nach rechts dargestellt
  axDens = plt.axes(rect_dens)
  im = axDens.imshow(M, interpolation='nearest',
      aspect=(xticks[0]-xticks[-1])/(yticks[0]-yticks[-1]),
      extent=[xticks[0], xticks[-1], yticks[0], yticks[-1]])
  plt.xlabel(xlabel)
  plt.ylabel(ylabel)
  # plot projections
  axHistx = plt.axes(rect_histx)
  axHistx.xaxis.set_major_formatter(NullFormatter())
  axHistx.yaxis.set_major_formatter(NullFormatter())
  axHistx.set_xlim( axDens.get_xlim() )
  binwidth = (xmax-xmin)/100.0
  bins = np.arange(xmin, xmax, binwidth)
  axHistx.hist(x, bins=bins)
  axHisty = plt.axes(rect_histy)
  axHisty.xaxis.set_major_formatter(NullFormatter())
  axHisty.yaxis.set_major_formatter(NullFormatter())
  axHisty.set_ylim( axDens.get_ylim() )
  binwidth = (ymax-ymin)/100.0
  bins = np.arange(ymin, ymax, binwidth)
  axHisty.hist(y, bins=bins, orientation='horizontal')
  # plt.colorbar()
  return plt

