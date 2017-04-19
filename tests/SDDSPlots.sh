#usr/bin/bash

#sddsplot -columnNames=z,'(x,y,px,py,pz,Ax,Ay,Az)' -separate -graphic=symbol,vary -legend trajectory.sdds
sddsenvelope trajectory.sdds traj.env -copy=z -mean=x,y,z,px,py,pz,Ax,Ay,Az,gamma -standarddeviation=x,y,z,px,py,pz,gamma
#sddsplot -columnNames=z,'(x,y)' -split=pages -separate -graphic=dot,vary -legend time-trajectory.sdds
sddsplot -columnNames=z,'(xMean,yMean,zMean,pxMean,pyMean,pzMean,gammaMean,AxMean,AyMean,AzMean,xStDev,yStDev,zStDev,pxStDev,pyStDev,pzStDev,gammaStDev)' -separate -graphic=line,vary -legend traj.env
sddsplot -columnNames=t,'(Ex,Ey,Ez,Bx,By,Bz,PoyntingVector)' -separate -graphic=line,vary -legend  radiation.sdds
sddsfft -column=t radiation.sdds -psdOutput  fft-radiation.sdds
sddsplot -columnNames=f,'(FFTEx,FFTEy,FFTEz)' -separate -scales=0.5e12,12e12,1,1 fft-radiation.sdds
sddsinteg -integrate=PSDEx -versus=f -printFinal fft-radiation.sdds Power.sdds
sddscontour radiation@grid.sdds -xyz=z,x,Ex -topLine=@time -shade -mapShade=0,1


