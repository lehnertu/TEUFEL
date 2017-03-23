#usr/bin/bash

sddsplot -columnNames=z,'(x,y,px,py,pz,Ax,Ay,Az)' -separate -graphic=symbol,vary -legend trajectory.sdds
sddsenvelope trajectory.sdds traj.env -copy=z -mean=x,y,z,px,py,pz,Ax,Ay,Az -standarddeviation=x,y,z,px,py,pz
sddsplot -columnNames=z,'(x,y)' -split=page -separate -graphic=dot,vary -legend time-trajectory.sdds
sddsplot -columnNames=z,'(xMean,yMean,zMean,pxMean,pyMean,pzMean,xStDev,yStDev,zStDev,pxStDev,pyStDev,pzStDev)' -separate -graphic=line,vary -legend traj.env
