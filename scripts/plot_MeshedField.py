#!/usr/bin/env python3
# coding=UTF-8

import os.path
import numpy as np
from scipy import constants
from MeshedFields import *
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('file', help='the file name of the watch point HDF5 file')

args = parser.parse_args()

radfile = args.file
radOK = os.path.isfile(radfile)
if not radOK:
  print("file not found")
  sys.exit()

screen = MeshedField.ReadMeshedField(radfile)

print("%d points" % len(screen.points))
print("%d triangles" % len(screen.triangles))
area = screen.MeshArea()
print("total mesh area = %7.3f cm²" % (1.0e4*np.sum(area)))
normals = screen.MeshNormals()
average = np.sum(normals, axis=0)/screen.Np
print("screen normal = %s" % average)

area = screen.MeshArea()
S = [np.linalg.norm(screen.EnergyFlowVector(i)) for i in range(screen.Np)]
Pz = [screen.NormalEnergyFlow(i) for i in range(screen.Np)]

print("peak energy density = %.6g J/m²" % np.max(S))
print("total pulse energy = %.6g µJ" % (1e6*np.dot(area,Pz)))

print()
print("VTK display")
print("-----------")
print("rotate scene with left mouse button + drag")
print("shift view with middle mouse button + drag")
print("zoom with mouse wheel")
print("click on cell to plot waveform")
print("close waveform window before next interaction")
print()

def pick(id):
    if id>0 and id<screen.Np:
        print("cell No. %d pos=%s" % (id,screen.pos[id]))
        print("pointing vector S=%s" % screen.EnergyFlowVector(id))
        screen.ShowFieldTrace(id)

screen.ShowMeshedField(scalars=Pz,scalarTitle="Pz",pickAction=pick,showGrid=False)


