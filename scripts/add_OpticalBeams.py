#!/usr/bin/env python3
# coding=UTF-8

import copy
import os.path
import argparse
import numpy as np
import matplotlib.pyplot as plt
from screen import *

parser = argparse.ArgumentParser()
parser.add_argument('files', nargs='+', default=[], help='a list of files with recorded screen fields')
args = parser.parse_args()

is_first = True

for filename in args.files:
    if os.path.isfile(filename):
    
        print("reading ",filename)
        screen = TeufelScreen.read(filename)
        print()
        
        if is_first:
            # create a screen exactly like the one from the first file
            # but with zero fields in all arrays
            (Nx,Ny) = screen.shape()
            (EVec, BVec) = screen.getFieldTrace(0,0)
            output = copy.deepcopy(screen)
            NullVec = np.zeros_like(EVec)
            for ix in range(Nx):
                for iy in range(Ny):
                    output.setFieldTrace(ix,iy,NullVec,NullVec)
            # do the zeroing only once
            is_first = False
            
        # sum up the fields
        for ix in range(Nx):
            for iy in range(Ny):
                (ESum, BSum) = output.getFieldTrace(ix,iy)
                (EVec, BVec) = screen.getFieldTrace(ix,iy)
                output.setFieldTrace(ix,iy,EVec+ESum,BVec+BSum)
                

