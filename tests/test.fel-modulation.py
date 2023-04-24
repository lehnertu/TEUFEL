#!/usr/bin/env python3

import os, sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate
import scipy.constants as sc

TEUFEL_ROOT = os.getcwd()
sys.path.append(TEUFEL_ROOT+'/scripts')

from teufel import *

mec2 = sc.m_e * sc.c**2 / sc.e
λ = 1e-4
ω = 2*sc.pi*sc.c/λ
k = ω/sc.c
λU = 0.1
kU = 2*sc.pi/λU
E_kin = 30.0e6
γ = E_kin / mec2 + 1
print(f'γ = {γ:.3f}')
K = np.sqrt(2.0*(λ/λU * 2*γ**2 -1.0))
print(f'K = {K:.3f}')

# particle velocity components inside the undulator
βx = lambda z: K/γ * np.sin(kU*z)
βz = lambda z: np.sqrt(1.0 - 1.0/γ**2 - βx(z)**2)
βz_star = 1.0 - λ/λU
# time inside an undulator period
f = lambda z: 1.0/βz(z)
ct = lambda z: scipy.integrate.quad(f, 0, z)[0]

# time dependency of the wave electric field
Ex = lambda z: np.sin(k*ct(z) - k*z)
# overlap between particle and wave
f2 = lambda z: np.sin(kU*z)*Ex(z)
FF = scipy.integrate.quad(f2, 0, λU)[0] / λU

E0 = 1e6
Δγ = sc.e/(sc.m_e*sc.c**2)*E0*K/γ*λU*FF
print(f'Δγ = {Δγ:.6f}')

initial = TeufelWatch.read('fel-modulation_start.h5')
print(f'start at z = {np.mean(initial.z):.6f} m')
final = TeufelWatch.read('fel-modulation_stop.h5')
print(f'stop at z = {np.mean(final.z):.6f} m')
dg = 0.1*(final.gamma-initial.gamma)
mod = 0.5*(np.max(dg) - np.min(dg))
print(f'Δγ = {mod:.6f} deviating by {100.0*np.abs(mod-Δγ)/Δγ:.2f}%')

errors = 0

if np.abs(mod-Δγ)/Δγ > 0.01: errors += 1
    
if errors==0:
    print("\033[1;32mOK.\033[0m")
else:
    print("\033[1;31mERROR.\033[0m")

sys.exit(errors)
