#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  9 09:11:54 2019

@author: dori
"""

import numpy as np
import matplotlib.pyplot as plt

p = {} # dictionary of parameters taken from Takahashi et al. 1991

"""
I had to relax the condition of maximum overall time for all of the crystals to
be the same to allow smooth interpolation
"""

p[ -3.7] = [1.3e-8,  1.4, 30]
p[ -5.3] = [8.4e-9,  2.0, 30]
p[ -8.6] = [2.5e-8, 1.52, 30]
p[-10.6] = [3.3e-8, 1.51, 30]
p[-12.2] = [6.3e-8, 1.28,  7, 1.1e-8, 2.16, 30]
p[-14.4] = [6.4e-8, 1.63,  4, 2.3e-8, 2.37, 30]
p[-16.5] = [8.1e-8, 1.36,  5, 3.5e-8, 1.91, 30]
p[-18.2] = [4.2e-8, 1.49, 30]
p[-20.1] = [3.8e-8, 1.44, 30]
p[-22.0] = [3.4e-8, 1.42, 30]

k = np.array(list(p.keys()))

def growT(t, a0, b0, t0, a1=None, b1=None, t1=None):
  if t <= t0:
    return a0*np.power(t,b0)
  elif t1 is not None:
    if t <= t1:
      return a1*np.power(t,b1)
    else:
      return np.nan
  else:
    return np.nan

def growth(T, t):
  """
  T is in degrees Celsius
  t is in minutes
  returns mass in grams
  """
  try:
    b = k[k<=T].max()
    u = k[k>=T].min()
    if b == u:
      return growT(t, *p[b])
  except:
    return np.nan
  gb, gu = growT(t, *p[b]), growT(t, *p[u])
#  if abs(b-T)<abs(u-T):
#    return gb
#  else:
#    return gu
  return np.interp(T, [b, u], [gb, gu])

vgrowth = np.vectorize(growth)

if __name__ is '__main__':
  temps = np.linspace(-2.0,-24, 1000)
  mins = np.array([3, 5, 7, 10, 12, 15, 20, 25, 30])
  xt, xm = np.meshgrid(temps, mins)
  masses = vgrowth(xt, xm)
  plt.semilogy(temps, masses.T)
  plt.legend([str(i) for i in mins], title="time [min]")
  plt.xlabel('Temperature [deg C]')
  plt.ylabel('Ice crystal mass   [g]')
  plt.title('Reproduction of plots from Takahashi 1991')
  plt.grid()
  plt.xlim([-2, -24])