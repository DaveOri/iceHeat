import numpy as np
from meteoSI import g, Rair, Rvapor, Cp, Tnull, T_virt_rh, moist_rho_rh

Cw = 4182. # specific heat of water J/K/kg
Ci = 2093. # specific heat of ice J/K/kg
Lc = 2.501e6 # latent heat of condensation J/kg
Ls = 2.834e6 # latent heat of sublimation J/kg
Lf = Ls - Lc # latent heat of fusion


def integrator_Tconst(p, Tc, RH, Nt, q_hydro=0.0, v0=0.0, Tstop=1800, dt=0.001):
    T = -Tnull - Tc # temperature [K]
    rh = RH/100. # 1.1 # relative humidity [Pa/Pa]

    Tv = T_virt_rh(T, rh, p)
    rho_a = moist_rho_rh(p, T, rh, q_hydro)

    time = np.array([0., 3, 5, 10, 12, 15, 25, 30])*60. # seconds
    mass = np.array([0.,.5, .8,  5, 9, 15, 50, 70])*1e-6 # grams

    t = np.arange(0.,Tstop, dt)
    #m1 = np.power(10.0, np.interp(np.log10(t), np.log10(time), np.log10(mass)))
    m = np.interp(t, time, mass)

    dm = m[1:]-m[:-1]
    #dm1 = m1[1:]-m1[:-1]

    dM = dm*Nt/dt # mass of ice condensed from vapour per second in a cubic meter [g m-3 s-1]
    dQv = 1e-3*dM*Ls # Heat released per second in a cubic meter due to ice condensation [J m-3 s-1]
    dQm = dQv/rho_a # Same as above, but for 1 kg of air [J kg-1 s-1]
    dTemp = dQm/Cp # temperature variation due to air heating [K s-1]
    Temp = T + dTemp.cumsum()*dt # temperature evolution over time [K]

    Tvp = T_virt_rh(Temp, rh, p) # evolution of virtual temperature [K]
    ab = -g*(Tv-Tvp)/Tvp # evolution of buoyant acceleration [m s-2]
    v = v0 + ab.cumsum()*dt # evolution of speed [m s-1]
    s = 0.0 + v.cumsum()*dt + 0.5*ab.cumsum()*dt*dt # evolution of space [m]
    #s1 = v.cumsum()*dt # evolution of space [m]

    return t, m, Temp, Tvp, ab, v, s