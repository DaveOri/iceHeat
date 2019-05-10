import numpy as np
from meteoSI import g, Rair, Rvapor, Cp, Tnull, T_virt_rh, moist_rho_rh
from Takahashi1991 import vgrowth

Cw = 4182. # specific heat of water J/K/kg
Ci = 2093. # specific heat of ice J/K/kg
Lc = 2.501e6 # latent heat of condensation J/kg
Ls = 2.834e6 # latent heat of sublimation J/kg
Lf = Ls - Lc # latent heat of fusion


def integrator_Tconst(p, Tc, RH, Nt, q_hydro_0=0.0, v0=0.0, Tstart=0.0, Tstop=1800, dt=0.001):
    T = -Tnull - Tc # temperature [K]
    rh = RH/100. # 1.1 # relative humidity [Pa/Pa]

    Tv = T_virt_rh(T, rh, p)
    rho_a = moist_rho_rh(p, T, rh, q_hydro_0)

    t = np.arange(Tstart, Tstop, dt)
    m = vgrowth(Tc, t/60.) # single particle mass over time [grams]
    dm = m[1:]-m[:-1] # single particle mass variation over time [grams]
    q_hydro = q_hydro_0 + m*Nt*1e-3 # mass concentration of hydro [kg/m3]
    q_hydro = q_hydro/rho_a # mass mixing ratio [kg/kg]

    dM = dm*Nt/dt # mass of ice condensed from vapour per second in a cubic meter [g m-3 s-1]
    dQv = 1e-3*dM*Ls # Heat released per second in a cubic meter due to ice condensation [J m-3 s-1]
    dQm = dQv/rho_a # Same as above, but for 1 kg of air [J kg-1 s-1]
    #dTemp = dQm/Cp # temperature variation due to air heating [K s-1]
    dTemp = dQm/(Cp + Ci*q_hydro) # temperature variation due to air heating [K s-1] considering hydro absorbing part of heat
    
    Temp = T + dTemp.cumsum()*dt # temperature evolution over time [K]

    Tvp = T_virt_rh(Temp, rh, p) # evolution of virtual temperature [K]
    ab = g*((Tvp-Tv)/Tv - 1e-3*dM.cumsum()*dt) # evolution of buoyant acceleration [m s-2]
    v = v0 + ab.cumsum()*dt # evolution of speed [m s-1]
    s = 0.0 + v.cumsum()*dt + 0.5*ab.cumsum()*dt*dt # evolution of space [m]

    return t, m, Temp, Tvp, ab, v, s