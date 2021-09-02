#--------------------------------------------------------------
#                        exoasteroids
#               Model exoplanet and exo-asteroid transits!
#                     Aug 2021, Natalia Guerrero
#---------------------------------------------------------------

import batman
import numpy as np
import matplotlib.pyplot as plt
import sys
import astropy.units as u

#---------------------------------------------------------------
#             Set up planet parameters
#---------------------------------------------------------------

# TODO: be fancy and pull from tango input.py
rstar = 1.

# Keep P, a/Rstar, e, omega constant, vary t0, Rp/Rstar, inc.; Pjup = 100d

params = batman.TransitParams()
params.t0 = 0.                      #time of inferior conjunction
params.per = 100                     #orbital period (Kepler-10 c: 45.2942679 d)
params.rp = 0.1 * rstar             #planet radius (in units of stellar radii)
params.a = 3                         #semi-major axis (in units of stellar radii)
params.inc = 90                      #orbital inclination (in degrees)
params.ecc = 0.                      #eccentricity
params.w = 90.                       #longitude of periastron (in degrees)
params.u = [0.1, 0.3]            #limb darkening coefficients [u1, u2]
params.limb_dark = "quadratic"       #limb darkening model

t = np.arange(0, 1000, 0.0001)
# make a data point every 2.4 hours

m = batman.TransitModel(params, t)    #initializes model
plan_flux = m.light_curve(params) - 1.          #calculates light curve


#---------------------------------------------------------------
#             Set up asteroid ensemble parameters
#               Greeks: leading
#               Trojans: trailing
#---------------------------------------------------------------

n_trojans = 3  #100     # this is a starting point
n_greeks = 3  #100

dpd = float(params.per) / 360.  # days per degree

rarstar = 1e-4   # r_asteroid / r_star 

trojans = []
greeks = []

for trojan in range(n_trojans):
  t0_scatter = np.random.uniform(-1,1,1)[0]*20*dpd
  #random scatter in inclination
  inc_scatter = np.random.uniform(-1,1,1)[0]*13.7 # 13.7 deg is mean for Jupiter trojans 
  
  t_params = batman.TransitParams()
  # t0 = t0_planet - trojan phase +/- random scatter
  t_params.t0 = params.t0 - 60*dpd + t0_scatter
  t_params.per = 100.
  t_params.rp = rarstar * rstar # * 1e-4
  t_params.a = 3.   # 3 rstar
  t_params.inc = 90. + inc_scatter
  t_params.ecc = 0.
  t_params.w = 90.
  t_params.u = [0.1,0.3]            #limb darkening coefficients [u1, u2]
  t_params.limb_dark = "quadratic"       #limb darkening model
  t_m = batman.TransitModel(t_params, t)
  t_flux = t_m.light_curve(t_params) - 1.  # set baseline flux = 0
  trojans.append(t_flux)
  # print(f"t0 = {t_params.t0}d = {params.t0} - {60*dpd} trojan shift + {t0_scatter} scatter")
  t_b = (np.sin(t_params.inc)*t_params.a) / rstar     # impact parameter
  print(f"trojan {trojan} inclination: {t_params.inc}, b = {t_b}")
  plt.plot(t,t_flux)

for greek in range(n_greeks):
  t0_scatter = np.random.uniform(-1,1,1)[0]*20*dpd
  #random scatter in inclination
  inc_scatter = np.random.uniform(-1,1,1)[0]*13.7 # 13.7 deg is mean for Jupiter trojans 
  
  g_params = batman.TransitParams()
  # t0 = t0_planet + greek phase +/- random scatter
  g_params.t0 = params.t0 + 60*dpd + t0_scatter
  g_params.per = 100.
  g_params.rp = rarstar * rstar # * 1e-5 is 10 km
  g_params.a = 3.   # 3 rstar
  g_params.inc = 90. + inc_scatter
  g_params.ecc = 0.
  g_params.w = 90.
  g_params.u = [0.1,0.3]            #limb darkening coefficients [u1, u2]
  g_params.limb_dark = "quadratic"       #limb darkening model
  g_m = batman.TransitModel(g_params, t)
  g_flux = g_m.light_curve(g_params) - 1.  # set baseline flux = 0
  greeks.append(g_flux)
  # print(f"t0 = {t_params.t0}d = {params.t0} - {60*dpd} trojan shift + {t0_scatter} scatter")
  g_b = (np.sin(g_params.inc)*g_params.a) / rstar     # impact parameter
  print(f"greek {greek} inclination: {g_params.inc}, b = {g_b}")
  plt.plot(t,g_flux)

asteroids = np.append(greeks,trojans,axis=0)
everything = np.append(asteroids,[plan_flux],axis=0)
lc_everything = np.sum(everything, axis=0)
  
#for i in np.arange(len(trojans)):
plt.plot(t,lc_everything)
plt.ticklabel_format(style='sci',scilimits=[-1,1]) 
plt.ylim(-100e-6,10e-6)
plt.xlim(50,150)
plt.xlabel("Time [d]")
plt.ylabel("Relative flux")
plt.text(100,0,f"R_asteroid / Rstar = {rarstar}")
plt.show()


#---------------------------------------------------------------
#             Set up asteroid ensemble parameters
#               Greeks: leading
#               Trojans: trailing
#---------------------------------------------------------------


"""
# pulled from tango.py
if is_plot_model:
    from pytransit import QuadraticModel
    #Let us use PyTransit to compute the transits
    xtr_model = np.arange(min(tvec)-size_time,max(tvec)+size_time,0.0001)
    fluxtr_model = [1]*len(xtr_model)
    for o in range(npl):
        tm = QuadraticModel(interpolate=False)
        tm.set_data(xtr_model,exptimes=t_cad,nsamples=n_cad)
        fluxtr_planet = tm.evaluate(k=rp[o], ldc=[u1,u2], t0=T0[o], p=P[o], a=a[o], i=inclination[o])
        #Avoid errors because of occultations calculated by pytransit
        phase = abs(((xtr_model-T0[o])%P[o])/P[o])
        phase[phase>0.5] -= 1
        fluxtr_planet[abs(phase)>0.125] = 1.
        fluxtr_model *= fluxtr_planet
    #Change the model to percentage
    fluxtr_model *= 100
"""
