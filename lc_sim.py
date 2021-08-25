#--------------------------------------------------------------
#                        exoasteroids
#               Model exoplanet and exo-asteroid transits!
#                     Aug 2021, Natalia Guerrero
#---------------------------------------------------------------

import batman
import numpy as np
import matplotlib.pyplot as plt
import sys


#---------------------------------------------------------------
#             Set up planet parameters
#---------------------------------------------------------------

# TODO: be fancy and pull from tango input.py
rstar = 1.056

params = batman.TransitParams()
params.t0 = 138.6776401              #time of inferior conjunction
params.per = 45.2942679              #orbital period
params.rp = 0.0219198 * rstar        #planet radius (in units of stellar radii)
params.a = 34.82371319132805 * rstar #semi-major axis (in units of stellar radii)
params.inc = 88.7850952              #orbital inclination (in degrees)
params.ecc = 0.                      #eccentricity
params.w = 90.                       #longitude of periastron (in degrees)
params.u = [0.455, 0.221]            #limb darkening coefficients [u1, u2]
params.limb_dark = "quadratic"       #limb darkening model

t = np.arange(131.5123148808634, 1591.0016788740177, 0.0001)
# make a data point every 8.64 seconds

m = batman.TransitModel(params, t)    #initializes model
plan_flux = m.light_curve(params)          #calculates light curve


#---------------------------------------------------------------
#             Set up asteroid ensemble parameters
#               Greeks: leading
#               Trojans: trailing
#---------------------------------------------------------------

n_trojans = 1 #100     # this is a starting point
n_greeks = 100

dpd = params.per / 360.  # days per degree


for trojan in range(trojans):
  t_params = batman.TransitParams()
  t_params.t0 = params.t0 - 60*


plt.plot(t, flux)
plt.xlabel("Time from central transit")
plt.ylabel("Relative flux")
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
