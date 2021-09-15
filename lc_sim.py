#--------------------------------------------------------------
#                        exoasteroids
#               Model exoplanet and exo-asteroid transits!
#                     Aug 2021, Natalia Guerrero
#---------------------------------------------------------------

import batman
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import sys
import astropy.units as u

#---------------------------------------------------------------
#             Set up planet parameters
#---------------------------------------------------------------

# TODO: be fancy and pull from tango input.py
rstar = 1.

# Keep P, a/Rstar, e, omega constant, vary t0, Rp/Rstar, inc.; Pjup = 100d

# FYI 5 AU = 1120 Rsun

params = batman.TransitParams()
params.t0 = 0.                      #time of inferior conjunction
params.per = 100                     #orbital period (Kepler-10 c: 45.2942679 d)
params.rp = 0.1 * rstar             #planet radius (in units of stellar radii)
params.a = 10                        #semi-major axis (in units of stellar radii)
params.inc = 90 - 1.305              #orbital inclination (in degrees)
params.ecc = 0.                      #eccentricity
params.w = 90.                       #longitude of periastron (in degrees)
params.u = [0.1, 0.3]            #limb darkening coefficients [u1, u2]
params.limb_dark = "quadratic"       #limb darkening model

t = np.arange(0, 1000, 0.0001)
# make a data point every 2.4 hours

m = batman.TransitModel(params, t)      # initializes model
plan_flux = m.light_curve(params) - 1.  # calculates light curve, start from 0 baseline


#---------------------------------------------------------------
#             Set up asteroid ensemble parameters
#               Greeks: leading
#               Trojans: trailing
#---------------------------------------------------------------

fig = plt.figure(figsize=(10,8))
gs = gridspec.GridSpec(2,2)
# axs = plt.subplots(2,2,sharey=False, sharex=False,gridspec_kw={'wspace':0.1, 'hspace':0.1}, figsize=(10,8))

ax = fig.add_subplot(gs[1,:])
ax_t0 = fig.add_subplot(gs[0,0])
ax_inc = fig.add_subplot(gs[0,1])

n_trojans = 3  #100     # this is a starting point
n_greeks = 3  #100

dpd = float(params.per) / 360.  # days per degree

rarstar = 0.0024973   # r_asteroid [Rsun] / r_star [Rsun] (r_moon = 0.27 * r_earth)

# roughly 25 deg spread in phase
t0_scatter = np.random.uniform(-12.5*dpd,12.5*dpd,1000)  # days
t_t0_choices = params.t0 - 60*dpd + t0_scatter
g_t0_choices = params.t0 + 60*dpd + t0_scatter

#random scatter in inclination
inc_choices = np.random.normal(params.inc,10,1000) # FYI 13.7 deg is mean for Jupiter trojans 

a,b,c = ax_t0.hist(np.abs(t_t0_choices),30,density=True)
ax_t0.set_title("|t0 phase shift| distribution [days]")
d,e,f = ax_inc.hist(inc_choices,30,density=True)
ax_inc.set_title("inclination distribution [deg]")
ax_inc.vlines(params.inc,0,0.05,linestyle='dotted')

trojans = []  
greeks = []

for trojan in range(n_trojans):
  t_params = batman.TransitParams()
  # t0 = t0_planet - trojan phase +/- random scatter
  t_params.t0 = np.random.choice(t_t0_choices,1,replace=False)[0]
  
  t_params.per = 100.
  t_params.rp = rarstar * rstar 
  t_params.a = params.a   
  t_params.inc = np.random.choice(inc_choices,1,replace=False)[0]
  #t_params.inc = params.inc - np.random.choice(inc_scatter,1,replace=False)[0]
  t_params.ecc = params.ecc
  t_params.w = params.w
  t_params.u = [0.1,0.3]            #limb darkening coefficients [u1, u2]
  t_params.limb_dark = "quadratic"       #limb darkening model
  t_m = batman.TransitModel(t_params, t)
  t_flux = t_m.light_curve(t_params) - 1.  # set baseline flux = 0
  trojans.append(t_flux)
  # print(f"t0 = {t_params.t0}d = {params.t0} - {60*dpd} trojan shift + {t0_scatter} scatter")
  #t_b = (np.sin(t_params.inc)*t_params.a) / rstar     # impact parameter
  t_b = (np.cos(np.deg2rad(t_params.inc))*t_params.a) / rstar
  print(f"trojan {trojan} inc: {t_params.inc}, a: {t_params.a}, b = {t_b}")
  ax.plot(t,t_flux,label=f"t {trojan} b = {t_b:0.3f}")

for greek in range(n_greeks):
  
  g_params = batman.TransitParams()
  # t0 = t0_planet + greek phase +/- random scatter
  g_params.t0 = np.random.choice(g_t0_choices,1,replace=False)[0]
  g_params.per = 100.
  g_params.rp = rarstar * rstar # * 1e-5 is 10 km
  g_params.a = params.a  
  g_params.inc = np.random.choice(inc_choices,1,replace=False)[0]
  g_params.ecc = params.ecc
  g_params.w = params.w
  g_params.u = [0.1,0.3]            #limb darkening coefficients [u1, u2]
  g_params.limb_dark = "quadratic"       #limb darkening model
  g_m = batman.TransitModel(g_params, t)
  g_flux = g_m.light_curve(g_params) - 1.  # set baseline flux = 0
  greeks.append(g_flux)
  # print(f"t0 = {t_params.t0}d = {params.t0} - {60*dpd} trojan shift + {t0_scatter} scatter")
  #g_b = (np.sin(g_params.inc)*g_params.a) / rstar     # impact parameter
  g_b = (np.cos(np.deg2rad(g_params.inc))*g_params.a) / rstar
  print(f"greek {greek} inc: {g_params.inc}, a: {g_params.a}, b = {g_b}")
  ax.plot(t,g_flux,label=f"g {greek} b = {g_b:0.2f}")

asteroids = np.append(greeks,trojans,axis=0)
everything = np.append(asteroids,[plan_flux],axis=0)
lc_everything = np.sum(everything, axis=0)

#plt.plot(t,lc_everything,label="summed lc")
ax.plot(t,plan_flux,label="planet")
ax.ticklabel_format(style='sci',scilimits=[-1,1]) 
# transit limits
ax.vlines([100 - 60*dpd-12.5*dpd, 100 - 60*dpd+12.5*dpd],-5e-5,1e-5,linestyle='dotted',alpha=0.5)
ax.vlines([100 + 60*dpd-12.5*dpd, 100 + 60*dpd+12.5*dpd],-5e-5,1e-5,linestyle='dotted',alpha=0.5)
# L4, L5 crossings
ax.vlines([100 - 60*dpd, 100 + 60*dpd],-5e-5,1e-5,linestyle='dotted',alpha=0.9)
ax.text(110, -1e-5,"L4 transit")
ax.text(80, -1e-5, "L5 transit")

ax.set_ylim(-50e-6,10e-6)
ax.set_xlim(50,150) #150)
ax.set_xlabel("Time [d]")
ax.set_ylabel("Relative flux")
ax.text(100,0.5e-5,f"R_asteroid / Rstar = {rarstar}")
ax.legend(loc="lower left")

#ax1.hist()

plt.show()


#---------------------------------------------------------------
#             fsklafkjlsdaflksjaflks
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
