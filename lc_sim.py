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
from astropy import units as u
from astropy import constants as const

#---------------------------------------------------------------
#             Set up planet parameters
#---------------------------------------------------------------
print("Setting up star, planet model")

# TODO: be fancy and pull from tango input.py
rstar = const.R_sun
mstar = const.M_sun # solar mass

# Keep P, a/Rstar, e, omega constant, vary t0, Rp/Rstar, inc.; Pjup = 100d

def get_a(p):
  """
  Calculate semi-major axis
  
  p: period [days], Quantity
  
  Returns a in AU, Quantity
  
  """
  G = const.G # AU^3 Msun^-1 days^-2
  a = ((G * (mstar + 0.1*mstar) * (p.to(u.s)**2)) / (4 * (np.pi)**2))**(1/3)
  return a.to(u.au)

def get_b(inc, a):
  """
  Calculate impact parameter
  
  inc: inclination [deg], Quantity
  a: semi-major axis [AU], Quantity 
  
  Returns b, dimensionless Quantity
  """
  b = (np.cos(np.deg2rad(inc))*a) / rstar
  return b.decompose()

def get_tdur(p, rp, a, inc):
  '''
  Calculate transit duration (in days)
  
  p: period [days]
  rp: planet radius [R_sun]
  a: semi-major axis [AU]
  inc: inclination [deg]
  
  '''
  dur = ((p/np.pi)*(np.sqrt((rstar + rp)**2 - (get_b(inc,a)*rstar)**2) / a))
  return dur.decompose().to(u.hr) 

# FYI 5 AU = 1120 Rsun

def get_quantity_value(q):
  return(q.value)


params = batman.TransitParams()
params.t0 = get_quantity_value(0. * u.day )             #time of inferior conjunction
params.per = get_quantity_value(100 * u.day)            #orbital period (Kepler-10 c: 45.2942679 d)
params.rp = get_quantity_value((0.1 * rstar).to(u.R_sun))             #planet radius (in units of stellar radii)
params.a = get_quantity_value(get_a(params.per * u.day).to(u.R_sun))        #semi-major axis (in units of stellar radii)
params.inc = 90 # get_quantity_value((90 - 1.305) * u.deg)   #orbital inclination (in degrees)
params.ecc = 0.                     #eccentricity
params.w = get_quantity_value(90. * u.deg)              #longitude of periastron (in degrees)
params.u = [0.1, 0.3]               #limb darkening coefficients [u1, u2]
params.limb_dark = "quadratic"      #limb darkening model

t = np.arange(0,1000, 0.0001) # * u.day
# make a data point every 2.4 hours

m = batman.TransitModel(params, t)      # initializes model
plan_flux = m.light_curve(params) - 1.  # calculates light curve, start from 0 baseline

p_tdur = get_tdur(params.per * u.day, params.rp * u.R_sun, params.a * u.R_sun, params.inc * u.deg)

print(f"planet a: {get_a(params.per * u.day):0.5}, planet b: {get_b(params.inc * u.deg, params.a * u.R_sun):0.5}")
print(f"Planet T_dur should be {p_tdur:0.5}")

#---------------------------------------------------------------
#             Set up asteroid ensemble parameters
#               Greeks: leading
#               Trojans: trailing
#---------------------------------------------------------------

n_trojans = 100  #100     # this is a starting point
n_greeks = 100  #100

print(f"Setting up asteroid model for {n_trojans} trailing, {n_greeks} leading")

dpd = params.per / 360.  # days per degree

rarstar = (0.27 * const.R_earth)/const.R_sun  # r_asteroid [Rsun] / r_star [Rsun] (r_moon = 0.27 * r_earth)

# roughly 25 deg spread in phase
t0_scatter = np.random.uniform((-12.5)*dpd, (12.5 )*dpd,1000) # days
t_t0_choices = params.t0 - (60)*dpd + t0_scatter
g_t0_choices = params.t0 + (60)*dpd + t0_scatter

#random scatter in inclination
inc_choices = np.random.normal(params.inc,10,1000) # FYI 13.7 deg is mean for Jupiter trojans 

trojans = []  
greeks = []

all_t_b = []
all_g_b = []

debug = False   # if True, prints inc, b for each asteroid and adds them to lc plot

t_transits = 0
g_transits = 0

print("Creating trojans")
for trojan in range(n_trojans):
  t_params = batman.TransitParams()
  # t0 = t0_planet - trojan phase +/- random scatter
  t_params.t0 = np.random.choice(t_t0_choices,1,replace=False)[0]
  t_params.per = 100. 
  t_params.rp = get_quantity_value(rarstar * rstar.to(u.R_sun))
  t_params.a = params.a   
  t_params.inc = np.random.choice(inc_choices,1,replace=False)[0]
  #t_params.inc = params.inc - np.random.choice(inc_scatter,1,replace=False)[0]
  t_params.ecc = params.ecc
  t_params.w = params.w
  t_params.u = [0.1,0.3]            #limb darkening coefficients [u1, u2]
  t_params.limb_dark = "quadratic"       #limb darkening model
  t_m = batman.TransitModel(t_params, t)
  
  t_b = get_b(t_params.inc * u.deg, t_params.a * u.R_sun) 
  all_t_b.append(t_b)     # save impact parameter
  
  if np.abs(t_b) < 1:
    t_flux = t_m.light_curve(t_params) - 1.  # set baseline flux = 0
    trojans.append(t_flux)
    t_transits +=1
  else:
    trojans.append(np.zeros(len(t)))
  
  if debug:
    print(f"t0 = {t_params.t0}d = {params.t0} - {60*dpd} trojan shift + {t0_scatter} scatter")
    print(f"trojan {trojan} inc: {t_params.inc}, a: {t_params.a}, b = {t_b}")

print(f"{t_transits}/{n_trojans} trailing asteroids have |b|<1")
print("Saving trojans")
np.savez('trojan_array',trojans = trojans, all_t_b = all_t_b, t_t0_choices = t_t0_choices, n_t_visible = [t_transits, n_trojans])

print("Creating greeks")
for greek in range(n_greeks):
  g_params = batman.TransitParams()
  # t0 = t0_planet + greek phase +/- random scatter
  g_params.t0 = np.random.choice(g_t0_choices,1,replace=False)[0]
  g_params.per = 100.
  g_params.rp = get_quantity_value(rarstar * rstar.to(u.R_sun))
  g_params.a = params.a  
  g_params.inc = np.random.choice(inc_choices,1,replace=False)[0]
  g_params.ecc = params.ecc
  g_params.w = params.w
  g_params.u = [0.1,0.3]            #limb darkening coefficients [u1, u2]
  g_params.limb_dark = "quadratic"       #limb darkening model
  g_m = batman.TransitModel(g_params, t)
  
  g_b = get_b(g_params.inc * u.deg, g_params.a * u.R_sun) # impact parameter
  all_g_b.append(g_b)
  
  if np.abs(g_b) < 1:
    g_flux = g_m.light_curve(g_params) - 1.  # set baseline flux = 0
    greeks.append(g_flux)
    g_transits +=1
  else:
    greeks.append(np.zeros(len(t)))  

  if debug:
    print(f"t0 = {t_params.t0}d = {params.t0} - {60*dpd} trojan shift + {t0_scatter} scatter")
    print(f"greek {greek} inc: {g_params.inc}, a: {g_params.a}, b = {g_b}")

print(f"{g_transits}/{n_greeks} leading asteroids have |b|<1")

print("Saving greeks")
np.savez('greek_array',greeks = greeks, all_g_b = all_g_b, g_t0_choices = g_t0_choices, n_g_visible = [g_transits, n_greeks])

print("Putting it all together")
asteroids = np.append(greeks,trojans,axis=0)
everything = np.append(asteroids,[plan_flux],axis=0)
lc_everything = np.sum(everything, axis=0)

np.savez('lc_flux',t=t,plan_flux=plan_flux,everything=lc_everything)

np.savez('inc_choices',inc_choices=inc_choices)

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
