#--------------------------------------------------------------
#                        plot_lc_sim.py
#               Plot the models we made
#                     Sep *16* 2021, Natalia Guerrero
#---------------------------------------------------------------

import batman
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import sys
from astropy import units as u
from astropy import constants as const
# import lc_sim

debug = False 

# NEEDS TO MATCH PLANET PERIOD PERIODT
dpd = 100 / 360.  # days per degree

# Other parameters
rarstar = (0.27 * const.R_earth)/const.R_sun   # r_asteroid [Rsun] / r_star [Rsun] (r_moon = 0.27 * r_earth)
pl_inc = 90. # 90 - 1.305

# 3 diagnostic histograms on top, 
# 1 lc window on botom
print("setting up plot")
fig = plt.figure(figsize=(10,8))
gs = gridspec.GridSpec(2,3)
# axs = plt.subplots(2,2,sharey=False, sharex=False,gridspec_kw={'wspace':0.1, 'hspace':0.1}, figsize=(10,8))

ax = fig.add_subplot(gs[1,:])
ax_t0 = fig.add_subplot(gs[0,0])
ax_inc = fig.add_subplot(gs[0,1])
ax_b = fig.add_subplot(gs[0,2])

# load in data
print("loading in data")
lc_data = np.load('lc_flux.npz')
g_data = np.load('greek_array.npz')
t_data = np.load('trojan_array.npz')
inc_data = np.load('inc_choices.npz')

# set up variables
t = np.arange(0, 1000, 0.0001)
#all_g_b = g_data['all_g_b']
all_t_b = t_data['all_t_b']
t_t0_choices = t_data['t_t0_choices']
sim_g_inc = inc_data['sim_g_inc']
print(f"len(sim_g_inc) = {len(sim_g_inc)}")
#inc_choices = inc_data['inc_choices']
lc_everything = lc_data['everything']

print("plotting")

# t0 histogram, just for trojans, but representative of greeks also
a,b,c = ax_t0.hist(np.abs(t_t0_choices),30)
ax_t0.set_title("|t0 phase shift| distribution [days]")

# impact parameters
i,j,k = ax_b.hist(all_t_b,50) #np.append(all_g_b, all_t_b, axis=0)
ax_b.set_title("impact parameter distribution")
ax_b.vlines(1,0, 7000,linestyle='dotted')
ax_b.vlines(-1,0, 7000 ,linestyle='dotted')

# inclination histogram
d,e,f = ax_inc.hist(sim_g_inc,20) #inc_choices
ax_inc.set_title("inclination distribution [deg]")
ax_inc.vlines(pl_inc,0,1,linestyle='dotted')

#ax.plot(t,t_flux,label=f"t {trojan} b = {t_b:0.3f}")
#ax.plot(t,g_flux,label=f"g {greek} b = {g_b:0.2f}")



ax.plot(t,lc_everything,label="summed lc")

if debug:
  ax.plot(t,lc_data['plan_flux'],label="planet")

  ax.ticklabel_format(style='sci',scilimits=[-1,1]) 
# transit limits
ax.vlines([100 - 60*dpd-12.5*dpd, 100 - 60*dpd+12.5*dpd],-5e-5,1e-5,linestyle='dotted',alpha=0.5)
ax.vlines([100 + 60*dpd-12.5*dpd, 100 + 60*dpd+12.5*dpd],-5e-5,1e-5,linestyle='dotted',alpha=0.5)
# L4, L5 crossings
ax.vlines([100 - 60*dpd, 100 + 60*dpd],-5e-5,1e-5,linestyle='dotted',alpha=0.9)
ax.text(110, -1e-5,"L4 transit")
ax.text(80, -1e-5, "L5 transit")

# how many of each transited
#ax.text(80, -6e-5, f"{t_data['n_t_visible'][0]}/{t_data['n_t_visible'][1]} trailing asteroids have |b|<1")
#ax.text(105, -6e-5, f"{g_data['n_g_visible'][0]}/{g_data['n_g_visible'][1]} leading asteroids have |b|<1")

ax.set_ylim(lc_everything.min() * 0.01, 10e-6)
ax.set_xlim(100 - 25, 100 + 25) #150)
ax.set_xlabel("Time [d]")
ax.set_ylabel("Relative flux")
ax.text(100,0.5e-5,f"R_asteroid / Rstar = {rarstar:0.5}")
ax.legend(loc="lower left")

plt.show()
