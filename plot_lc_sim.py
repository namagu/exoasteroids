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
import astropy.units as u

debug = False 

# NEEDS TO MATCH PLANET PERIOD PERIODT
dpd = 100 / 360.  # days per degree

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

#a,b,c = ax_t0.hist(np.abs(t_t0_choices),30,density=True)
#ax_t0.set_title("|t0 phase shift| distribution [days]")
#d,e,f = ax_inc.hist(inc_choices,30,density=True)
#ax_inc.set_title("inclination distribution [deg]")
#ax_inc.vlines(params.inc,0,0.05,linestyle='dotted')

# load in data
print("loading in data")
lc_data = np.load('lc_flux.npz')
g_b_data = np.load('greek_array.npz')
t_b_data = np.load('trojan_array.npz')

# set up variables
t = np.arange(0, 1000, 0.0001)
all_g_b = g_b_data['all_g_b']
all_t_b = t_b_data['all_t_b']
lc_everything = np.sum(lc_data['everything'],axis=0)

#ax.plot(t,t_flux,label=f"t {trojan} b = {t_b:0.3f}")
#ax.plot(t,g_flux,label=f"g {greek} b = {g_b:0.2f}")

print("plotting")
i,j,k = ax_b.hist(np.append(all_g_b, all_t_b, axis=0),30,density=True)
ax_b.set_title("impact parameter distribution")
ax_b.vlines(1,0,1,linestyle='dotted')

ax.plot(t,lc_everything,label="summed lc")

if debug:
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

plt.show()
