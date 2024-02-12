# Import respiration data
import pandas as pd

# Import experimental respiration data and list of site IDs to be used for model parameterization
df = pd.read_csv('../../data/incubations/multiplexer/processed-respiration-data.csv')
sites = pd.read_csv('data/site-ID-to-use-for-calibrations.csv')
IDs = sites['site']
df=df[df["site.id"].isin(IDs)]


# Plot experimental respiration data from Jack pine sites following the 'No burn' treatment. 
from matplotlib import pyplot

df=df[df['vegetation']=='Pinus_banksiana']

fig,ax=pyplot.subplots(nrows=1)
no_burn = df[df["burn.trtmt"]=='control']
x_noburn = no_burn['whole.days.since.wet.up']
y_noburn = no_burn['new.cum_CO2C_g']
ax.scatter(x_noburn, y_noburn, s=3, color = 'grey', label = "no burn")
ax.set_xlabel('Time (days)')
ax.set_ylabel('Cumulative C-CO2 g respired')
ax.set_title('No burn', y=1, pad=-15)


# fig,ax=pyplot.subplots(nrows=3)
# no_burn = df[df["burn.trtmt"]=='control']
# low_severity = df[df["burn.trtmt"]=='low severity']
# high_severity = df[df["burn.trtmt"]=='high severity']
# x_noburn = no_burn['whole.days.since.wet.up']
# y_noburn = no_burn['frac.C']
# x_low = low_severity['whole.days.since.wet.up']
# y_low = low_severity['frac.C']
# x_burn = high_severity['whole.days.since.wet.up']
# y_burn = high_severity['frac.C']
# ax[0].scatter(x_noburn, y_noburn, s=3, color = 'black', label = "no burn")
# ax[1].scatter(x_burn, y_burn, s=4, color = 'grey', label = "burn")
# ax[2].scatter(x_low, y_low, s=4, color = 'red', label = "burn")
# ax[0].set_xlabel('Time (days)')
# ax[1].set_xlabel('Time (days)')
# ax[2].set_xlabel('Time (days)')
# ax[0].set_ylabel('Fraction of initial\nC remaining')
# ax[1].set_ylabel('Fraction of initial\nC remaining')
# ax[2].set_ylabel('Fraction of initial\nC remaining')
# ax[0].set_title('No burn', y=1, pad=-15)
# ax[1].set_title('Low severity burned soil', y=1, pad=-15)
# ax[2].set_title('High severity burned soil', y=1, pad=-15)

pyplot.show()


# Import CORPSE simulation results.

# This section plots the results
# Each set of results should have the same set of pools as the initial values structure from the beginning of the simulation
import Whitman_sims
import CORPSE_array

fig,ax=pyplot.subplots(nrows=2,ncols=1,clear=True,num='CORPSE results')

for sim in Whitman_sims.results:
    totalC=CORPSE_array.sumCtypes(Whitman_sims.results[sim][0], 'u')+CORPSE_array.sumCtypes(Whitman_sims.results[sim][0], 'p')
    ax[0].plot(Whitman_sims.t*365,Whitman_sims.results[sim][0]['CO2'].diff()/totalC[0]*100,label=sim)

    # ax[1].plot(t*365,results[sim][0]['uFastC'],label='Simple')
    # ax[1].plot(t*365,results[sim][0]['uSlowC'],label='Complex')
    # ax[1].plot(t*365,results[sim][0]['uNecroC'],label='Necromass')

    ax[1].plot(Whitman_sims.t*365,Whitman_sims.results[sim][0]['livingMicrobeC']/totalC[0]*100)
    #ax[2].plot(t*365,totalC) # How to plot total C in all four functional groups combined?

ax[0].set_xlabel('Time (days)')
ax[1].set_xlabel('Time (days)')
# ax[2].set_xlabel('Time (days)')
ax[0].set_ylabel('CO$_2$ flux rate \n(% initial C/day)')
ax[1].set_ylabel('Microbial biomass \n(% initial C)')
# ax[1].set_ylabel('SOM pools')
# ax[1].legend()
ax[0].legend(fontsize='small')
ax[0].set_title('CO$_2$ fluxes', y=1, pad=-15)
ax[1].set_title('Microbial biomass', y=1, pad=-15)

pyplot.subplots_adjust(hspace = 0.2)

pyplot.show()


# Put the model and experimental results on the same figure.
fig,ax=pyplot.subplot(nrows=1)

# But i don't want to have the graphs beside each other. i want them on top of each other?



