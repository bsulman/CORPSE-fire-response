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
y_noburn = no_burn['g.CO2C.per.initial.g.C.per.day']*100
ax.scatter(x_noburn, y_noburn, s=3, color = 'grey', label = "no burn")
ax.set_xlabel('Time (days)')
ax.set_ylabel('CO$_2$ flux rate \n(% initial C/day)')
ax.set_title('Laboratory incubation results', y=1, pad=-15)
ax.legend(loc = 'center right')


pyplot.show()


# Import CORPSE simulation results.

# This section plots the results
# Each set of results should have the same set of pools as the initial values structure from the beginning of the simulation
import Whitman_sims
import CORPSE_array

fig,ax=pyplot.subplots(nrows=2,ncols=1,clear=True,num='CORPSE results')

for sim in Whitman_sims.results:
    totalC=CORPSE_array.sumCtypes(Whitman_sims.results[sim][0], 'u')+CORPSE_array.sumCtypes(Whitman_sims.results[sim][0], 'p')
    ax[0].plot(Whitman_sims.t*365,Whitman_sims.results[sim][0]['CO2'].diff()/totalC[0]*100,label=sim) # % of initial C respired each day
#   ax[0].plot(Whitman_sims.t*365,Whitman_sims.results[sim][0]['CO2']/totalC[0]*100, label = sim) # Cumulative % of initial C respired
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


fig,ax=pyplot.subplots(nrows=1,ncols=1,clear=True,num='CORPSE results')
for sim in Whitman_sims.results:
    totalC=CORPSE_array.sumCtypes(Whitman_sims.results[sim][0], 'u')+CORPSE_array.sumCtypes(Whitman_sims.results[sim][0], 'p')
    ax.plot(Whitman_sims.t*365, Whitman_sims.results[sim][0]['uFastC']/totalC[0]*100, color = 'green', label = 'uFastC')
    #ax.plot(Whitman_sims.t*365, Whitman_sims.results[sim][0]['uSlowC']/totalC[0]*100, color = 'blue', label = 'uSlowC')
    ax.plot(Whitman_sims.t*365, Whitman_sims.results[sim][0]['uNecroC']/totalC[0]*100, color = 'orange', label = 'uNecroC')
    ax.plot(Whitman_sims.t*365, Whitman_sims.results[sim][0]['uPyC']/totalC[0]*100, color = 'purple', label = 'uPyC')
ax.set_xlabel('Time (days)')
ax.set_ylabel('Percent of total C')
ax.legend(fontsize='small')
ax.set_title('Percent of total C in various pools', y=1, pad=-15)

pyplot.subplots_adjust(hspace = 0.2)

pyplot.show()



# Put the model and experimental results on the same figure.
fig,ax=pyplot.subplots(nrows=1)

for sim in Whitman_sims.results:
    totalC=CORPSE_array.sumCtypes(Whitman_sims.results[sim][0], 'u')+CORPSE_array.sumCtypes(Whitman_sims.results[sim][0], 'p')
    ax.plot(Whitman_sims.t*365,Whitman_sims.results[sim][0]['CO2']/totalC[0]*100, color = 'black', label = 'Model results: '+sim) # Cumulative % of initial C respired

no_burn = df[df["burn.trtmt"]=='control']
x_noburn = no_burn['whole.days.since.wet.up']
y_noburn = no_burn['new.cum_CO2C_g']/no_burn['initial.total.C.g']*100
ax.scatter(x_noburn, y_noburn, s=3, color = 'dimgrey', label = "Lab results: No burn")

ax.set_xlabel('Time (days)')
ax.set_ylabel('Cumulative C-CO$_2$ respired\n(% of initial C)')
#ax.set_title('Laboratory incubation results', y=1, pad=-15)
ax.legend(loc = 'upper left')

pyplot.show()
   
    

# Save results file
df_output = Whitman_sims.results['No burn'][0]
df_output.to_csv('No_burn.csv', index_label='Time')

