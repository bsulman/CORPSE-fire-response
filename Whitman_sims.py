import pandas
import CORPSE_array
import CORPSE_solvers
from numpy import array,arange
import copy

# Initial pools for soil organic matter simulation
# There are three kinds of chemically-defined C (Fast, slow, and microbial necromass). "Fast" has higher maximum decomposition rate and microbial CUE
# Each C type can be in a protected or unprotected state. When protected, it is not subject to microbial decomposition
SOM_init={'CO2': array(0.0),    # Cumulative CO2 from microbial respiration
 'livingMicrobeC': array(0.06), # Active, living microbial biomass
 'pFastC': array(1.97),         # Protected fast-decomposing C
 'pNecroC': array(22.1),        # Protected microbial necromass C
 'pSlowC': array(0.61),         # Protected slow-decomposing C
 'uFastC': array(5.0),          # Unprotected fast-decomposing C
 'uNecroC': array(0.19),        # Unprotected microbial necromass C
 'uSlowC': array(8.25)}         # Unprotected slow-decomposing C
 
# Parameters controlling the model
params={
    'vmaxref':{'Fast':9.0,'Slow':0.25,'Necro':4.5}, #  Relative maximum enzymatic decomp rates for each C type (year-1)
    'Ea':{'Fast':5e3,'Slow':30e3,'Necro':5e3},      # Activation energy (controls T dependence)
    'kC':{'Fast':0.01,'Slow':0.01,'Necro':0.01},    # Michaelis-Menton half saturation parameter (g microbial biomass/g substrate)
    'gas_diffusion_exp':0.6,  # Determines suppression of decomposition at high soil moisture
    'substrate_diffusion_exp':1.5,   # Controls suppression of decomp at low soil moisture
    'minMicrobeC':1e-3,       # Minimum microbial biomass (fraction of total C). Prevents microbes from going extinct and allows some slow decomposition under adverse conditions
    'Tmic':0.25,              # Microbial biomass mean lifetime (years)
    'et':0.6,                 # Fraction of microbial biomass turnover (death) that goes to necromass instead of being immediately mineralized to CO2
    'eup':{'Fast':0.6,'Slow':0.05,'Necro':0.6},     # Microbial carbon use efficiency for each substrate type (fast, slow, necromass). This amount is converted to biomass during decomposition, and the remainder is immediately respired as CO2
    'tProtected':75.0,        # Protected C turnover time (years). This is the time scale for which protected C is released back to unprotected state.
    'protection_rate':{'Fast':0.3,'Slow':0.001,'Necro':1.5}, # Protected carbon formation rate (year-1). Higher number means more will become protected. Can be modified as a function of soil texture/mineralogy to represent different sorption potentials
    'new_resp_units':True,   # At some point I changed the units of vmaxref to be normalized for other factors so they are actually in year-1 units. Leave this as True values that are easier to interpret.
}

# This makes an array of all the model time steps (in units of years). In this case, it starts at zero, ends at 120 days, and has a time step of one day
t=arange(0,120/365,1.0/365)

# This section is setting up different initial values and parameters for different simulations representing microbial community traits
# Here we set up an empty python dictionary to hold the different sets of parameters and initial values
initvals={}
paramsets={}

# Resistant: Low biomass loss, fast growth
initvals['Fast-growing survivor']=copy.deepcopy(SOM_init)      # Makes a copy of the default initial values. Need to use deepcopy so changing the value here doesn't change it for every simulation
initvals['Fast-growing survivor']['livingMicrobeC']=array(1.0) # Start with a high amount of initial microbial biomass. Assumes this simulation starts right after the fire
paramsets['Fast-growing survivor']=copy.deepcopy(params)
paramsets['Fast-growing survivor']['vmaxref']['Fast']=10.0     # Higher potential decomposition rate for more labile C



# Susceptible: High biomass loss, slow growth
initvals['Fire susceptible']=copy.deepcopy(SOM_init)
initvals['Fire susceptible']['livingMicrobeC']=array(0.05)    # Low initial microbial biomass
paramsets['Fire susceptible']=copy.deepcopy(params)
paramsets['Fire susceptible']['vmaxref']['Fast']=6.0          # Slower maximum decomposition rate for labile C



# Recovering: High biomass loss, fast growth
initvals['Post-fire rebounder']=copy.deepcopy(SOM_init)
initvals['Post-fire rebounder']['livingMicrobeC']=array(0.05) # Low initial microbial biomass
paramsets['Post-fire rebounder']=copy.deepcopy(params)
paramsets['Post-fire rebounder']['vmaxref']['Fast']=18.0      # Very fast potential decomposition/growth rate



# Resilient: Low biomass loss, slow growth
initvals['Slow-growing survivor']=copy.deepcopy(SOM_init)
initvals['Slow-growing survivor']['livingMicrobeC']=array(1.0) # High initial microbial biomass
paramsets['Slow-growing survivor']=copy.deepcopy(params)
paramsets['Slow-growing survivor']['vmaxref']['Fast']=6.0      # Slower decomposition/growth rate


# Set up a data structure to hold the results of the different simulations
results={}
# Goes through each functional type and runs a simulation using the appropriate set of parameters and initial values
# Simulations are assuming a constant temperature of 20 C and constant moisture of 60% of saturation
# Inputs are empty because this is running as an incubation without any constant inputs of C
for functype in initvals:
    results[functype] = CORPSE_solvers.run_models_ODE(Tmin=20.0,Tmax=20.0,thetamin=0.6,thetamax=0.6,
                                            times=t,inputs={},clay=20.0,initvals=initvals[functype],params=paramsets[functype])


# This section plots the results
# Each set of results should have the same set of pools as the initial values structure from the beginning of the simulation
from matplotlib import pyplot

fig,ax=pyplot.subplots(nrows=2,ncols=1,clear=True,num='CORPSE results')
for sim in results:
    totalC=CORPSE_array.sumCtypes(results[sim][0], 'u')+CORPSE_array.sumCtypes(results[sim][0], 'p')
    ax[0].plot(t*365,results[sim][0]['CO2'].diff()/totalC[0]*100,label=sim)

    # ax[1].plot(t*365,results[sim][0]['uFastC'],label='Simple')
    # ax[1].plot(t*365,results[sim][0]['uSlowC'],label='Complex')
    # ax[1].plot(t*365,results[sim][0]['uNecroC'],label='Necromass')

    ax[1].plot(t*365,results[sim][0]['livingMicrobeC']/totalC[0]*100)

ax[0].set_xlabel('Time (days)')
ax[1].set_xlabel('Time (days)')
# ax[2].set_xlabel('Time (days)')
ax[0].set_ylabel('CO$_2$ flux rate (% initial C/day)')
ax[1].set_ylabel('Microbial biomass (% initial C)')
# ax[1].set_ylabel('SOM pools')
# ax[1].legend()
ax[0].legend(fontsize='medium')
ax[0].set_title('CO$_2$ fluxes')
ax[1].set_title('Microbial biomass')

pyplot.show()

