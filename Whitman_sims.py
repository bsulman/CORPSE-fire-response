import pandas
import CORPSE_solvers
from numpy import array,arange
import copy


# Initial pools for soil organic matter simulation
# There are three kinds of chemically-defined C (Fast, slow, and microbial necromass). "Fast" has higher maximum decomposition rate and microbial CUE
# Each C type can be in a protected or unprotected state. When protected, it is not subject to microbial decomposition
SOM_init={'CO2': array(0.0),    # Cumulative CO2 from microbial respiration
 'livingMicrobeC': array(0.012), # Active, living microbial biomass
 'pFastC': array(0.18),         # Protected fast-decomposing C
 'pNecroC': array(0.2),        # Protected microbial necromass C
 'pSlowC': array(0.6),         # Protected slow-decomposing C
 'uFastC': array(.1),          # Unprotected fast-decomposing C
 'uNecroC': array(.1),        # Unprotected microbial necromass C
 'uSlowC': array(99.0)}         # Unprotected slow-decomposing C
 
# Parameters controlling the model
params={
    'vmaxref':{'Fast':19.0,'Slow':1.25,'Necro':10}, #  Relative maximum enzymatic decomp rates for each C type (year-1)
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
t=arange(0,90/365,1/365)

# This section is setting up different initial values and parameters for different simulations representing microbial community traits
# Here we set up an empty python dictionary to hold the different sets of parameters and initial values
initvals={}
paramsets={}

# Microbial community, no burn
initvals['No burn']=copy.deepcopy(SOM_init)      # Makes a copy of the default initial values. Need to use deepcopy so changing the value here doesn't change it for every simulation
paramsets['No burn']=copy.deepcopy(params)

# # Microbial community, low severity burn
# initvals['Low severity burn']=copy.deepcopy(SOM_init)      # Makes a copy of the default initial values. Need to use deepcopy so changing the value here doesn't change it for every simulation
# initvals['Low severity burn']['livingMicrobeC']=array(0.01) # Start with low initial microbial biomass. Assumes this simulation starts right after the fire
# paramsets['Low severity burn']=copy.deepcopy(params)
# paramsets['Low severity burn']['vmaxref']['Fast']=20.0     # Higher potential decomposition rate for more labile C
# paramsets['Low severity burn']['eup']['Fast'] = 0.2        # Lower CUE for ufastC and necromass
# paramsets['Low severity burn']['eup']['Necro'] = 0.2


# # Resistant: Low biomass loss, fast growth
# # Fast-growing survivor
# initvals['Fast-growing survivor']=copy.deepcopy(SOM_init)      # Makes a copy of the default initial values. Need to use deepcopy so changing the value here doesn't change it for every simulation
# initvals['Fast-growing survivor']['livingMicrobeC']=array(0.01) # Start with a high amount of initial microbial biomass. Assumes this simulation starts right after the fire
# paramsets['Fast-growing survivor']=copy.deepcopy(params)
# paramsets['Fast-growing survivor']['vmaxref']['Fast']=20.0     # Higher potential decomposition rate for more labile C
# paramsets['Fast-growing survivor']['eup']['Fast'] = 0.2
# paramsets['Fast-growing survivor']['eup']['Necro'] = 0.2


# # Susceptible: High biomass loss, slow growth
# # Fire susceptible
# initvals['Fire susceptible']=copy.deepcopy(SOM_init)
# initvals['Fire susceptible']['livingMicrobeC']=array(0.001)    # Low initial microbial biomass
# paramsets['Fire susceptible']=copy.deepcopy(params)
# paramsets['Fire susceptible']['vmaxref']['Fast']=0.2          # Slower maximum decomposition rate for labile C
# paramsets['Fire susceptible']['eup']['Fast'] = 0.5
# paramsets['Fire susceptible']['eup']['Necro'] = 0.5


# # Recovering: High biomass loss, fast growth
# # Post-fire rebounder
# initvals['Post-fire rebounder']=copy.deepcopy(SOM_init)
# initvals['Post-fire rebounder']['livingMicrobeC']=array(0.001) # Low initial microbial biomass
# paramsets['Post-fire rebounder']=copy.deepcopy(params)
# paramsets['Post-fire rebounder']['vmaxref']['Fast']=20.0      # Very fast potential decomposition/growth rate
# paramsets['Post-fire rebounder']['eup']['Fast'] = 0.2
# paramsets['Post-fire rebounder']['eup']['Necro'] = 0.2


# # Resilient: Low biomass loss, slow growth
# # Slow-growing survivor
# initvals['Slow-growing survivor']=copy.deepcopy(SOM_init)
# initvals['Slow-growing survivor']['livingMicrobeC']=array(0.01) # High initial microbial biomass
# paramsets['Slow-growing survivor']=copy.deepcopy(params)
# paramsets['Slow-growing survivor']['vmaxref']['Fast']=0.2      # Slower decomposition/growth rate
# paramsets['Slow-growing survivor']['eup']['Fast'] = 0.5
# paramsets['Slow-growing survivor']['eup']['Necro'] = 0.5


# Set up a data structure to hold the results of the different simulations
results={}
# Goes through each functional type and runs a simulation using the appropriate set of parameters and initial values
# Simulations are assuming a constant temperature of 20 C and constant moisture of 60% of saturation
# Inputs are empty because this is running as an incubation without any constant inputs of C
for functype in initvals:
    results[functype] = CORPSE_solvers.run_models_ODE(Tmin=5.0,Tmax=20.0,thetamin=0.4,thetamax=0.9,
                                            times=t,inputs={},clay=2.5,initvals=initvals[functype],params=paramsets[functype])
