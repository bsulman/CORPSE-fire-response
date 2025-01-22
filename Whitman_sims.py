import pandas as pd
import CORPSE_solvers
from numpy import array,arange
import copy


# Initial pools for soil organic matter simulation
# There are three kinds of chemically-defined C (Fast, slow, and microbial necromass). "Fast" has higher maximum decomposition rate and microbial CUE
# Each C type can be in a protected or unprotected state. When protected, it is not subject to microbial decomposition
SOM_init={'CO2': array(0.0),    # Cumulative C-CO2 from microbial respiration
 'MBC_1': array(4.5),           # Active, living microbial biomass
 'MBC_2': array(0.045),         # MBC_2 = 1% of MBC_1
 'MBC_3': array(0.00),
 'MBC_4': array(0.00), 
 'pFastC': array(0.0),          # Protected fast-decomposing C
 'pNecroC': array(0.0),         # Protected microbial necromass C
 'pSlowC': array(0.0),          # Protected slow-decomposing C
 'pPyC': array(0.0),            # Protected PyC
 'uFastC': array(3.0),          # Unprotected fast-decomposing C, Changed value from 10 to 11 based on calibration results
 'uNecroC': array(3.0),         # Unprotected microbial necromass C, 
 'uSlowC': array(88),           # Unprotected slow-decomposing C, Changed value from 86 to 60.2 based on calibration results
 'uPyC': array(4.0)}            # Unprotected PyC

 
# Parameters controlling the model (sandy soil, no burn)
params={
    'vmaxref':{'MBC_1':{'Fast':6.9,'Slow':0.11,'Necro':7.0, 'Py':0.1},      #  Relative maximum enzymatic decomp rates for each C type (year-1)
               'MBC_2':{'Fast':19.2,'Slow':0.0064,'Necro':45.0, 'Py':0.01},
               'MBC_3':{'Fast':0.0,'Slow':0.0,'Necro':0.0, 'Py':0.0}, 
               'MBC_4':{'Fast':0.0,'Slow':0.0,'Necro':0.0, 'Py':0.0}},
    'Ea':{'Fast':5e3,'Slow':30e3,'Necro':5e3, 'Py':35e3},      # Activation energy (controls T dependence)
    'kC':{'MBC_1':{'Fast':0.01,'Slow':0.01,'Necro':0.01, 'Py':0.01},    # Michaelis-Menton half saturation parameter (g microbial biomass/g substrate)
          'MBC_2':{'Fast':0.01,'Slow':0.04,'Necro':0.01, 'Py':0.04},    # Michaelis-Menton half saturation parameter (g microbial biomass/g substrate)
          'MBC_3':{'Fast':0.1,'Slow':0.01,'Necro':0.1, 'Py':0.01},    # Michaelis-Menton half saturation parameter (g microbial biomass/g substrate)
          'MBC_4':{'Fast':0.1,'Slow':0.01,'Necro':0.1, 'Py':0.01}},    # Michaelis-Menton half saturation parameter (g microbial biomass/g substrate)
    'gas_diffusion_exp':0.6,  # Determines suppression of decomposition at high soil moisture
    'substrate_diffusion_exp':1.5,   # Controls suppression of decomp at low soil moisture
    'minMicrobeC':{'MBC_1':1e-3,'MBC_2':1e-5,'MBC_3':1e-3,'MBC_4':1e-3}, # Minimum microbial biomass (fraction of total C). Prevents microbes from going extinct and allows some slow decomposition under adverse conditions
    'Tmic':{'MBC_1':0.5,'MBC_2':0.15,'MBC_3':0.25,'MBC_4':0.25}, # Microbial biomass mean lifetime (years)
    'et':{'MBC_1':0.8,'MBC_2':0.8,'MBC_3':0.6,'MBC_4':0.6},                  # Fraction of microbial biomass turnover (death) that goes to necromass instead of being immediately mineralized to CO2
    'eup':{'MBC_1':{'Fast':0.6,'Slow':0.3,'Necro':0.65, 'Py':0.15},      # Microbial carbon use efficiency for each substrate type (fast, slow, necromass). This amount is converted to biomass during decomposition, and the remainder is immediately respired as CO2
            'MBC_2':{'Fast':0.36,'Slow':0.1,'Necro':0.3, 'Py':0.05},     
            'MBC_3':{'Fast':0.5,'Slow':0.3,'Necro':0.6, 'Py':0.1},
            'MBC_4':{'Fast':0.5,'Slow':0.3,'Necro':0.6, 'Py':0.1}},
    'tProtected':75.0,        # Protected C turnover time (years). This is the time scale for which protected C is released back to unprotected state.
    'protection_rate':{'Fast':0.0,'Slow':0.001,'Necro':0, 'Py':0}, # Protected carbon formation rate (year-1). Higher number means more will become protected. Can be modified as a function of soil texture/mineralogy to represent different sorption potentials
    'new_resp_units':True,   # At some point I changed the units of vmaxref to be normalized for other factors so they are actually in year-1 units. Leave this as True values that are easier to interpret.
}


envir_vals={'thetamin': array(0.5),
              'thetamax': array(0.7),
              'porosity': array(0.4)}



import CORPSE_array
CORPSE_array.check_params(params)

# This makes an array of all the model time steps (in units of years). In this case, it starts at zero, ends at 120 days, and has a time step of one day
t=arange(0,70/365,1/365)

# This section is setting up different initial values and parameters for different simulations representing microbial community traits
# Here we set up an empty python dictionary to hold the different sets of parameters and initial values
initvals={}
paramsets={}
envir_params={}

# Microbial community, no burn
initvals['no burn sandy soil']=copy.deepcopy(SOM_init)      # Makes a copy of the default initial values. Need to use deepcopy so changing the value here doesn't change it for every simulation
paramsets['no burn sandy soil']=copy.deepcopy(params)
envir_params['no burn sandy soil']=copy.deepcopy(envir_vals)

# Simulation: Incubation of sandy soil following high severity burn
initvals['high sev burn sandy soil']=copy.deepcopy(SOM_init)       # Makes a copy of the default initial values. Need to use deepcopy so changing the value here doesn't change it for every simulation
initvals['high sev burn sandy soil']['MBC_1']=array(1) # Start with low initial microbial biomass. Assumes this simulation starts right after the fire
initvals['high sev burn sandy soil']['MBC_2']=array(0.3)
initvals['high sev burn sandy soil']['MBC_3']=array(0.0)
initvals['high sev burn sandy soil']['MBC_4']=array(0.0)
initvals['high sev burn sandy soil']['uFastC']=array(0.75) 
initvals['high sev burn sandy soil']['uNecroC']=array(0.5) 
initvals['high sev burn sandy soil']['uSlowC']=array(80) # Changed value from 80 to 100 based on calibration results
initvals['high sev burn sandy soil']['uPyC']=array(15)
paramsets['high sev burn sandy soil']=copy.deepcopy(params)
paramsets['high sev burn sandy soil']['vmaxref']['Fast']=75.0     # 
paramsets['high sev burn sandy soil']['vmaxref']['Slow']=0.15
paramsets['high sev burn sandy soil']['vmaxref']['Necro']=75.0
paramsets['high sev burn sandy soil']['vmaxref']['Py']=0.05
paramsets['high sev burn sandy soil']['eup']['Fast'] = 0.5        # 
paramsets['high sev burn sandy soil']['eup']['Slow'] = 0.3
paramsets['high sev burn sandy soil']['eup']['Necro'] = 0.5
paramsets['high sev burn sandy soil']['eup']['Py'] = 0.1
envir_params['high sev burn sandy soil']=copy.deepcopy(envir_vals)
envir_params['high sev burn sandy soil']['thetamin']=array(0.5)
envir_params['high sev burn sandy soil']['thetamax']=array(0.7)
envir_params['high sev burn sandy soil']['porosity']=array(0.4)


# Set up a data structure to hold the results of the different simulations
results={}
# Goes through each functional type and runs a simulation using the appropriate set of parameters and initial values
# Simulations are assuming a constant temperature of 20 C and constant moisture of 60% of saturation
# Inputs are empty because this is running as an incubation without any constant inputs of C

for functype in initvals:
    results[functype] = CORPSE_solvers.run_models_ODE(Tmin=18.0,Tmax=24.0,thetamin=envir_params[functype]['thetamin'],
                                                      thetamax=envir_params[functype]['thetamax'],
                                            times=t,inputs={},clay=2.5,initvals=initvals[functype],params=paramsets[functype])


# Tally total number of microbial pools being used in simulation
from numpy import where
num_micro_pools=0
num_micro_pools=where(SOM_init['MBC_1']>0, num_micro_pools+1, num_micro_pools+0)   
num_micro_pools=where(SOM_init['MBC_2']>0, num_micro_pools+1, num_micro_pools+0)   
num_micro_pools=where(SOM_init['MBC_3']>0, num_micro_pools+1, num_micro_pools+0)   
num_micro_pools=where(SOM_init['MBC_4']>0, num_micro_pools+1, num_micro_pools+0)   


# This section plots the results
# Each set of results should have the same set of pools as the initial values structure from the beginning of the simulation
from matplotlib import pyplot

# Plot CO2 fluxes
fig,ax=pyplot.subplots(nrows=1,ncols=1,clear=True,num='CORPSE results')

for sim in results:
    totalC=CORPSE_array.sumCtypes(results[sim][0], 'u')+CORPSE_array.sumCtypes(results[sim][0], 'p')
    ax.plot(t*365,results[sim][0]['CO2'].diff()/totalC[0]*100,label=sim)
    # ax[1].plot(t*365,results[sim][0]['uFastC'],label='Simple')
    # ax[1].plot(t*365,results[sim][0]['uSlowC'],label='Complex')
    # ax[1].plot(t*365,results[sim][0]['uNecroC'],label='Necromass')
ax.set_xlabel('Time (days)')
ax.set_ylabel('CO$_2$ flux rate (% initial C/day)')
ax.legend(fontsize='small')
ax.set_title('CO$_2$ fluxes')

pyplot.show()
    


# Plot microbial pool sizes
nrows=int(num_micro_pools)

fig,ax=pyplot.subplots(nrows=nrows,ncols=1,clear=True,num='CORPSE results')
for sim in results:
    if nrows == 1: 
        ax.plot(t*365,results[sim][0]['MBC_1']/totalC[0]*100)
        ax.set_ylabel('MBC 1')
        ax.set_xlabel('Time (days)')
        ax.legend(fontsize='small')
        ax.set_title('Microbial biomass C pool size (% of initial C)')
    elif nrows == 2: 
        ax[0].plot(t*365,results[sim][0]['MBC_1']/totalC[0]*100)
        ax[0].set_ylabel('MBC 1')
        ax[1].plot(t*365,results[sim][0]['MBC_2']/totalC[0]*100)  
        ax[1].set_ylabel('MBC 2')
        ax[1].set_xlabel('Time (days)')
    elif nrows ==3: 
        ax[0].plot(t*365,results[sim][0]['MBC_1']/totalC[0]*100)
        ax[0].set_ylabel('MBC 1')
        ax[1].plot(t*365,results[sim][0]['MBC_2']/totalC[0]*100)  
        ax[1].set_ylabel('MBC 2')
        ax[2].plot(t*365,results[sim][0]['MBC_3']/totalC[0]*100)
        ax[2].set_ylabel('MBC 3')
        ax[2].set_xlabel('Time (days)')
    elif nrows == 4:
        ax[0].plot(t*365,results[sim][0]['MBC_1']/totalC[0]*100)
        ax[0].set_ylabel('MBC 1')
        ax[1].plot(t*365,results[sim][0]['MBC_2']/totalC[0]*100)  
        ax[1].set_ylabel('MBC 2')
        ax[2].plot(t*365,results[sim][0]['MBC_3']/totalC[0]*100)
        ax[2].set_ylabel('MBC 3')
        ax[3].plot(t*365,results[sim][0]['MBC_4']/totalC[0]*100)
        ax[3].set_ylabel('MBC 4')
        ax[3].set_xlabel('Time (days)')

ax[0].legend(fontsize='small')
ax[0].set_title('Microbial biomass C pool size (% of initial C)')
# ax[1].set_title('Microbial biomass')

pyplot.show()