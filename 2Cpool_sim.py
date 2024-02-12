# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 14:16:52 2023

@author: danab
"""
from numpy import arange
from numpy import exp

# Initial pools for soil organic matter simulation: Fast and slow C
# "Fast" has higher maximum decomposition rate
params = {
    "initC":1.0,
    "fracC":{'fast':0.05,'slow':0.95},
    "DecompRate":{'fast':0.1,"slow":0.00001}
    } 

# Array of model time steps (in units of years). 
# Start = 0, End = 120 days, timestep = 1 day
t=arange(0,60/365,1.0/365)

# Names of the C types. 
chem_types = ['fast','slow']



# Modelling total C remaining
def TotalC_remaining(t,params):

    # Objective: Calculate C respired (1 - Mt) at given t
    # Mt = fractional C remaining = M1*e^(-k1*t) + M2*e^(-k2*t)
    #C_remaining = [(params['fracC'][i]*exp(-params['DecompRate'][i]*t)) for i in chem_types];
    C_remaining = [(params['fracC']['fast']*exp(-params['DecompRate']['fast']*t)) + (params['fracC']['slow']*exp(-params['DecompRate']['slow']*t))]
    C_resp = [1-C_remaining[0]]
    return C_resp
    # Doesn't work yet.

TotalC_remaining=TotalC_remaining(t, params)



# Plot results
from matplotlib import pyplot

fig,ax=pyplot.plot(t*365,TotalC_remaining[0])

# Why is this linear? Because I've created a linear model in which the fractional
# size of the fast and slow C pools doesn't depend on the previous timestep.

# I should take what I learned hear and attempt to modify the CORPSE model to just
# take fast and slow C pools as inputs. Maybe also temperature.


