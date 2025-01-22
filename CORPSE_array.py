# This file holds the definition and key functions of the actual CORPSE model

# List of parameters the model needs
expected_params={'vmaxref': 'Relative maximum enzymatic decomp rates for each microbe pool (length 3)',
        	'Ea':	'Activation energy (length 3)',
        	'kC':	'Michaelis-Menton parameter (length 3)',
        	'gas_diffusion_exp': 'Determines suppression of decomp at high soil moisture',
            'substrate_diffusion_exp':'Determines suppression of decomp at low soil moisture',
        	'minMicrobeC':	   'Minimum microbe biomass (fraction of total C)',
        	'Tmic': 'microbe lifetime at 20C (years)',
        	'et':  'Fraction of turnover from first microbial pool not converted to CO2',
        	'eup': 'Carbon uptake efficiency (length 3)',
        	'tProtected':	'Protected C turnover time (years)',
        	'protection_rate':'Protected carbon formation rate (year-1) (length 3)',
            'new_resp_units':True,
            }


# Names of the C types. Can edit this to change the number and name of pools in edited model simulations
chem_types = ['Fast','Slow','Necro', 'Py']

# set number of microbial pools (up to 4) to use - would be nice to integrate this into Microbial_sims.py
microbial_pools = ['MBC_1','MBC_2','MBC_3','MBC_4'] 

# Makes a list of the pools that should actually be in the model, including both protected and unprotected states
expected_pools = ['u'+t+'C' for t in chem_types]+\
                 ['p'+t+'C' for t in chem_types]+\
                 [m for m in microbial_pools]+\
                 ['CO2','originalC']

# Function for updating protected C formation rate for different soil textures, based on Mayes et al (2012) paper
# All soils: slope=0.4833,intercept=2.3282
# Alfisols: slope=0.5945, intercept=2.2788
def prot_clay(claypercent,slope=0.4833,intercept=2.3282,BD=1.15,porosity=0.4):
    ''' Calculate protection rate as a function of clay content, based on sorption isotherms from Mayes et al (2012) Table 3
    Calculates Qmax in mgC/kg soil from Mayes et al 2012, converted to g/m3 using bulk density
    Typically used as relative value for calculating protection_rate parameter.
    claypercent: Soil % clay (out of 100)
    slope: Either for all soils or a soil order, from Mayes paper
    intercept: Either for all soils or a soil order, from Mayes paper
    BD: Soil bulk density in g/cm3
    '''
    from numpy import log10
    prot=1.0*(10**(slope*log10(claypercent)+intercept)*BD*1e-6)
    return prot

# Check if the parameters sent to the model included the correct set of parameters and raise an error if not
def check_params(params):
    '''params: dictionary containing parameter values. Should contain these fields (showing reasonable default values):
             vmaxref=[2500,600,2000]; Relative maximum enzymatic decomp rates
             Ea=[37e3,54e3,50e3];     Activation energy
             kC=[0.01,0.01,0.01];     Michaelis-Menton parameter
             gas_diffusion_exp=2.5;   Determines suppression of decomp at high soil moisture
             minMicrobeC=1e-3;       Minimum microbe biomass (fraction of total C)
             Tmic=0.15;        microbe turnover rate
             et=0.5;           Fraction of turnover not converted to CO2
             eup=[0.6,0.05,0.6];  Carbon uptake efficiency - slow growing microbes
             tProtected=75.0;     Protected C turnover time (years)
             protection_rate=[1.0,0.0,1.0];  Protected carbon formation rate (year-1)'''

    #from numpy import iterable,array
    unused_params=expected_params.copy()
    for k in params.keys():
        if k not in expected_params:
            raise ValueError('Parameter set contains unexpected parameter %s'%k)
        unused_params.pop(k)
        #if iterable(params[k]):
         #   params[k]=array(params[k])
    if len(unused_params)>0:
        for k in unused_params.keys():
            print ('Missing parameter: %s [%s]'%(k,unused_params[k]))
        raise ValueError('Missing parameters: %s'%unused_params.keys())

# The main model function. Given the current state of the model along with temperature, moisture, and parameters, it calculates the rate of change of all pools
from numpy import zeros,size,where,atleast_1d
def CORPSE_deriv(SOM,T,theta,params,claymod=1.0):
    '''Calculate rates of change for all CORPSE pools
       T: Temperature (K)
       theta: Soil water content (fraction of saturation)

       Returns same data structure as SOM'''

    # convert input into array
    theta=atleast_1d(theta)
    T=atleast_1d(T)

    # Constraint theta to 0 < theta < 1
    theta[theta<0]=0.0
    theta[theta>1]=1.0

    # Calculate maximum potential C decomposition rate of each C pool by each microbial group
    decomp=decompRate(SOM,T,theta,params)

    # Create empty variables for dead MBC and microbial turnover rate
    deadmic_C_production=0
    microbeTurnover={}
    maintenanceResp={}
    total_maintenanceResp=0
    CO2prod=0
    microbeGrowth={}
    
    # Calculate microbial turnover and dead MBC production for each microbial pool. If MBC pool > 0, calculate microbial turnover. If not, set turnover equal to 0.
    microbeTurnover=dict([(m, where((SOM[m]>0), (SOM[m]-params['minMicrobeC'][m]*(sumCtypes(SOM,'u')))/params['Tmic'][m],0)) for m in microbial_pools])
    
    for m in microbial_pools:    
        # Ensure that microbial turnover is > 0
        if isinstance(microbeTurnover[m],float):
            microbeTurnover[m]=max(0.0,microbeTurnover[m])
        else:
            microbeTurnover[m][microbeTurnover[m]<0.0]=0.0
            
        deadmic_C_production=deadmic_C_production+microbeTurnover[m]*params['et'][m]
        maintenanceResp[m]=microbeTurnover[m]*(1.0-params['et'][m])
        total_maintenanceResp=total_maintenanceResp+maintenanceResp[m]
        CO2prod=CO2prod+maintenanceResp[m]
    
        microbeGrowth[m]=0
        for t in chem_types: 
            CO2prod=CO2prod + decomp[m][t]*(1.0-params['eup'][m][t])
            microbeGrowth[m]=microbeGrowth[m]+decomp[m][t]*params['eup'][m][t]




    # Update protected carbon
    protectedCturnover = dict([(t,SOM['p'+t+'C']/params['tProtected']) for t in chem_types])
    protectedCprod =     dict([(t,SOM['u'+t+'C']*params['protection_rate'][t]*claymod) for t in chem_types])
        # Some slow, protected C is being formed, but none of the pools are undergoing turnover > 0. Where is the protected C coming from??
        
        
    #### Assign new values to SOM pools (based on calculated rates of change)
    derivs=SOM.copy() # Create a copy of the original SOM dictionary
    for k in derivs.keys():
        derivs[k]=0.0 # Set all the C pools to 0
        
    # Fill in dictionary with new MBC pool size based on microbial growth and turnover 
    for m in microbial_pools:
        derivs[m]=microbeGrowth[m]-microbeTurnover[m]

    derivs['CO2']=CO2prod

    total_decomp={}   # create empty dictionary 
    for t in chem_types:
        total_decomp['u'+t+'C']=0 # Assign an initial value of 0 to for decomposition rate of C pool
        total_decomp['p'+t+'C']=0
        for m in microbial_pools:
            # Calculate cumulative decomposition of each pool by each microbial group
            total_decomp['u'+t+'C']=total_decomp['u'+t+'C']+decomp[m][t]
        # Fill in derivs dictionary with new C pool (Fast, Slow, Necro, Py) size based on decomposition rate, protected C formation, and protected C turnover
        derivs['u'+t+'C']=-total_decomp['u'+t+'C']+protectedCturnover[t]-protectedCprod[t]
        derivs['p'+t+'C']=-total_decomp['p'+t+'C']+protectedCprod[t]-protectedCturnover[t]

    # Add new dead MBC to the necromass pool 
    derivs['uNecroC']=derivs['uNecroC']+deadmic_C_production
    return derivs


# Decomposition rate
def decompRate(SOM,T,theta,params):

    # This only really needs to be calculated once
    # This corrects the units of vmaxref for the moisture function, so the 
    # units can be in actual 1/time. 
    # This essentially sets the maximum value of the moisture function to 1
    if params['new_resp_units']:
        theta_resp_max=params['substrate_diffusion_exp']/(params['gas_diffusion_exp']*(1.0+params['substrate_diffusion_exp']/params['gas_diffusion_exp']))
        aerobic_max=theta_resp_max**params['substrate_diffusion_exp']*(1.0-theta_resp_max)**params['gas_diffusion_exp']

    else:
        aerobic_max=1.0
 
    vmax=dict([(m, Vmax(m,T,params)) for m in microbial_pools])
        
    # Skip the decomposition calculation if there is no carbon or no microbe biomass (to avoid dividing by zero)
    dodecomp=dict([(m,(sumCtypes(SOM,'u')!=0.0)&(theta!=0.0)&(SOM[m]!=0.0)) for m in microbial_pools])

    # Create empty dictionary to store decomposition rates
    decompRate={}
    
    for m in microbial_pools:
        decompRate[m]=dict([(t, where(dodecomp[m],vmax[m][t]*theta**params['substrate_diffusion_exp']*(SOM['u'+t+'C'])*SOM[m]/(sumCtypes(SOM,'u')*params['kC'][m][t]+(SOM['MBC_1']+SOM['MBC_2']+SOM['MBC_3']+SOM['MBC_4']))*(1.0-theta)**params['gas_diffusion_exp']/aerobic_max,0.0)) for t in chem_types])

    return decompRate # output is the decomposition rates of each C pool by given MBC group


def Vmax(Micro_pool, T,params):
    '''Vmax function, normalized to Tref=293.15
    T is in K'''

    Tref=293.15;
    Rugas=8.314472;

    from numpy import exp

    # Calculate temperature adjusted vmax for each chem type
    Vmax=dict([(t,params['vmaxref'][Micro_pool][t]*exp(-params['Ea'][t]*(1.0/(Rugas*T)-1.0/(Rugas*Tref)))) for t in chem_types]);

    return Vmax

# Add together the C types. prefix is for specifying protected or unprotected (p or u)
#  Doesn't include living MBC or CO2
def sumCtypes(SOM,prefix):
    out=SOM[prefix+chem_types[0]+'C']
    if len(chem_types)>1:
        for t in chem_types[1:]:
            out=out+SOM[prefix+t+'C']

    return out
