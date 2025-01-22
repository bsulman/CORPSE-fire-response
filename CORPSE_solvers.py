import CORPSE_array as CORPSE_deriv
import numpy
fields = CORPSE_deriv.expected_pools
fields.remove('originalC')

# This is a function that translates the CORPSE model pools to/from the format that the equation solver expects
# The solver will call it multiple times and passes it a list of parameters that needs to be converted to a named "dictionary" that CORPSE expects
def fsolve_wrapper(SOM_list,T,theta,inputs,clay,params):
    from numpy import atleast_1d,concatenate

    # Make an empty dictionary and fill it with the right values
    # Values are in the correct order because we use the same "fields" list of pools
    SOM_dict={}
    for n in range(len(fields)):
        SOM_dict[fields[n]]=atleast_1d(SOM_list[n])
        #SOM_dict[fields[n]]=atleast_1d(SOM_list[fields[n]])


    # Call the CORPSE model function that returns the derivative (with time) of each pool
    deriv=CORPSE_deriv.CORPSE_deriv(SOM_dict,atleast_1d([T]),atleast_1d([theta]),params,claymod=CORPSE_deriv.prot_clay(clay)/CORPSE_deriv.prot_clay(20))

    # Since we have carbon inputs, these also need to be added to those rates of change with time
    for pool in inputs.keys():
        deriv[pool]+=inputs[pool]


    # Put the CORPSE pools back into a list that the equation solver can deal with
    vals=list(concatenate([deriv[f] for f in fields]))

    if len(SOM_list)==len(fields)*2:
        vals=vals+vals_Ohorizon

    return vals

# The ordinary differential equation (ODE) integrating function also wants to send the current time to the function it's integrating
# Our model doesn't have an explicit dependence on time, but we need a separate function that can deal with the extra argument.
# We just ignore the time argument and pass the rest to the same function we used for the numerical solver
def ode_wrapper(SOM_list,time,Tmax,Tmin,thetamax,thetamin,*args,**kwargs):
    from numpy import cos,pi
    T=(cos(time*2*pi)+1)*(Tmax-Tmin)/2+Tmin
    theta=(cos(time*2*pi)+1)*(thetamax-thetamin)/2+thetamin
    return fsolve_wrapper(SOM_list,T,theta,*args,**kwargs)

# Uses an alternate method: Iterating through time steps but loading all points into a vector for more efficient calculation
# May run faster for large number of points, but potentially less accurate depending on time step
def vector_iterate(SOM_init,params,T,theta,inputs,clay,times):
    from numpy import zeros,atleast_1d
    # totaltime and dt in units of years
    nsteps=len(times)
    if len(T.shape)>1:
        npoints=T.shape[1]
    else:
        npoints=len(T)
    nrecords=nsteps

    # Set up pools
    SOM={}
    SOM_out={}
    for field in SOM_init.keys():
        SOM_out[field]=zeros((npoints,nrecords))
        if len(atleast_1d(SOM_init['uFastC'])) == 1: 
            SOM[field]=zeros(npoints)+SOM_init[field]
        else:
            SOM[field]=SOM_init[field].values

    # Iterate through simulations
    for step in range(nsteps):
        if step==nsteps-1:
            dt=times[step]-times[step-1]
        else:
            dt=times[step+1]-times[step]
        if len(T.shape)>1:
            T_step=T[step,:]
        else:
            T_step=T
        if len(theta.shape)>1:
            theta_step=theta[step,:]
        else:
            theta_step=theta
        # In this case, T, theta, clay, and all the pools in SOM are vectors containing one value per geographical location
        deriv=CORPSE_deriv.CORPSE_deriv(SOM,T_step,theta_step,params,claymod=CORPSE_deriv.prot_clay(clay.values)/CORPSE_deriv.prot_clay(2.5))

        # Since we have carbon inputs, these also need to be added to those rates of change with time
        for pool in inputs.keys():
            deriv[pool]+=inputs[pool]

        for field in deriv.keys():
            SOM[field]=SOM[field]+deriv[field]*dt

        if (step*dt)%10==0:
            print('Time = %d'%(step*dt))
        for field in SOM.keys():
            SOM_out[field][:,step]=SOM[field]

    return SOM_out

# This function runs an actual simulation using the ODE solver
# Tmin and Tmax allow a sinusoidal temperature variation. Similar for thetamin and thetamax. Set min and max equal for constant state
# times is an array of all the time steps for the simulation
# clay is clay content (%) used to adjust protected C formation rates
# initvals is a dictionary of the initial values for all pools
# params is a dictionary of all the parameter values
# inputs is a dictionary of the C input rates of all pools. Assumes zero rate for pools not in the inputs data structure (so it can be empty for no inputs)
def run_models_ODE(Tmin,Tmax,thetamin,thetamax,times,inputs,params,clay,initvals):
    import time,pandas
    from numpy import atleast_1d
    t0=time.time()
    from scipy.integrate import odeint
    # Again, start with an empty list to hold the output
    SOM_out_ODE=[]
    print('ODE integrator')


    def get_initvals(initvals,point):
        # Set initial values into a list to give the solver
        # I'm using a convenient piece of Python syntax for making lists
        if isinstance(initvals,dict):
            ivals=[initvals[f] for f in fields]
        elif isinstance(initvals, pandas.DataFrame):
            ivals=[initvals.iloc[point][f] for f in fields]
        else:
            if is_numlike(initvals[0]):
                ivals=initvals
            else:
                init_sim=initvals[point]
                ivals=[init_sim.iloc[-1][f] for f in fields]

        return ivals

    for point in range(len(atleast_1d(clay))):
        print('Point %d of %d'%(point,len(atleast_1d(clay))))

        ivals=get_initvals(initvals,point)

        # Runs the ODE integrator
        result=odeint(ode_wrapper,ivals,times,
            args=(atleast_1d(Tmax)[point]+273.15,atleast_1d(Tmin)[point]+273.15,atleast_1d(thetamax)[point],atleast_1d(thetamin)[point],inputs,atleast_1d(clay)[point],params))
        # Store the output in a pandas DataFrame (similar to R's dataframes)
        result_df=pandas.DataFrame(result[:,:len(fields)],columns=fields,index=times)

        # Add it to the list of the output from each point
        SOM_out_ODE.append(result_df)


    print('Time elapsed: %1.1f s'%(time.time()-t0))

    return SOM_out_ODE

# # Run a simulation using the explicit iterator instead of the ODE solver. Can edit this function to allow more complex temperature and moisture patterns, among other things
# def run_models_iterator(Tmin,Tmax,thetamin,thetamax,times,inputs,params,clay,initvals):
#     # Iterate explicitly
#     import time
#     from pandas import DataFrame
#     from numpy import arange
#     t0=time.time()
    
#     from numpy import cos,pi
#     T=(cos(times[:,None]*2*pi)+1)*(Tmax.values-Tmin.values)/2+Tmin.values
#     theta=(cos(times[:,None]*2*pi)+1)*(thetamax.values-thetamin.values)/2+thetamin.values
    
#     dt=times[1]-times[0]
#     result_iterator=vector_iterate(initvals,params,T+273.15,theta,inputs,clay,times)
#     SOM_out_iterator=[]
#     for point in range(len(Tmin)):
#         df=DataFrame(index=times,columns=fields)
#         for field in result_iterator.keys():
#             df[field]=result_iterator[field][point,:]
#         # df['livingMicrobeN']=df['livingMicrobeC']/params['CN_microbe']
#         SOM_out_iterator.append(df)

#     print('Time elapsed: %1.1f s'%(time.time()-t0))
#     return SOM_out_iterator

# Functions for adding together all the C pools. They work on either dictionary or dataframe data types because both have the same names for the pools
def totalCarbon(SOM, microbial_pools):
    totalMBC=0
    for m in microbial_pools:
        totalMBC=totalMBC+SOM[m]
    return SOM['uFastC']+SOM['uSlowC']+SOM['uNecroC']+SOM['uPyC']+SOM['pFastC']+SOM['pSlowC']+SOM['pNecroC']+SOM['pPyC']+totalMBC
    
