
function eulers_method(my_func,x_t,params,h)
    """
    input: 
        my_func<function>: a function that specifies the time derivative of a dynamical system
        x_t: the state of the system
        h: the Δt to use
    """
    dX = my_func(x_t,params)
    x_t_new = x_t+dX.*h#the next value of the sytem is equal to the previous value times the derivatve at that point times a time step 
    return x_t_new
end


function heuns_method(my_func,x_t,params,h)
    """
    input: 
        my_func<function>: a function that specifies the time derivative of a dynamical system
        x_t: the state of the system
        h: the Δt to use
    """
    x_t_prime = eulers_method(my_func,x_t,params,h)#the next value of the sytem is equal to the previous value times the derivatve at that point times a time step 
    dX = my_func(x_t,params)
    dX_prime = my_func(x_t_prime,params)#perform eulers method using the output of eulers method as an input
    x_t_new = x_t + h/2*(dX+dX_prime)#combine for heuns method
    return x_t_new
end

function ODESolver(equation,initial_conditions,params,max_time,num_steps=1000,integrator_method=heuns_method) 
    """
    This is the function that is called, given an ODE in the form of a function and a set of initial condittions will return the time series of the solution as a vector where each element is a vector storing the values of each variable at each time step
    inputs: 
        equation<function>: this function should take in two parameters X, a vector of the variables in the differential equations and p, a vector of the paramters, returns dX, a vector of the differenital increments to be added to each variable
        inital_conditions<vector>: a vector of the inital coniditions for each vairable, our X₀. 
        params<vector>: a vector of the paramteres passed into <equation>
        max_time<float>: the maximium time at which to stop the integration. NOTE: this is not the number of runs of the integrator but the maximum value of the time variable. 
        num_steps<float>: the number of integrator steps to take over this time. The step size is then <max_time>/<num_steps>
        integrator_method<function>: either eulers_method or heuns_method, this is the method used to numerically solve the ODE
    returns: 
        time_range<vector>: the values of the time variable for each step
        x_values<vector>: a vector where each element of a vector of the values of the variables in our system of differntial equations at each time step. 
    """
    ΔT= max_time/num_steps#calculate the time interval of each step
    time_range = collect(0:ΔT:max_time)
    x_values = []

    x = initial_conditions 
    for i in time_range
        #global x
        push!(x_values,x)
        x = integrator_method(equation,x,params,ΔT)
    end

    return time_range,x_values    
end

