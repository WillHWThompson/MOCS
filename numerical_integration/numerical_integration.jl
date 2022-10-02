
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
    ΔT= max_time/num_steps#calculate the time interval of each step
    time_range = collect(0:num_steps)
    x_values = []

    x = initial_conditions 
    for i in time_range
        #global x
        push!(x_values,x)
        x = integrator_method(equation,x,params,ΔT)
    end

    return time_range,x_values    
end

