using JLD2
using JLD2: Group
using Oceananigans
using Oceananigans.Units: Time

function set_data_on_disk!(original_data_filename; direction = "backward", T_start = nothing, T_end = nothing)

    # First check that a valid direction has been given
    if (direction !== "backward") && (direction !== "forward")
        error("Invalid direction: $direction")
    end
    
    # Open the file
    jldopen(original_data_filename,"r+") do file

        # Checking if this is the first time we've edited the file
        if !haskey(file, "direction")
            # It should currently be forward, because direction hasn't been assigned yet
            file["direction"] = "forward"

            # This is also the first time we've edited the file, so we'll create some data to store original simulation times
            iterations = parse.(Int, keys(file["timeseries/t"]))
            times = [file["timeseries/t/$iter"] for iter in iterations]
            g = Group(file, "timeseries/t_simulation")
            for (i,iter) in enumerate(iterations)
                g["$iter"] = times[i]
            end
        end    

        # Now load in current direction and simulation times
        current_direction = file["direction"]
        println("Current direction is $current_direction")
        iterations = parse.(Int, keys(file["timeseries/t"]))
        t_simulation = [file["timeseries/t_simulation/$iter"] for iter in iterations]
        Nt = length(t_simulation)

        # If already the right direction, we just need to make sure the times correspond to the filter period we want
        if current_direction == direction
            println("No need to reverse order of data")
            # We might need to update times
            if (current_direction == "forward") 

                if isnothing(T_start) # If we're not given a starting time, set it to simulation start time 
                    T_start = t_simulation[1]
                    
                end

                if isnothing(T_end) 
                    T_end = t_simulation[Nt] # If we're not given a end time, set it to simulation end time 
                    
                end

                # Make a new time coordinate that is zero at T_start
                Base.delete!(file, "timeseries/t")
                g = Group(file, "timeseries/t")
                for (i,iter) in enumerate((iterations))
                    g["$iter"] = t_simulation[i] - T_start
                end
                
            else # current direction is backward

                if isnothing(T_start) # If we're not given a starting time, set it to simulation start time 
                    T_start = t_simulation[Nt]
                    
                end

                if isnothing(T_end) 
                    T_end = t_simulation[1] # If we're not given a end time, set it to simulation end time 
                    
                end

                # Make a new time coordinate that is zero at T_end
                Base.delete!(file, "timeseries/t")
                g = Group(file, "timeseries/t")
                for (i,iter) in enumerate(iterations)
                    g["$iter"] = T_end - t_simulation[i]
                end 
            end

        else # current direction is wrong, so we need to reverse everything and also update times
            Base.delete!(file, "direction")
            file["direction"] = direction
            println("Reversing order of data")
            for var in keys(file["timeseries"])
                if (var == "t_simulation") || (var == "t") # Then no serialized entry
                    backward_iterations = reverse(keys(file["timeseries/$var"])) 
                else
                    backward_iterations = reverse(keys(file["timeseries/$var"])[2:end]) # don't include serialized
                end
                # for velocities, we reverse all entries and negate them
                if (var == "u") || (var == "v") || (var == "w")
                    for iter in backward_iterations 
                        data = file["timeseries/$var/$iter"] # Load in data
                        Base.delete!(file, "timeseries/$var/$iter") # Delete the entry
                        file["timeseries/$var/$iter"] = -data # Write it again, but negative                  
                    end
                    println("Reversed order of and switched sign of $var")

                # For time, we need to reverse the data and also make sure its correctly aligned with start and end points
                elseif var == "t"
                    if current_direction == "forward"
                        if isnothing(T_start) # If we're not given a starting time, set it to simulation start time 
                            T_start = t_simulation[1]
                        end
        
                        if isnothing(T_end) 
                            T_end = t_simulation[Nt] # If we're not given a end time, set it to simulation end time                             
                        end

                        for (i, iter) in enumerate(backward_iterations)
                            Base.delete!(file, "timeseries/$var/$iter") # Delete the entry
                            file["timeseries/$var/$iter"] = T_end - t_simulation[Nt+1-i]
                        end

                    else # current_direction is backward

                        if isnothing(T_start) # If we're not given a starting time, set it to simulation start time 
                            T_start = t_simulation[Nt]
                            println("New T_start is $T_start")
                        end
        
                        if isnothing(T_end) 
                            T_end = t_simulation[1] # If we're not given a end time, set it to simulation end time 
                            println("New T_end is $T_end")
                        end

                        for (i, iter) in enumerate(backward_iterations)
                            Base.delete!(file, "timeseries/$var/$iter") # Delete the entry
                            file["timeseries/$var/$iter"] = t_simulation[Nt+1-i] - T_start
                        end
                    end
                    println("Reversed order of and shifted $var")
                # if var is t_simulation, the data doesn't change
                elseif var == "t_simulation"
                    for iter in backward_iterations 
                        data = file["timeseries/$var/$iter"] # Load in data
                        Base.delete!(file, "timeseries/$var/$iter") # Delete the entry
                        file["timeseries/$var/$iter"] = data # Write it again 
                    end
                    println("Reversed order of $var")
                # if var is a scalar, we need to make sure that it's on the correct (center) grid   
                else
                    location = file["timeseries/$var/serialized/location"]
                    if location !== (Center, Center, Center)
                        @warn "$var is on grid $location, but should be output on (Center, Center, Center). Continuing anyway."
                    end 
                    for iter in backward_iterations 
                        data = file["timeseries/$var/$iter"] # Load in data
                        Base.delete!(file, "timeseries/$var/$iter") # Delete the entry
                        file["timeseries/$var/$iter"] = data # Write it again 
                    end
                    println("Reversed order of $var")
                end
            end
            println("New direction is $direction")
        end
        T_filter = T_end - T_start # Update total filter time
        return T_filter 
    end
      
end

function load_data(original_data_filename, original_var_names, velocity_names; architecture=CPU(), backend= InMemory())
    velocity_timeseries = ()
    for vel_name in velocity_names
        vel_ts = FieldTimeSeries(original_data_filename, vel_name; architecture=architecture, backend=backend)
        velocity_timeseries = (velocity_timeseries...,vel_ts)
    end
   
    var_timeseries = ()
    for var_name in original_var_names
        var_ts = FieldTimeSeries(original_data_filename, var_name; architecture=architecture, backend=backend)
        var_timeseries = (var_timeseries...,var_ts)
    end
    grid = velocity_timeseries[1].grid

    return velocity_timeseries, var_timeseries, grid
end

function set_BW_filter_params(;N=1,freq_c=1) 
    N_coeffs = 2^(N-1)
    filter_params = NamedTuple()
    for i in 1:N_coeffs
        
        a = (freq_c/2^N)*sin(pi/(2^(N+1))*(2*i-1))
        b = (freq_c/2^N)*cos(pi/(2^(N+1))*(2*i-1))
        c = freq_c*sin(pi/(2^(N+1))*(2*i-1))
        d = freq_c*cos(pi/(2^(N+1))*(2*i-1))

        temp_params = NamedTuple{(Symbol("a$i"), Symbol("b$i"),Symbol("c$i"),Symbol("d$i"))}([a,b,c,d])
        filter_params = merge(filter_params,temp_params)
    end
    
    return filter_params
end

function create_original_vars(original_var_names, grid)
    # Creates auxiliary fields to store the saved variables
    vars = Dict()
    for var_name in original_var_names
        vars[Symbol(var_name)] = CenterField(grid)
    end
    return NamedTuple(vars)
end

function create_filtered_vars(original_var_names, velocity_names, filter_params; map_to_mean=true)
    N_coeffs = Int(length(filter_params)/4)
    gC = ()  # Start with an empty tuple
    gS = ()  # Start with an empty tuple
    for var_name in original_var_names
        for i in 1:N_coeffs
            new_gC = Symbol(var_name,"C", i) 
            new_gS = Symbol(var_name,"S", i) 
            gC = (gC..., new_gC)  
            gS = (gS..., new_gS)  
        end
    end

    # May also need xi maps. We need one for every velocity dimension, so lets use the velocity names to name them
    if map_to_mean
        for velocity_name in velocity_names
            for i in 1:N_coeffs
                new_gC = Symbol("xi_",velocity_name,"_C", i) 
                new_gS = Symbol("xi_",velocity_name,"_S", i) 
                gC = (gC..., new_gC)  
                gS = (gS..., new_gS)  
            end
        end
    end

    return (gC...,gS...)
end

function make_gC_forcing_func(i)
    c_index = (i-1)*4 + 3
    d_index = (i-1)*4 + 4
    # Return a new function. In 3D, this would be (x,y,z,t,gC,gS,p) -> -p[c_index]*gC - p[d_index]*gS
    return (args...) -> -args[end][c_index]*args[end-2] - args[end][d_index]*args[end-1]
end

function make_gS_forcing_func(i)
    c_index = (i-1)*4 + 3
    d_index = (i-1)*4 + 4
    # Return a new function. In 3D, this would be (x,y,z,t,gC,gS,p) -> -p[c_index]*gS + p[d_index]*gC
    return (args...) -> -args[end][c_index]*args[end-1] + args[end][d_index]*args[end-2]
end

function make_xiC_forcing_func(i)
    c_index = (i-1)*4 + 3
    d_index = (i-1)*4 + 4
    # Return a new function. In 3D, this would be (x,y,z,t,u,p) -> -p[c_index]/(p[c_index]^2 + p[d_index]^2)*u
    return (args...) -> -args[end][c_index]/(args[end][c_index]^2 + args[end][d_index]^2)*args[end-1]
end

function make_xiS_forcing_func(i)
    c_index = (i-1)*4 + 3
    d_index = (i-1)*4 + 4
    # Return a new function. In 3D, this would be (x,y,z,t,u,p) -> -p[d_index]/(p[c_index]^2 + p[d_index]^2)*u
    return (args...) -> -args[end][d_index]/(args[end][c_index]^2 + args[end][d_index]^2)*args[end-1]
end

function create_forcing(filtered_vars, original_var_names, velocity_names, filter_params)

    N_coeffs = Int(length(filter_params)/4)

    # Initialize dictionary
    gCdict = Dict()
    gSdict = Dict()

    # Make forcing function for original data term. Final arg is the field dependence.
    original_var_forcing_func(args...) = args[end]

    for var_name in original_var_names
        original_var = Symbol(var_name)
        for i in 1:N_coeffs
            
            gCkey = Symbol(var_name,"C",i)   # Dynamically create a Symbol for the key
            gSkey = Symbol(var_name,"S",i)   # Dynamically create a Symbol for the key
            
            # Store in dictionary
            gC_forcing_i = Forcing(make_gC_forcing_func(i), parameters = filter_params, field_dependencies = (gCkey,gSkey))
            gC_original_var_forcing = Forcing(original_var_forcing_func,field_dependencies= (;original_var))
            gCdict[gCkey] = (gC_forcing_i, gC_original_var_forcing)

            gS_forcing_i = Forcing(make_gS_forcing_func(i), parameters = filter_params, field_dependencies = (gCkey,gSkey))
            gSdict[gSkey] = gS_forcing_i
            
        end
    end

    # We might need xi forcing too, forced by the model's velocities
    if length(filtered_vars) > N_coeffs*length(original_var_names)*2
        for velocity_name in velocity_names
            vel_key = Symbol(velocity_name)
            for i in 1:N_coeffs
                gCkey = Symbol("xi_",velocity_name,"_C", i) 
                gSkey = Symbol("xi_",velocity_name,"_S", i) 
                
                # Store in dictionary
                gC_forcing_i = Forcing(make_gC_forcing_func(i), parameters =filter_params, field_dependencies = (gCkey,gSkey))
                xiC_forcing_i = Forcing(make_xiC_forcing_func(i),parameters =filter_params, field_dependencies = (;vel_key))
                gCdict[gCkey] = (xiC_forcing_i,gC_forcing_i)

                gS_forcing_i = Forcing(make_gS_forcing_func(i), parameters =filter_params, field_dependencies = (gCkey,gSkey))
                xiS_forcing_i = Forcing(make_xiS_forcing_func(i),parameters =filter_params, field_dependencies = (;vel_key))
                gSdict[gSkey] = (xiS_forcing_i,gS_forcing_i)
            end
        end
    end
    forcing = (; NamedTuple(gCdict)..., NamedTuple(gSdict)...)

    return forcing
end

function create_output_fields(model, original_var_names, velocity_names, filter_params)
    N_coeffs = Int(length(filter_params)/4)
    N_filtered_vars = length(model.tracers)
    half_N_filtered_vars = Int(N_filtered_vars/2)
    outputs_dict = Dict()
    for (i_var, original_var) in enumerate(original_var_names)

        g_total = (filter_params[1]*model.tracers[(i_var-1)*N_coeffs + 1] + filter_params[2]*model.tracers[half_N_filtered_vars + (i_var-1)*N_coeffs + 1])
    
        for i in 2:N_coeffs
            a_index = (i-1)*4 + 1
            b_index = (i-1)*4 + 2
            gC_index = (i_var-1)*N_coeffs + i
            gS_index = half_N_filtered_vars + (i_var-1)*N_coeffs + i
            g_total = g_total + (filter_params[a_index]*model.tracers[gC_index] + filter_params[b_index]*model.tracers[gS_index])
        end
        outputs_dict[original_var*"_filtered"] = g_total
    end

    # May need to output maps too
    if N_filtered_vars > N_coeffs*length(original_var_names)*2
        for (i_vel, vel) in enumerate(velocity_names)
            i_var = i_vel + length(original_var_names)
            g_total = (filter_params[1]*model.tracers[(i_var-1)*N_coeffs + 1] + filter_params[2]*model.tracers[half_N_filtered_vars + (i_var-1)*N_coeffs + 1])
        
            for i in 2:N_coeffs
                a_index = (i-1)*4 + 1
                b_index = (i-1)*4 + 2
                gC_index = (i_var-1)*N_coeffs + i
                gS_index = half_N_filtered_vars + (i_var-1)*N_coeffs + i
                g_total = g_total + (filter_params[a_index]*model.tracers[gC_index] + filter_params[b_index]*model.tracers[gS_index])
            end
            outputs_dict["xi_"*vel] = g_total
        end
    end

    # Let's also add the saved vars for comparison
    for original_var in original_var_names
        outputs_dict[original_var] = getproperty(model.auxiliary_fields,Symbol(original_var))
    end
    return outputs_dict
end

function update_input_data!(sim, input_data)
    velocity_timeseries = input_data.velocities
    original_var_timeseries = input_data.original_vars
    model = sim.model
    t = sim.model.clock.time
    
    kwargs = (; (Symbol(vel_fts.name) => vel_fts[Time(t)] for vel_fts in velocity_timeseries)...)
    set!(model; kwargs...)          

    # We also update the saved variables to be used for forcing - these are auxiliary fields so need to be set separately
    for original_var_fts in original_var_timeseries
        set!(getproperty(model.auxiliary_fields, Symbol(original_var_fts.name)), original_var_fts[Time(t)])
    end
end


# function update_saved_data!(sim, fts_velocities)
#     model = sim.model
#     time = sim.model.clock.time
#     u_fts, v_fts = fts_velocities
#     set!(model, u=u_fts[Time(time)], v= v_fts[Time(time)])
#     return nothing
# end

# function update_vorticity!(sim, fts_vorticity)
#     model = sim.model
#     time = sim.model.clock.time
#     ω_fts = fts_vorticity
#     set!(model.auxiliary_fields.ω, ω_fts[Time(time)])
#     return nothing
# end

