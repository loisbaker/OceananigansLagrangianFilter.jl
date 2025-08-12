
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

#TODO the following shouldn't be done with indexing, do it by keys instead
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
        # halo regions get filled automatically
    end
end


function sum_forward_backward_contributions!(combined_output_filename,forward_output_filename,backward_output_filename,T,velocity_names, original_var_names)
# Combine the forward and backward simulations by summing them into a single file

    # Make a copy of the forward file to fill in with combined data
    cp(forward_output_filename, combined_output_filename,force=true)

    # List the names of the fields that we will combine
    data_field_names = vcat(["xi_" * vel for vel in velocity_names], [var * "_filtered" for var in original_var_names])

    # Then, we open the combined data and add backward scalars and maps at the correct timesteps
    jldopen(combined_output_filename,"r+") do file
        forward_iterations = parse.(Int, keys(file["timeseries/t"]))
        for field in data_field_names
            # Open the backward data as a FieldTimeSeries
            fts_backward = FieldTimeSeries(backward_output_filename, field)

            # Loop over forward times and add the backward data
            for iter in forward_iterations
                forward_time = file["timeseries/t/$iter"]
                forward_data = file["timeseries/$field/$iter"] # Load in data
                Base.delete!(file, "timeseries/$field/$iter") # Delete the entry
                # Write it again, adding the backward data using FieldTimeSeries interpolation. parent is used to strip offset from the backward data
                file["timeseries/$field/$iter"] = forward_data .+ parent(fts_backward[Time(T-forward_time)].data) 
            end
        end
    end
    println("Combined forward and backward contributions into $combined_output_filename")
end

function remove_halos(data,grid)
    Hx = grid.Hx
    Hy = grid.Hy
    Hz = grid.Hz
    data = data[
                Hx != 0 ? (Hx+1:end-Hx) : (:),
                Hy != 0 ? (Hy+1:end-Hy) : (:),
                Hz != 0 ? (Hz+1:end-Hz) : (:)]
    return data
end

function _create_coords(grid)
    full_grid_size = (grid.Nx + 2*grid.Hx, grid.Ny+ 2*grid.Hy, grid.Nz+ + 2*grid.Hz)
    x = xnodes(grid, Center(), with_halos = true)
    y = ynodes(grid, Center(), with_halos = true)
    z = znodes(grid, Center(), with_halos = true)


    all_coords = (x=x,y=y,z=z)
    coord_dict = Dict()
    coord_dict["full_grid_size"] = full_grid_size
    for (coord_name, coord) in pairs(all_coords)
        if coord !== nothing
            coord_dict[String(coord_name)] = parent(coord)
            if coord_name == :x
                coord_dict["x_mesh"] = parent(coord) .+ zeros(full_grid_size)
            elseif coord_name == :y
                coord_dict["y_mesh"] = parent(coord)' .+ zeros(full_grid_size)
            elseif coord_name == :z
                coord_dict["z_mesh"] = reshape(parent(coord),1,1,length(parent(coord))) .+ zeros(full_grid_size)
            end
        end
    end
    return coord_dict
end

function regrid_to_mean_position!(combined_output_filename, original_var_names, velocity_names, npad = 5)

    jldopen(combined_output_filename,"r+") do file
        iterations = parse.(Int, keys(file["timeseries/t"]))
        grid = file["serialized/grid"]
        coord_dict = _create_coords(grid)
        # Work out the periodic directions
        test_var = original_var_names[1]
        BCs = file["timeseries/$test_var/serialized/boundary_conditions"]
        periodic_dimensions = []
        if BCs.west == PeriodicBoundaryCondition()
            push!(periodic_dimensions,"x")
        end
        if BCs.south == PeriodicBoundaryCondition()
            push!(periodic_dimensions,"y")
        end
        if BCs.top == PeriodicBoundaryCondition()
            push!(periodic_dimensions,"z")
        end
        
        # First add the necessary serialized entry for each new variable
        for var in original_var_names 
            new_path = "timeseries/$var"*"_filtered_regrid/serialized"
            if haskey(file, new_path)
                Base.delete!(file, new_path) #incase we already tried to write this
            end
            g = Group(file, new_path)
            for property in keys(file["timeseries/$var/serialized"])
                g[property] = file["timeseries/$var/serialized/$property"]
            end
        end
        for iter in ProgressBar(iterations)
            Xi_list = []
            regular_coord_mesh = []
            n_true_dims = 0
            for dim in ("x","y","z")
                if dim in keys(coord_dict) # This limits to only the non-singleton dimensions
                    n_true_dims +=1
                    push!(regular_coord_mesh, coord_dict["$(dim)_mesh"])
                    if (dim == "x") && ("u" in velocity_names)
                        Xi_u = coord_dict["x_mesh"] .+ file["timeseries/xi_u/$iter"]
                        # Lose the halo regions (they don't help with fixed boundaries as they're zero, or with periodic as its repeated information)
                        Xi_u = remove_halos(Xi_u,grid)
                        # move xi points outside of domain into domain
                        if "x" in periodic_dimensions                           
                            Xi_u .-= floor.((Xi_u .- coord_dict["x"][grid.Hx+1])./grid.Lx) .* grid.Lx
                        end
                        push!(Xi_list,vec(Xi_u))


                    elseif (dim == "x") && !("u" in velocity_names)
                        Xi_u = coord_dict["x_mesh"]
                        # Lose the halo regions 
                        Xi_u = remove_halos(Xi_u,grid)
                        push!(Xi_list,vec(Xi_u))

                    elseif (dim == "y") && ("v" in velocity_names)
                        Xi_v =  coord_dict["y_mesh"] .+ file["timeseries/xi_v/$iter"]
                        # Lose the halo regions 
                        Xi_v = remove_halos(Xi_v,grid)
                        # move xi points outside of domain + halo regions into domain
                        if "y" in periodic_dimensions
                            Xi_v .-= floor.((Xi_v .- coord_dict["y"][grid.Hy+1])./grid.Ly) .* grid.Ly
                        end
                        push!(Xi_list,vec(Xi_v))

                    elseif (dim == "y") && !("v" in velocity_names)
                        Xi_v = coord_dict["y_mesh"]
                        # Lose the halo regions 
                        Xi_v = remove_halos(Xi_v,grid)
                        push!(Xi_list,vec(Xi_v))

                    elseif (dim == "z") && ("w" in velocity_names)
                        Xi_w = coord_dict["z_mesh"] .+ file["timeseries/xi_w/$iter"]
                        # Lose the halo regions 
                        Xi_w = remove_halos(Xi_w,grid)
                        # move xi points outside of domain + halo regions into domain
                        if "z" in periodic_dimensions
                            Xi_w .-= floor.((Xi_w .- coord_dict["z"][grid.Hz+1])./grid.Lz) .* grid.Lz
                        end
                        push!(Xi_list,vec(Xi_w))

                    elseif (dim == "z") && !("w" in velocity_names)
                        Xi_w = coord_dict["z_mesh"] .+ zeros(coord_dict["original_size"])
                        # Lose the halo regions 
                        Xi_w = remove_halos(Xi_w,grid)
                        push!(Xi_list,vec(Xi_w))

                    else
                        error("Something's wrong")
                    end
                end
            
            end
            # Now we do some padding on the periodic dimensions, introducing new elements to the list near the periodic boundaries
            # First construct a matrix that contains the coordinates and the fields to interpolate
            for var in original_var_names
                var_data = file["timeseries/$var"*"_filtered/$iter"]
                # Lose the halo regions 
                var_data = remove_halos(var_data,grid)

                push!(Xi_list, vec(var_data))
            end

            data_tuple = Tuple(Xi_list)
            data_array = hcat(data_tuple...) # This is a matrix where the rows are the data points and the columns are the coordinates then the variables

            # Then we take the array and repeat rows as necessary to add extra padding data
            column_number = 1
            n_columns = size(data_array,2)

            for dim in periodic_dimensions
                if dim == "x"
                    max_x = coord_dict["x"][end-grid.Hx]
                    min_x = coord_dict["x"][grid.Hx+1] 
                    xpad = npad*grid.Δxᶜᵃᵃ
                    Xi_x_vec = data_array[:,column_number] 
                    mask_max_x = (Xi_x_vec .> max_x - xpad)
                    mask_min_x = (Xi_x_vec .< min_x + xpad)
                    Xi_x_to_repeat = Xi_x_vec[mask_max_x .| mask_min_x]
                    # Remove or add Lx
                    Xi_x_to_repeat[(Xi_x_to_repeat .> max_x - xpad)] .-= grid.Lx
                    Xi_x_to_repeat[(Xi_x_to_repeat .< min_x + xpad)] .+= grid.Lx
                    extra_padding = fill(NaN, (length(Xi_x_to_repeat), n_columns))
                    extra_padding[:,column_number] = Xi_x_to_repeat
                    # And fill in the rest of the columns with straightforward repeateded data
                    for i in 1:n_columns
                        if i != column_number # Don't overwrite the x coordinate
                            extra_padding[:,i] = data_array[mask_max_x .| mask_min_x,i]
                        end
                    end
                    # Now join it on to the data array
                    data_array = vcat(data_array, extra_padding)
                    # Move along the rows to the next coordinate
                    column_number += 1

                elseif dim == "y"
                    max_y = coord_dict["y"][end-grid.Hy]
                    min_y = coord_dict["y"][grid.Hy+1]
                    ypad = npad*grid.Δyᵃᶜᵃ
                    Xi_y_vec = data_array[:,column_number]
                    mask_max_y = (Xi_y_vec .> max_y - ypad)
                    mask_min_y = (Xi_y_vec .< min_y + ypad)
                    Xi_y_to_repeat = Xi_y_vec[mask_max_y .| mask_min_y] 
                    # Remove or add Ly
                    Xi_y_to_repeat[(Xi_y_to_repeat .> max_y - ypad)] .-= grid.Ly
                    Xi_y_to_repeat[(Xi_y_to_repeat .< min_y + ypad)] .+= grid.Ly

                    extra_padding = fill(NaN, (length(Xi_y_to_repeat), n_columns))
                    extra_padding[:,column_number] = Xi_y_to_repeat
                    # And fill in the rest of the columns with straightforward repeateded data
                    for i in 1:n_columns
                        if i != column_number # Don't overwrite the y coordinate
                            extra_padding[:,i] = data_array[mask_max_y .| mask_min_y,i]
                        end
                    end
                    # Now join it on to the data array
                    data_array = vcat(data_array, extra_padding)
                    # Move along to the next coordinate
                    column_number += 1
                elseif dim == "z"
                    max_z = coord_dict["z"][end-grid.Hz]
                    min_z = coord_dict["z"][grid.Hz+1]
                    zpad = npad*grid.Δz.cᵃᵃᶜ
                    Xi_z_vec = data_array[:,column_number]
                    mask_max_z = ((Xi_z_vec .> max_z - zpad))
                    mask_min_z = ((Xi_z_vec .< min_z + zpad))
                    Xi_z_to_repeat = Xi_z_vec[mask_max_z .| mask_min_z] 

                    # Remove or add Lz
                    Xi_z_to_repeat[(Xi_z_to_repeat .> max_z - zpad)] .-= grid.Lz
                    Xi_z_to_repeat[(Xi_z_to_repeat .< min_z + zpad)] .+= grid.Lz

                    extra_padding = fill(NaN, (length(Xi_z_to_repeat), n_columns))
                    extra_padding[:,column_number] = Xi_z_to_repeat
                    # And fill in the rest of the columns with straightforward repeateded data
                    for i in 1:n_columns
                        if i != column_number # Don't overwrite the z coordinate
                            extra_padding[:,i] = data_array[mask_max_z .| mask_min_z,i]
                        end
                    end
                    # Now join it on to the data array
                    data_array = vcat(data_array, extra_padding)

                end
            end
            
            coords = data_array[:,1:n_true_dims] # The first columns are the coordinates
            var_data = data_array[:,(n_true_dims+1):end] # The rest is the data
            
            for (ivar,var) in enumerate(original_var_names)
                
                values = var_data[:,ivar]
                interpolator = scipy_interpolate.LinearNDInterpolator(coords, values)
                interp_data = pyconvert(Array,interpolator(regular_coord_mesh...))
                new_var_loc = "timeseries/$var"*"_filtered_regrid/$iter"
                if haskey(file, new_var_loc)
                    Base.delete!(file, new_var_loc) #incase we already tried to write this variable
                end
                file[new_var_loc] = interp_data
                
            end
        end
    end
end

# function regrid_to_mean_position_GeoStats!(combined_output_filename, original_var_names, velocity_names, interpolation_model=IDW(), maxneighbors=10, npad = 10)
#     # use velocity names to make this good in multiple dims
#     jldopen(combined_output_filename,"r+") do file
#         iterations = parse.(Int, keys(file["timeseries/t"]))
#         original_grid = file["serialized/grid"]
#         dim_dict, new_grid = _create_generalized_regular_grid(original_grid)
#         Hx, Hy, Hz = original_grid.Hx, original_grid.Hy, original_grid.Hz
        
#         # Work out the periodic directions
#         test_var = original_var_names[1]
#         BCs = file["timeseries/$test_var/serialized/boundary_conditions"]
#         periodic_dimensions = []
#         if BCs.west == PeriodicBoundaryCondition()
#             push!(periodic_dimensions,"x")
#         end
#         if BCs.south == PeriodicBoundaryCondition()
#             push!(periodic_dimensions,"y")
#         end
#         if BCs.top == PeriodicBoundaryCondition()
#             push!(periodic_dimensions,"z")
#         end
        
#         # First add the necessary serialized entry for each new variable
#         for var in original_var_names 
#             new_path = "timeseries/$var"*"_filtered_regrid/serialized"
#             if haskey(file, new_path)
#                 Base.delete!(file, new_path) #incase we already tried to write this
#             end
#             g = Group(file, new_path)
#             for property in keys(file["timeseries/$var/serialized"])
#                 g[property] = file["timeseries/$var/serialized/$property"]
#             end
#         end
        
#         for iter in ProgressBar(iterations)
            
#             Xi_list = []
#             n_true_dims = 0
#             for dim in ("x","y","z")
#                 if dim in keys(dim_dict) # This limits to only the non-singleton dimensions
#                     n_true_dims +=1
#                     if (dim == "x") && ("u" in velocity_names)
#                         Xi_u = dim_dict["x"] .+ file["timeseries/xi_u/$iter"]
#                         # Lose the halo regions (they don't help with fixed boundaries as they're zero, or with periodic as its repeated information)
#                         Xi_u = Xi_u[
#                             Hx != 0 ? (Hx+1:end-Hx) : (:),
#                             Hy != 0 ? (Hy+1:end-Hy) : (:),
#                             Hz != 0 ? (Hz+1:end-Hz) : (:)]

#                         # move xi points outside of domain into domain
#                         if "x" in periodic_dimensions                           
#                             Xi_u .-= floor.((Xi_u .- dim_dict["x"][Hx+1])./original_grid.Lx) .* original_grid.Lx
#                         end
#                         push!(Xi_list,vec(Xi_u))


#                     elseif (dim == "x") && !("u" in velocity_names)
#                         Xi_u = dim_dict["x"] .+ zeros(dim_dict["original_size"])
#                         # Lose the halo regions 
#                         Xi_u = Xi_u[
#                             Hx != 0 ? (Hx+1:end-Hx) : (:),
#                             Hy != 0 ? (Hy+1:end-Hy) : (:),
#                             Hz != 0 ? (Hz+1:end-Hz) : (:)]
#                         push!(Xi_list,vec(Xi_u))
#                     elseif (dim == "y") && ("v" in velocity_names)
                        
#                         Xi_v =  dim_dict["y"]' .+ file["timeseries/xi_v/$iter"]
#                         # Lose the halo regions 
#                         Xi_v = Xi_v[
#                             Hx != 0 ? (Hx+1:end-Hx) : (:),
#                             Hy != 0 ? (Hy+1:end-Hy) : (:),
#                             Hz != 0 ? (Hz+1:end-Hz) : (:)]
#                         # move xi points outside of domain + halo regions into domain
#                         if "y" in periodic_dimensions                    
#                             Xi_v .-= floor.((Xi_v .- dim_dict["y"][Hy+1])./original_grid.Ly) .* original_grid.Ly
#                         end
#                         push!(Xi_list,vec(Xi_v))
#                     elseif (dim == "y") && !("v" in velocity_names)
#                         Xi_v = dim_dict["y"]' .+ zeros(dim_dict["original_size"])
#                         # Lose the halo regions 
#                         Xi_v = Xi_v[
#                             Hx != 0 ? (Hx+1:end-Hx) : (:),
#                             Hy != 0 ? (Hy+1:end-Hy) : (:),
#                             Hz != 0 ? (Hz+1:end-Hz) : (:)]
#                         push!(Xi_list,vec(Xi_v))
#                     elseif (dim == "z") && ("w" in velocity_names)
#                         Xi_w = reshape(dim_dict["z"],1,1,length(dim_dict["z"])) .+ file["timeseries/xi_w/$iter"]
#                         # Lose the halo regions 
#                         Xi_w = Xi_w[
#                             Hx != 0 ? (Hx+1:end-Hx) : (:),
#                             Hy != 0 ? (Hy+1:end-Hy) : (:),
#                             Hz != 0 ? (Hz+1:end-Hz) : (:)]
#                         # move xi points outside of domain + halo regions into domain
#                         if "z" in periodic_dimensions
#                             Xi_w .-= floor.((Xi_w .- dim_dict["z"][Hz+1])./original_grid.Lz) .* original_grid.Lz
#                         end
#                         push!(Xi_list,vec(Xi_w))
#                     elseif (dim == "z") && !("w" in velocity_names)
#                         Xi_w = reshape(dim_dict["z"],1,1,length(dim_dict["z"])) .+ zeros(dim_dict["original_size"])
#                         # Lose the halo regions 
#                         Xi_w = Xi_w[
#                             Hx != 0 ? (Hx+1:end-Hx) : (:),
#                             Hy != 0 ? (Hy+1:end-Hy) : (:),
#                             Hz != 0 ? (Hz+1:end-Hz) : (:)]
#                         push!(Xi_list,vec(Xi_w))
#                     else
#                         error("Something's wrong")
#                     end
#                 end
            
#             end

#             # Now we do some padding on the periodic dimensions, introducing new elements to the list near the periodic boundaries
#             # First construct a matrix that contains the coordinates and the fields to interpolate
#             for var in original_var_names
#                 var_data = file["timeseries/$var"*"_filtered/$iter"]
#                 # Lose the halo regions 
#                 var_data = var_data[
#                     Hx != 0 ? (Hx+1:end-Hx) : (:),
#                     Hy != 0 ? (Hy+1:end-Hy) : (:),
#                     Hz != 0 ? (Hz+1:end-Hz) : (:)]
#                 push!(Xi_list, vec(var_data))
#             end

#             data_tuple = Tuple(Xi_list)
#             data_array = hcat(data_tuple...) # This is a matrix where the rows are the data points and the columns are the coordinates then the variables

#             # Then we take the array and repeat rows as necessary to add extra padding data
            
#             column_number = 1
#             n_columns = size(data_array,2)

#             for dim in periodic_dimensions
#                 if dim == "x"
#                     max_x = dim_dict["x"][end-Hx]
#                     min_x = dim_dict["x"][Hx+1] 
#                     xpad = npad*original_grid.Δxᶜᵃᵃ
#                     Xi_x_vec = data_array[:,column_number] 
#                     mask_max_x = ((Xi_x_vec .< max_x ) .& (Xi_x_vec .> max_x - xpad))
#                     mask_min_x = ((Xi_x_vec .> min_x ) .& (Xi_x_vec .< min_x + xpad))
#                     Xi_x_to_repeat = Xi_x_vec[mask_max_x .| mask_min_x]
#                     # Remove or add Lx
#                     Xi_x_to_repeat[(Xi_x_to_repeat.< max_x).& (Xi_x_to_repeat .> max_x - xpad)] .-= original_grid.Lx
#                     Xi_x_to_repeat[(Xi_x_to_repeat .> min_x ) .& (Xi_x_to_repeat .< min_x + xpad)] .+= original_grid.Lx

#                     extra_padding = fill(NaN, (length(Xi_x_to_repeat), n_columns))
#                     extra_padding[:,column_number] = Xi_x_to_repeat
#                     # And fill in the rest of the columns with straightforward repeateded data
#                     for i in 1:n_columns
#                         if i != column_number # Don't overwrite the x coordinate
#                             extra_padding[:,i] = data_array[mask_max_x .| mask_min_x,i]
#                         end
#                     end
#                     # Now join it on to the data array
#                     data_array = vcat(data_array, extra_padding)
#                     # Move along the rows to the next coordinate
#                     column_number += 1
#                 elseif dim == "y"
#                     max_y = dim_dict["y"][end-Hy]
#                     min_y = dim_dict["y"][Hy+1]
#                     ypad = npad*original_grid.Δyᵃᶜᵃ
#                     Xi_y_vec = data_array[:,column_number]
#                     mask_max_y = ((Xi_y_vec .< max_y ) .& (Xi_y_vec .> max_y - ypad))
#                     mask_min_y = ((Xi_y_vec .> min_y ) .& (Xi_y_vec .< min_y + ypad))
#                     Xi_y_to_repeat = Xi_y_vec[mask_max_y .| mask_min_y] 

#                     # Remove or add Ly
#                     Xi_y_to_repeat[(Xi_y_to_repeat.< max_y).& (Xi_y_to_repeat .> max_y - ypad)] .-= original_grid.Ly
#                     Xi_y_to_repeat[(Xi_y_to_repeat .> min_y ) .& (Xi_y_to_repeat .< min_y + ypad)] .+= original_grid.Ly

#                     extra_padding = fill(NaN, (length(Xi_y_to_repeat), n_columns))
#                     extra_padding[:,column_number] = Xi_y_to_repeat
#                     # And fill in the rest of the columns with straightforward repeateded data
#                     for i in 1:n_columns
#                         if i != column_number # Don't overwrite the y coordinate
#                             extra_padding[:,i] = data_array[mask_max_y .| mask_min_y,i]
#                         end
#                     end
#                     # Now join it on to the data array
#                     data_array = vcat(data_array, extra_padding)
#                     # Move along to the next coordinate
#                     column_number += 1
#                 elseif dim == "z"
#                     max_z = dim_dict["z"][end-Hz]
#                     min_z = dim_dict["z"][Hz+1]
#                     zpad = npad*original_grid.Δz.cᵃᵃᶜ   
#                     Xi_z_vec = data_array[:,column_number]
#                     mask_max_z = ((Xi_z_vec .< max_z ) .& (Xi_z_vec .> max_z - zpad))
#                     mask_min_z = ((Xi_z_vec .> min_z ) .& (Xi_z_vec .< min_z + zpad))
#                     Xi_z_to_repeat = Xi_z_vec[mask_max_z .| mask_min_z] 

#                     # Remove or add Lz
#                     Xi_z_to_repeat[(Xi_z_to_repeat.< max_z).& (Xi_z_to_repeat .> max_z - zpad)] .-= original_grid.Lz
#                     Xi_z_to_repeat[(Xi_z_to_repeat .> min_z ) .& (Xi_z_to_repeat .< min_z + zpad)] .+= original_grid.Lz

#                     extra_padding = fill(NaN, (length(Xi_z_to_repeat), n_columns))
#                     extra_padding[:,column_number] = Xi_z_to_repeat
#                     # And fill in the rest of the columns with straightforward repeateded data
#                     for i in 1:n_columns
#                         if i != column_number # Don't overwrite the z coordinate
#                             extra_padding[:,i] = data_array[mask_max_z .| mask_min_z,i]
#                         end
#                     end
#                     # Now join it on to the data array
#                     data_array = vcat(data_array, extra_padding)

#                 end
#             end
            
             
#             # We use GeoStats interpolation in Julia, which isn't as good as scipy's 
#             coords = Tuple.(eachrow(data_array[:,1:n_true_dims])) # The first columns are the coordinates
#             var_data = data_array[:,(n_true_dims+1):end] # The rest is the data
            
#             for (ivar,var) in enumerate(original_var_names)
                
#                 table_var = (; var=var_data[:,ivar])
#                 geotable = georef(table_var, coords)
#                 interp = geotable |> InterpolateNeighbors(new_grid, model=interpolation_model,maxneighbors=maxneighbors)
#                 interp_data = reshape(interp.var,(dim_dict["original_size"]...))
#                 new_var_loc = "timeseries/$var"*"_filtered_regrid/$iter"
#                 if haskey(file, new_var_loc)
#                     Base.delete!(file, new_var_loc) #incase we already tried to write this variable
#                 end
#                 file[new_var_loc] = interp_data
                
#             end
#         end
#     end
# end


# Maybe we should keep this one elsewhere
using DataStructures: OrderedDict
using NCDatasets
function jld2_to_netcdf(jld2_filename,nc_filename)
    jldopen(jld2_filename, "r") do file
        
        iterations = parse.(Int, keys(file["timeseries/t"]))
        times = [file["timeseries/t/$iter"] for iter in iterations]
        dt = times[2] - times[1]
        grid = file["serialized/grid"]

        rm(nc_filename, force=true)
        ds = NCDataset(nc_filename,"c", attrib = OrderedDict(
            "Julia"                     => "This file was generated using Julia Version $VERSION",

        ))

        # Dimensions
        ds.dim["time"] = length(times)
        ds.dim["y_afa"] = try; length(grid.yᵃᶠᵃ) catch; 1 end
        ds.dim["x_faa"] = try; length(grid.xᶠᵃᵃ) catch; 1 end
        ds.dim["x_caa"] = try; length(grid.xᶜᵃᵃ) catch;1 end
        ds.dim["y_aca"] = try; length(grid.yᵃᶜᵃ) catch; 1 end
        ds.dim["z_aaf"] = try; length(grid.z.cᵃᵃᶠ) catch; 1 end
        ds.dim["z_aac"] = try; length(grid.z.cᵃᵃᶜ) catch; 1 end


        # Declare grid variables

        nctime = defVar(ds,"time", Float64, ("time",), attrib = OrderedDict(
            "units"                     => "seconds",
            "long_name"                 => "Time",
        ))
        nctime[:] = times

        ncy_afa = defVar(ds,"y_afa", Float32, ("y_afa",), attrib = OrderedDict(
            "units"                     => "m",
            "long_name"                 => "Cell face locations in the y-direction.",
        ))
        ncy_afa[:] = try; parent(grid.yᵃᶠᵃ) catch; 0 end

        ncx_faa = defVar(ds,"x_faa", Float32, ("x_faa",), attrib = OrderedDict(
            "units"                     => "m",
            "long_name"                 => "Cell face locations in the x-direction.",
        ))
        ncx_faa[:] = try; parent(grid.xᶠᵃᵃ) catch; 0 end

        ncx_caa = defVar(ds,"x_caa", Float32, ("x_caa",), attrib = OrderedDict(
            "units"                     => "m",
            "long_name"                 => "Cell center locations in the x-direction.",
        ))
        ncx_caa[:] = try; parent(grid.xᶜᵃᵃ) catch; 0 end

        ncy_aca = defVar(ds,"y_aca", Float32, ("y_aca",), attrib = OrderedDict(
            "units"                     => "m",
            "long_name"                 => "Cell center locations in the y-direction.",
        ))
        ncy_aca[:] = try; parent(grid.yᵃᶜᵃ) catch; 0 end

        ncz_aac = defVar(ds,"z_aac", Float32, ("z_aac",), attrib = OrderedDict(
            "units"                     => "m",
            "long_name"                 => "Cell center locations in the z-direction.",
        ))
        ncz_aac[:] = try; parent(grid.z.cᵃᵃᶜ) catch; 0 end

        ncz_aaf = defVar(ds,"z_aaf", Float32, ("z_aaf",), attrib = OrderedDict(
            "units"                     => "m",
            "long_name"                 => "Cell face locations in the z-direction.",
        ))
        ncz_aaf[:] = try; parent(grid.z.cᵃᵃᶠ) catch; 0 end

        ncdx_caa = defVar(ds,"dx_caa", Float32, (), attrib = OrderedDict(
            "units"                     => "m",
            "long_name"                 => "Spacing between cell faces (located at the cell centers) in the x-direction.",
        ))
        ncdx_caa[:] = try; parent(grid.Δxᶜᵃᵃ) catch; 0 end

        ncdx_faa = defVar(ds,"dx_faa", Float32, (), attrib = OrderedDict(
            "units"                     => "m",
            "long_name"                 => "Spacing between cell centers (located at the cell faces) in the x-direction.",
        ))
        ncdx_faa[:] = try; parent(grid.Δxᶠᵃᵃ) catch; 0 end

        ncdy_aca = defVar(ds,"dy_aca", Float32, (), attrib = OrderedDict(
            "units"                     => "m",
            "long_name"                 => "Spacing between cell faces (located at cell centers) in the y-direction.",
        ))
        ncdy_aca[:] = try; parent(grid.Δyᵃᶜᵃ) catch; 0 end

        ncdy_afa = defVar(ds,"dy_afa", Float32, (), attrib = OrderedDict(
            "units"                     => "m",
            "long_name"                 => "Spacing between cell centers (located at cell faces) in the y-direction.",
        ))
        ncdy_afa[:] = try; parent(grid.Δyᵃᶠᵃ) catch; 0 end

        ncdz_aac = defVar(ds,"dz_aac", Float32, (), attrib = OrderedDict(
            "units"                     => "m",
            "long_name"                 => "Spacing between cell faces (located at cell centers) in the z-direction.",
        ))
        ncdz_aac[:] = try; parent(grid.Δz.cᵃᵃᶜ) catch; 0 end

        ncdz_aaf = defVar(ds,"dz_aaf", Float32, (), attrib = OrderedDict(
            "units"                     => "m",
            "long_name"                 => "Spacing between cell centers (located at cell faces) in the z-direction.",
        ))
        ncdz_aaf[:] = try; parent(grid.Δz.cᵃᵃᶠ) catch; 0 end

        # Define dimensions
        ncNx = defVar(ds,"Nx", Int64, (), attrib = OrderedDict(
            "units"                     => "None",
            "long_name"                 => "Number of cells in the x-direction.",
        ))
        ncNx[:] = grid.Nx

        ncNy = defVar(ds,"Ny", Int64, (), attrib = OrderedDict(
            "units"                     => "None",
            "long_name"                 => "Number of cells in the y-direction.",
        ))
        ncNy[:] = grid.Ny


        ncNz = defVar(ds,"Nz", Int64, (), attrib = OrderedDict(
            "units"                     => "None",
            "long_name"                 => "Number of cells in the z-direction.",
        ))
        ncNz[:] = grid.Nz
        # Define halos

        ncHx = defVar(ds,"Hx", Int64, (), attrib = OrderedDict(
            "units"                     => "None",
            "long_name"                 => "Number of halo cells in the x-direction.",
        ))
        ncHx[:] = grid.Hx

        ncHy = defVar(ds,"Hy", Int64, (), attrib = OrderedDict(
            "units"                     => "None",
            "long_name"                 => "Number of halo cells in the y-direction.",
        ))
        ncHy[:] = grid.Hy

        ncHz = defVar(ds,"Hz", Int64, (), attrib = OrderedDict(
            "units"                     => "None",
            "long_name"                 => "Number of halo cells in the z-direction.",
        ))
        ncHz[:] = grid.Hz

        # Define dimension lengths

        ncLx = defVar(ds,"Lx", Float64, (), attrib = OrderedDict(
            "units"                     => "m",
            "long_name"                 => "Length of the domain in the x-direction.",
        ))
        ncLx[:] = grid.Lx
        ncLy = defVar(ds,"Ly", Float64, (), attrib = OrderedDict(
            "units"                     => "m",
            "long_name"                 => "Length of the domain in the y-direction.",
        ))
        ncLy[:] = grid.Ly   
        ncLz = defVar(ds,"Lz", Float64, (), attrib = OrderedDict(
            "units"                     => "m",
            "long_name"                 => "Length of the domain in the z-direction.",
        ))
        ncLz[:] = grid.Lz

        # And define variables 
        location_map = Dict(
            (Center,Center,Center) => ("x_caa","y_aca","z_aac"),
            (Face,Center,Center)   => ("x_faa","y_aca","z_aac"),
            (Center,Face,Center)   => ("x_caa","y_afa","z_aac"),
            (Face,Face,Center)     => ("x_faa","y_afa","z_aac"),
            (Center,Center,Face)   => ("x_caa","y_aca","z_aaf"),
            (Face,Center,Face)     => ("x_faa","y_aca","z_aaf"),
            (Center,Face,Face)     => ("x_caa","y_afa","z_aaf"),
            (Face,Face,Face)       => ("x_faa","y_afa","z_aaf"),
        )
        

        for varname in keys(file["timeseries"])
            if varname != "t"
                    
                # Create a variable in the NetCDF file
                bc_string = sprint(show, file["timeseries/xi_u/serialized/boundary_conditions"])
                ncv = defVar(ds, varname, Float64, (location_map[file["timeseries/$varname/serialized/location"]]...,"time"), attrib = OrderedDict(
                    "boundary conditions"                     => bc_string,
                ))
                for (it,iter) in enumerate(iterations)
                    # Assign data to the variable
                    ncv[:, :, :, it] = file["timeseries/$varname/$iter"]
                end
                
            end
        end

        close(ds)
        
    end
end
# TODO 
# Look at boundaries, 3D, etc in interpolation
# Add in a single exponential as a filter
# How can BCs be implemented? Might need fill_halo_regions. especially for a closed boundary where U neq 0
# Keyword args