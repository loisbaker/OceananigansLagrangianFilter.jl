using JLD2
using JLD2: Group
using Oceananigans
using Oceananigans.Fields: Center
using Oceananigans.BoundaryConditions: PeriodicBoundaryCondition
using Oceananigans.Units: Time
using DataStructures: OrderedDict
using NCDatasets
using ProgressBars

# Python import to use the LinearNDInterpolator for regridding
using PythonCall
const scipy_interpolate = PythonCall.pynew()
function __init__()
    PythonCall.pycopy!(scipy_interpolate, pyimport("scipy.interpolate"))
end


"""
    copy_file_metadata!(original_file::JLD2.JLDFile, new_file::JLD2.JLDFile, timeseries_vars_to_copy::Tuple{Vararg{String}})

Copies Oceananigans output JLD2 file metadata to a new file. 

Copies variables and their serialized entries as defined in `timeseries_vars_to_copy`, but does not copy any timeseries data, or times.

"""
function copy_file_metadata!(original_file::JLD2.JLDFile, new_file::JLD2.JLDFile, timeseries_vars_to_copy::Tuple{Vararg{String}})

    # Copy over metadata that isn't associated with variables
    for path in ("grid","serialized","coriolis")
        _copy_jld2_recursive!(original_file, new_file, path)
    end

    # Copy serialized entries for timeseries variables
    for var in timeseries_vars_to_copy
        _copy_jld2_recursive!(original_file, new_file, "timeseries/$var/serialized")
    end

end



"""
    _copy_jld2_recursive!(source::JLD2.JLDFile, dest::JLD2.JLDFile, path::String)

A helper function to recursively copy data and groups.
"""
function _copy_jld2_recursive!(source::JLD2.JLDFile, dest::JLD2.JLDFile, path::String)
    
    # Check if the path points to a JLD2 group
    if isa(source[path], JLD2.Group)
        # Recreate the group in the destination file
        dest_group = JLD2.Group(dest, path)

        # Recursively copy the contents of the group
        for key in keys(source[path])
            _copy_jld2_recursive!(source, dest, "$path/$key")
        end
    else
        # Simply copy the variable if it's not a group
        dest[path] = source[path]
    end
end

"""
    create_input_data_on_disk(original_data_filename::String, var_names_to_filter::Tuple{Vararg{String}}, velocity_names::Tuple{Vararg{String}}; direction::String="forward", T_start::Union{Real, Nothing}=nothing, T_end::Union{Real, Nothing}=nothing)

Creates a new temporary JLD2 file for filtering named <original_data_filename>_filter_input.jld2 with reorganised and reorientated data

This function is designed to handle simulation data stored on disk and prepare it for time-stepping.

It performs four main tasks:
1. Deletes any existing version of this file
2. Creates a new file with the original file's metadata, and adds `direction` and `t_simulation` as extra variables.
3. Changes the time coordinate `t` to be relative to `T_start` or `T_end` depending on `direction`.
4. Reverses the order and/or sign of data fields from the original data and saves to the new file.


# Arguments
- `original_data_filename`: The path to the original JLD2 file to be copied.
- `var_names_to_filter`: A tuple of tracer variable names to copy from the original data.
- `velocity_names`: A tuple of velocity variable names to copy from the original data.

# Keyword Arguments
- `direction::String`: The desired direction of the time-series. Can be `"forward"` or `"backward"`.
  Defaults to `"forward"`.
- `T_start::Union{Real, Nothing}`: The start time for the filtered period. If `nothing`, it defaults to the
  first time step in the simulation for `"forward"` direction, or the last time step for `"backward"`.
- `T_end::Union{Real, Nothing}`: The end time for the filtered period. If `nothing`, it defaults to the
  last time step in the simulation for `"forward"` direction, or the first time step for `"backward"`.

# Returns
- `T_filter`: The total duration of the filtered period.

# Examples
```julia
# Create a file named "my_data_filter_input.jld2" to run in the backward direction
create_input_data_on_disk("my_data.jld2", ("var1", "var2"), ("vel1", "vel2"); direction="backward")

# Reformat "my_data_filter_input.jld2" for a forward simulation from time 10 to 20
create_input_data_on_disk("my_data.jld2", ("var1", "var2"), ("vel1", "vel2"); direction="forward", T_start=10, T_end=20)
```
"""
function create_input_data_on_disk(original_data_filename::String, var_names_to_filter::Tuple{Vararg{String}}, velocity_names::Tuple{Vararg{String}}; direction::String="forward", T_start::Union{Real, Nothing}=nothing, T_end::Union{Real, Nothing}=nothing)
    
    # Check that a valid direction has been given
    if direction ∉ ("backward", "forward")
        error("Invalid direction: $direction. Must be 'backward' or 'forward'.")
    end

    # Check that the original file exists 
    if !isfile(original_data_filename)
        error("Source file not found: $original_data_filename")
    end

    # Create a new filename and delete any file that already has that name
    new_filename = original_data_filename[1:end-5]*"_filter_input.jld2"
    if isfile(new_filename)
        rm(new_filename)
    end

    # Open the original file for reading 
    jldopen(original_data_filename,"r") do original_file

        # Check if velocities and tracers given are in the original_file
        for var in (var_names_to_filter..., velocity_names...)
            if !haskey(original_file["timeseries"], var)
                error("Variable '$var' not found in original data file.")
            end
        end

        # Check if T_start and T_end are found in the original_file
        iterations = parse.(Int, keys(original_file["timeseries/t"]))
        times = [original_file["timeseries/t/$iter"] for iter in iterations]
        
        # If T_start is not specified, use the first simulation time
        T_start = isnothing(T_start) ? times[1] : T_start

        # If T_end is not specified, use the last simulation time
        T_end = isnothing(T_end) ? times[end] : T_end

        if T_start < times[1] || T_start > times[end]
            error("T_start=$T_start is outside the range of the original data: [$times[1], $times[end]].")
        end

        if T_end < times[1] || T_end > times[end]
            error("T_end=$T_end is outside the range of the original data: [$times[1], $times[end]].")
        end

        # Create a truncated iterations variable with just the times to copy
        iterations_truncated = iterations[(times .>= T_start) .& (times .<= T_end)]

        # Create a new filename for the filtered input data
        jldopen(new_filename, "w") do new_file

            # We need to copy the standard metadata from the old file to the new file
            copy_file_metadata!(original_file, new_file, (var_names_to_filter..., velocity_names...))

            # Add an indicator of the direction 
            new_file["direction"] = direction

            # Add an old (t_simulation) and a new (t) time variable 
            t_sim_group = JLD2.Group(new_file, "timeseries/t_simulation")
            t_group = JLD2.Group(new_file, "timeseries/t")

            if direction == "forward"

                # First write times, starting from T_start
                for iter in iterations_truncated
                    t_sim_group["$iter"] = original_file["timeseries/t/$iter"]
                    t_group["$iter"] = original_file["timeseries/t/$iter"] - T_start
                end

                # Then write in data, as in original file
                for var in (var_names_to_filter..., velocity_names...)
                    for iter in iterations_truncated
                        new_file["timeseries/$var/$iter"] = original_file["timeseries/$var/$iter"]
                    end
                end

            elseif direction == "backward"

                # First write times, starting from T_end
                for iter in reverse(iterations_truncated)
                    t_sim_group["$iter"] = original_file["timeseries/t/$iter"]
                    t_group["$iter"] = T_end - original_file["timeseries/t/$iter"]
                end

                # Then write in data (as in original file for tracers and negated for velocities)
                for var in var_names_to_filter
                    for iter in reverse(iterations_truncated)
                        new_file["timeseries/$var/$iter"] = original_file["timeseries/$var/$iter"]
                    end
                end

                for var in velocity_names
                    for iter in reverse(iterations_truncated)
                        new_file["timeseries/$var/$iter"] = -original_file["timeseries/$var/$iter"]
                    end
                end

                
            end

            # Calculate and return the total filter time
            T_filter = T_end - T_start
            return T_filter
        end
    end
end

function load_data(original_data_filename, var_names_to_filter, velocity_names; architecture=CPU(), backend= InMemory())

    input_data_filename = original_data_filename[1:end-5]*"_filter_input.jld2"
    velocity_timeseries = ()
    for vel_name in velocity_names
        vel_ts = FieldTimeSeries(input_data_filename, vel_name; architecture=architecture, backend=backend)
        velocity_timeseries = (velocity_timeseries...,vel_ts)
    end
   
    var_timeseries = ()
    for var_name in var_names_to_filter
        var_ts = FieldTimeSeries(input_data_filename, var_name; architecture=architecture, backend=backend)
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

function create_original_vars(var_names_to_filter, grid)
    # Creates auxiliary fields to store the saved variables
    vars = Dict()
    for var_name in var_names_to_filter
        vars[Symbol(var_name)] = CenterField(grid)
    end
    return NamedTuple(vars)
end

function create_filtered_vars(var_names_to_filter, velocity_names, filter_params; map_to_mean=true)
    N_coeffs = Int(length(filter_params)/4)
    gC = ()  # Start with an empty tuple
    gS = ()  # Start with an empty tuple
    for var_name in var_names_to_filter
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

function create_forcing(filtered_vars, var_names_to_filter, velocity_names, filter_params)

    N_coeffs = Int(length(filter_params)/4)

    # Initialize dictionary
    gCdict = Dict()
    gSdict = Dict()

    # Make forcing function for original data term. Final arg is the field dependence.
    original_var_forcing_func(args...) = args[end]

    for var_name in var_names_to_filter
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
    if length(filtered_vars) > N_coeffs*length(var_names_to_filter)*2
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
function create_output_fields(model, var_names_to_filter, velocity_names, filter_params)
    N_coeffs = Int(length(filter_params)/4)
    N_filtered_vars = length(model.tracers)
    half_N_filtered_vars = Int(N_filtered_vars/2)
    outputs_dict = Dict()
    for (i_var, original_var) in enumerate(var_names_to_filter)

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
    if N_filtered_vars > N_coeffs*length(var_names_to_filter)*2
        for (i_vel, vel) in enumerate(velocity_names)
            i_var = i_vel + length(var_names_to_filter)
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
    for original_var in var_names_to_filter
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


function sum_forward_backward_contributions!(combined_output_filename,forward_output_filename,backward_output_filename,T,velocity_names, var_names_to_filter)
# Combine the forward and backward simulations by summing them into a single file

    # Make a copy of the forward file to fill in with combined data
    cp(forward_output_filename, combined_output_filename,force=true)

    # List the names of the fields that we will combine
    data_field_names = vcat(["xi_" * vel for vel in velocity_names], [var * "_filtered" for var in var_names_to_filter])

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

function regrid_to_mean_position!(combined_output_filename, var_names_to_filter, velocity_names, npad = 5)

    jldopen(combined_output_filename,"r+") do file
        iterations = parse.(Int, keys(file["timeseries/t"]))
        grid = file["serialized/grid"]
        coord_dict = _create_coords(grid)
        # Work out the periodic directions
        test_var = var_names_to_filter[1]
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
        for var in var_names_to_filter 
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
                            Xi_u .-= floor.((Xi_u .- grid.xᶠᵃᵃ[1])./grid.Lx) .* grid.Lx
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
                            Xi_v .-= floor.((Xi_v .- grid.yᵃᶠᵃ[1])./grid.Ly) .* grid.Ly
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
                            Xi_w .-= floor.((Xi_w .- grid.z.cᵃᵃᶠ[1])./grid.Lz) .* grid.Lz
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
            for var in var_names_to_filter
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
                    max_x = grid.xᶠᵃᵃ[grid.Nx+1]
                    min_x = grid.xᶠᵃᵃ[1]
                    xpad = npad*grid.Δxᶜᵃᵃ
                    Xi_x_vec = data_array[:,column_number] 
                    mask_max_x = (Xi_x_vec .> max_x - xpad) .& (Xi_x_vec .< max_x)
                    mask_min_x = (Xi_x_vec .< min_x + xpad) .& (Xi_x_vec .> min_x)
                    Xi_x_to_repeat = Xi_x_vec[mask_max_x .| mask_min_x]
                    # Remove or add Lx
                    Xi_x_to_repeat[(Xi_x_to_repeat .> max_x - xpad) .& (Xi_x_to_repeat .< max_x)] .-= grid.Lx
                    Xi_x_to_repeat[(Xi_x_to_repeat .< min_x + xpad) .& (Xi_x_to_repeat .> min_x)] .+= grid.Lx
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
                    max_y = grid.yᵃᶠᵃ[grid.Ny+1]
                    min_y = grid.yᵃᶠᵃ[1]
                    ypad = npad*grid.Δyᵃᶜᵃ
                    Xi_y_vec = data_array[:,column_number]
                    mask_max_y = (Xi_y_vec .> max_y - ypad) .& (Xi_y_vec .< max_y )
                    mask_min_y = (Xi_y_vec .< min_y + ypad) .& (Xi_y_vec .> min_y)   
                    Xi_y_to_repeat = Xi_y_vec[mask_max_y .| mask_min_y] 

                    # Remove or add Ly
                    Xi_y_to_repeat[(Xi_y_to_repeat .> max_y - ypad) .& (Xi_y_to_repeat .< max_y)] .-= grid.Ly
                    Xi_y_to_repeat[(Xi_y_to_repeat .< min_y + ypad) .& (Xi_y_to_repeat .> min_y)] .+= grid.Ly

                    
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
                    max_z = grid.z.cᵃᵃᶠ[grid.Nz+1]
                    min_z = grid.z.cᵃᵃᶠ[1]
                    zpad = npad*grid.Δz.cᵃᵃᶜ
                    Xi_z_vec = data_array[:,column_number]
                    mask_max_z = (Xi_z_vec .> max_z - zpad) .& (Xi_z_vec .< max_z)
                    mask_min_z = (Xi_z_vec .< min_z + zpad) .& (Xi_z_vec .> min_z)
                    Xi_z_to_repeat = Xi_z_vec[mask_max_z .| mask_min_z]

                    # Remove or add Lz
                    Xi_z_to_repeat[(Xi_z_to_repeat .> max_z - zpad) .& (Xi_z_to_repeat .< max_z)] .-= grid.Lz
                    Xi_z_to_repeat[(Xi_z_to_repeat .< min_z + zpad) .& (Xi_z_to_repeat .> min_z)] .+= grid.Lz

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
            
            for (ivar,var) in enumerate(var_names_to_filter)
                
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

function get_weight_function(t,tstar,N,freq_c)
    
    N_coeffs = 2^(N-1)
    G = zeros(t)
    for i in 1:N_coeffs
        
        a = (freq_c/2^N)*sin(pi/(2^(N+1))*(2*i-1))
        b = (freq_c/2^N)*cos(pi/(2^(N+1))*(2*i-1))
        c = freq_c*sin(pi/(2^(N+1))*(2*i-1))
        d = freq_c*cos(pi/(2^(N+1))*(2*i-1))

        G += (a*cos(d*abs(t-tstar)) + b*sin(d*abs(t-tstar)))*exp(-c*abs(t-tstar))
    end

    return G
end