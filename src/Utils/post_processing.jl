# Post-processing of the filter output: combining the forward and backward passes, regridding to the mean
# position, and the Eulerian filter and time shift for comparison

const scipy_interpolate = PythonCall.pynew()
const numpy = PythonCall.pynew()
function __init__()
    PythonCall.pycopy!(scipy_interpolate, pyimport("scipy.interpolate"))
    PythonCall.pycopy!(numpy, pyimport("numpy"))
end

"""
    sum_forward_backward_contributions!(config::AbstractConfig)

Combines the output from the forward and backward filter simulations into a
single output file. This function performs the final step of the offline
filter algorithm by summing the contributions from each pass.

The function performs the following steps:
1.  **Initializes the combined file**: A new JLD2 file is created to store
    the final output.
2.  **Copies metadata and unfiltered data**: The file structure, metadata,
    and the original, unfiltered data are copied from the forward output file.
3.  **Sums filtered contributions**: For each filtered variable, the data
    from the backward output file is loaded as a `FieldTimeSeries`. The data
    is then interpolated to match the time steps of the forward simulation,
    and the two datasets are summed and written to the combined output file.

The outputs of each pass are that pass's contribution to the filtered fields (see
[`create_output_fields`](@ref)), so all filtered fields, including the mean velocities, are summed.

Arguments
=========
- `config`: An instance of `AbstractConfig` containing the file paths and
  variable names.
- `extra_filtered_var_names::Tuple{Vararg{String}}=()`: Optional tuple of additional
  filtered variable names that have been calculated by the filter and also need to be
  combined. The outputs of each pass should be defined as its contribution, so that they
  are summed (e.g. negating quantities that change sign under time reversal in the
  backward pass, as `create_output_fields` does for the mean velocities).
- `extra_original_data_names::Tuple{Vararg{String}}=()`: Optional tuple of additional
  original names that have been output and should be copied to the combined output file.
"""
function sum_forward_backward_contributions!(config::AbstractConfig; extra_filtered_var_names::Tuple{Vararg{String}}=(),
    extra_original_data_names::Tuple{Vararg{String}}=())
    # Combine the forward and backward simulations by summing them into a single file

    output_filename = config.output_filename
    forward_output_filename = config.forward_output_filename
    backward_output_filename = config.backward_output_filename
    T = config.T
    velocity_names = config.velocity_names
    var_names_to_filter = config.var_names_to_filter
    compute_mean_velocities = config.compute_mean_velocities
    
    # The names of the filtered fields output by each pass (see create_output_fields), which we will combine
    filtered_var_names = filtered_output_names(config)

    # The original velocities are copied if the mean velocities were computed
    vel_names_to_filter = compute_mean_velocities ? velocity_names : ()

    # There might be some extra filtered variables that the user defined that we should combine too
    filtered_var_names = Tuple(unique((filtered_var_names..., extra_filtered_var_names...)))

    jldopen(output_filename,"w") do combined_file
        jldopen(forward_output_filename,"r") do forward_file

            # It's possible the file doesn't contain some of the filtered variables - for example the user might have defined their own outputs. 
            # Let's check and only try to copy the variables that exist in the forward file
            forward_file_all_names = keys(forward_file["timeseries"])
            missing_var_names = [var for var in filtered_var_names if !(var in forward_file_all_names)]
            
            if length(missing_var_names) > 0
                @warn "The following filtered variable names were not found in the forward output file and will be skipped: $(missing_var_names)"
            end

            filtered_var_names = Tuple([var for var in filtered_var_names if var in forward_file_all_names])


            # First copy the forward file metadata and file structure
            if config.output_original_data
                names_to_copy = (var_names_to_filter..., vel_names_to_filter..., filtered_var_names...)
            else
                names_to_copy = filtered_var_names
            end

            # There might be some extra variables provided to copy too
            names_to_copy = Tuple(unique((names_to_copy..., extra_original_data_names...)))
            copy_file_metadata!(forward_file, combined_file, names_to_copy)

            forward_iterations = parse.(Int, keys(forward_file["timeseries/t"]))

            # Copy over the unfiltered field data

            if config.output_original_data
                original_data_names = (var_names_to_filter..., vel_names_to_filter..., extra_original_data_names...,"t")
            else
                original_data_names = ("t",extra_original_data_names...)
            end

            
            for var_name in original_data_names
                for iter in forward_iterations
                    combined_file["timeseries/$var_name/$iter"] = forward_file["timeseries/$var_name/$iter"]
                end
            end
            

            # Copy over the filtered data, combined with backward data
            for var_name in filtered_var_names
                
                # Open the backward data as a FieldTimeSeries, so we can interpolate to match times
                fts_backward = FieldTimeSeries(backward_output_filename, var_name)

                # Loop over forward times and add the backward data
                for iter in forward_iterations
                    forward_time = forward_file["timeseries/t/$iter"]
                    forward_data = forward_file["timeseries/$var_name/$iter"] # Load in data
                    
                    # Write it again, adding the backward data using FieldTimeSeries interpolation. parent is used to strip offset from the backward data
                    combined_file["timeseries/$var_name/$iter"] = forward_data .+ parent(fts_backward[Time(T-forward_time)].data)
                end
            end

        end
    end
    @info "Combined forward and backward contributions into $output_filename"

end

"""
    _remove_halos(data::AbstractArray, grid::AbstractGrid)

Removes the halo regions from a 3D data array based on the halo sizes specified
in the `grid` object.

Arguments
=========
- `data`: An `AbstractArray` representing the 3D data with halo regions.
- `grid`: An object containing the halo sizes `Hx`, `Hy`, and `Hz`.

Returns
=======
- A view of the `data` array with the halo regions removed.
"""
function _remove_halos(data::AbstractArray, grid::AbstractGrid)
    Hx = grid.Hx
    Hy = grid.Hy
    Hz = grid.Hz
    data = data[
                Hx != 0 ? (Hx+1:end-Hx) : (:),
                Hy != 0 ? (Hy+1:end-Hy) : (:),
                Hz != 0 ? (Hz+1:end-Hz) : (:)]
    return data
end

"""
    _fill_halos!(data::AbstractArray, grid::AbstractGrid, value::Real=0.0)

Fills the halo regions of a 3D array `data` with a specified `value`.

This function modifies the `data` array in-place. The dimensions of the halo regions
are determined by the `Hx`, `Hy`, and `Hz` fields of the `grid` object. 

# Arguments
- `data::AbstractArray`: The 3D array whose halo regions will be filled.
- `grid::AbstractGrid`: An object containing the dimensions of the halo regions (`Hx`, `Hy`, `Hz`).
- `value::Real=0.0`: The scalar value to use for filling the halo regions. Defaults to `0.0`.

# Returns
- `nothing`: This function does not return a value. It modifies the input array directly.

"""
function _fill_halos!(data::AbstractArray, grid::AbstractGrid, value::Real=0.0)
    Hx = grid.Hx
    Hy = grid.Hy
    Hz = grid.Hz
    data[1:Hx,:,:] .= value
    data[end-Hx+1:end,:,:] .= value
    data[:,1:Hy,:] .= value
    data[:,end-Hy+1:end,:] .= value
    data[:,:,1:Hz] .= value
    data[:,:,end-Hz+1:end] .= value
    return nothing
end

"""
    _create_coords(grid)

Creates a dictionary of coordinate arrays and mesh grids for a given `grid` object.

# Arguments
- `grid`: A grid object that defines the spatial dimensions and halo sizes.
  It must have fields `Nx`, `Ny`, `Nz`, `Hx`, `Hy`, and `Hz`.

# Returns
- `Dict`: A dictionary with the following keys:
    - `"full_grid_size"`: A tuple representing the dimensions of the full grid
      `(Nx + 2*Hx, Ny + 2*Hy, Nz + 2*Hz)`.
    - `"x"`, `"y"`, `"z"`: 1D `AbstractArray`s containing the `Center` coordinates along each axis,
      including halos.
    - `"x_mesh"`, `"y_mesh"`, `"z_mesh"`: 3D `AbstractArray`s representing the `Center` coordinate
      values at every point in the full grid. These are created by broadcasting the 1D
      coordinates to the full grid size.
"""
function _create_coords(grid::AbstractGrid)
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
                coord_dict["x_mesh"] = reshape(parent(coord),length(parent(coord)),1,1) .+ zeros(full_grid_size)
            elseif coord_name == :y
                coord_dict["y_mesh"] = reshape(parent(coord),1,length(parent(coord)),1) .+ zeros(full_grid_size)
            elseif coord_name == :z
                coord_dict["z_mesh"] = reshape(parent(coord),1,1,length(parent(coord))) .+ zeros(full_grid_size)
            end
        end
    end
    return coord_dict
end

"""
    _mask_immersed(data::AbstractArray, grid::AbstractGrid, Immersed::Bool, buffer_dz::Int = 3)
Applies masking to the input `data` array based on the presence of immersed boundaries in the `grid`. 
If `Immersed` is `true`, the function masks out regions below the bottom height of the immersed boundary, 
adding a buffer zone of `buffer_dz`*`dz` to avoid interpolation issues. RectilinearGrid is assumed.

Arguments
=========
- `data`: An instance of `AbstractArray` containing the array to be masked
- `grid`: An ImmersedBoundaryGrid grid object that defines the spatial dimensions and halo sizes.
- `data_has_halos`: A `Bool` indicating whether the input `data` includes halo regions. Defaults to `true`.
- `buffer_dz`: An `Int` specifying the number of grid cells to use as a buffer zone above the bottom height. Defaults to `3`.
"""
function _mask_immersed(data::AbstractArray, grid::ImmersedBoundaryGrid; fill_value::Float64 = NaN, data_has_halos::Bool = true, buffer_dz::Int = 0)
    # If there are immersed boundaries, we'll mask out the solid regions plus a small buffer zone with NaN (for field being regridded)
    # or 0 (for maps) to avoid bad interpolations
    if !data_has_halos # Halos have been removed already
        bottom_height = parent(grid.immersed_boundary.bottom_height)[grid.Hx+1:end-grid.Hx,grid.Hy+1:end-grid.Hy,:] # Remove halos
        z3D = _remove_halos(_create_coords(grid)["z_mesh"], grid)
        data[z3D .< bottom_height .+ buffer_dz*abs(grid.z.Δᵃᵃᶜ)] .= fill_value
    else # Data still has halos
        bottom_height = parent(grid.immersed_boundary.bottom_height)
        z3D = _create_coords(grid)["z_mesh"]
        data[z3D .< bottom_height .+ buffer_dz*abs(grid.z.Δᵃᵃᶜ)] .= fill_value
    end
    return data
end


"""
    regrid_to_mean_position!(config::AbstractConfig)

Regrids the filtered data to the mean position. This function reads the combined 
output file, interpolates the filtered variables to the mean position, and saves the
result in new variables within the same file.

The regridding process involves the following steps:
1.  **Extracts positions**: The mean positions (`xi_u`, `xi_v`, `xi_w`) and
    filtered variable data are extracted for each time step.
2.  **Handles periodicity**: For periodic dimensions (x, y, or z), the data is
    padded by repeating values near the boundaries to ensure accurate
    interpolation across the periodic boundaries.
3.  **Interpolates data**: A linear interpolator is used to map the filtered data
    from the irregular advected positions to the original, regular grid points.
4.  **Saves new fields**: The regridded data is saved as new variables in the
    combined output file, with a `_Lagrangian_filtered_at_mean` suffix.

Arguments
=========
- `config`: An instance of `AbstractConfig` containing the file paths, variable
  names, and grid information.
- `extra_vars_to_regrid::Tuple{Vararg{String}}=()`: Optional tuple of additional
  filtered variable names that have been calculated by the filter and also need to be
  regridded. Include velocities here if needed.
"""
function regrid_to_mean_position!(config::AbstractConfig; extra_vars_to_regrid::Tuple{Vararg{String}}=())
    #TODO split this function into smaller functions for readability
    output_filename = config.output_filename
    var_names_to_filter = config.var_names_to_filter
    compute_mean_velocities = config.compute_mean_velocities
    velocity_names = config.velocity_names
    npad = config.npad 
    label = config.label


    if !config.lagrangian
        error("Regridding to mean position is not meaningful for Eulerian filtering")
    end
    if compute_mean_velocities
        var_names_to_filter = (var_names_to_filter..., velocity_names...)
    end
    
    # Get names of labelled mean variables
    var_names_to_regrid = Tuple([var * label * "_Lagrangian_filtered" for var in var_names_to_filter])

    # There might be some extra variables to regrid that the user defined
    var_names_to_regrid = Tuple(unique((var_names_to_regrid..., extra_vars_to_regrid...)))

    jldopen(output_filename,"r+") do file
        iterations = parse.(Int, keys(file["timeseries/t"]))
        grid = file["serialized/grid"]
        coord_dict = _create_coords(grid)

        # Check if it is an immersed boundary grid:
        if isa(grid, ImmersedBoundaryGrid)
            @warn "Grid is an ImmersedBoundaryGrid, regridding to mean position may not be sensible. 
            Areas below bottom height will be masked with NaNs, assuming bottom height is a function of x and/or y."
            Immersed = true
        else
            Immersed = false
        end

        # Error if not a RectilinearGrid
        if !(grid isa RectilinearGrid) && !((grid isa ImmersedBoundaryGrid) && (grid.underlying_grid isa RectilinearGrid))
            error("Regridding to mean position is only implemented for RectilinearGrid")
        end

        # Work out the periodic and bounded directions
        test_var = var_names_to_regrid[1]
        BCs = file["timeseries/$test_var/serialized/boundary_conditions"]
        periodic_dimensions = []
        bounded_dimensions = []
        if BCs.west == PeriodicBoundaryCondition()
            push!(periodic_dimensions,"x")
        elseif !Immersed && !isnothing(BCs.west)
            # We only do the fixed boundary correction if we know this is a non singleton dimension
            push!(bounded_dimensions,"x")
            @info "Assuming velocities normal to x boundaries are zero (open boundaries not yet supported for regridding)"
        end
        if BCs.south == PeriodicBoundaryCondition()
            push!(periodic_dimensions,"y")
        elseif !Immersed && !isnothing(BCs.south)   
            push!(bounded_dimensions,"y")
            @info "Assuming velocities normal to y boundaries are zero (open boundaries not yet supported for regridding)"
        end
        if BCs.top == PeriodicBoundaryCondition()
            push!(periodic_dimensions,"z")
        elseif !Immersed && !isnothing(BCs.top)
            push!(bounded_dimensions,"z")
            @info "Assuming velocities normal to z boundaries are zero (open boundaries not yet supported for regridding)"
        end
        
        # First add the necessary serialized entry for each new variable
        for var in var_names_to_regrid
            new_path = "timeseries/$var"*"_at_mean/serialized"
            if haskey(file, new_path)
                Base.delete!(file, new_path) #incase we already tried to write this
            end
            g = Group(file, new_path)
            for property in keys(file["timeseries/$var/serialized"])
                g[property] = file["timeseries/$var/serialized/$property"]
            end
        end

        # Now loop over time steps, extract the positions and data, and interpolate to the mean position
        for iter in iterations
            Xi_list = []
            regular_coord_mesh = []
            n_true_dims = 0
            true_dims = []
            for dim in ("x","y","z")
                if dim in keys(coord_dict) # This limits to only the non-singleton dimensions
                    n_true_dims +=1
                    push!(true_dims, dim)
                    push!(regular_coord_mesh, coord_dict["$(dim)_mesh"])
                    if (dim == "x") && ("u" in velocity_names)
                        if !Immersed
                            Xi_u = coord_dict["x_mesh"] .+ file["timeseries/xi_u"*label*"/$iter"]
                        else
                            Xi_u = coord_dict["x_mesh"] .+ _mask_immersed(file["timeseries/xi_u"*label*"/$iter"], grid, fill_value = 0.)
                        end
                        # Lose the halo regions (they don't help with fixed boundaries as they're zero, or with periodic as its repeated information)
                        Xi_u = _remove_halos(Xi_u,grid)

                        # move xi points outside of domain into domain
                        if "x" in periodic_dimensions                           
                            Xi_u .-= floor.((Xi_u .- grid.xᶠᵃᵃ[1])./grid.Lx) .* grid.Lx
                        end
                        push!(Xi_list,vec(Xi_u))
                        # Accompany this vector with a vector of x-indices for this data
                        # Get the size of the array
                        indices_array = [I[1] for I in CartesianIndices(size(Xi_u))]
                        push!(Xi_list,vec(indices_array))

                    elseif (dim == "x") && !("u" in velocity_names)
                        Xi_u = coord_dict["x_mesh"]

                        # Lose the halo regions 
                        Xi_u = _remove_halos(Xi_u,grid)
                        push!(Xi_list,vec(Xi_u))
                        indices_array = [I[1] for I in CartesianIndices(size(Xi_u))]
                        push!(Xi_list,vec(indices_array))

                    elseif (dim == "y") && ("v" in velocity_names)
                        if !Immersed
                            Xi_v =  coord_dict["y_mesh"] .+ file["timeseries/xi_v"*label*"/$iter"]
                        else
                            Xi_v =  coord_dict["y_mesh"] .+ _mask_immersed(file["timeseries/xi_v"*label*"/$iter"], grid, fill_value = 0.)
                        end

                        # Lose the halo regions 
                        Xi_v = _remove_halos(Xi_v,grid)

                        # move xi points outside of domain + halo regions into domain
                        if "y" in periodic_dimensions
                            Xi_v .-= floor.((Xi_v .- grid.yᵃᶠᵃ[1])./grid.Ly) .* grid.Ly
                        end
                        push!(Xi_list,vec(Xi_v))
                        indices_array = [I[2] for I in CartesianIndices(size(Xi_v))]
                        push!(Xi_list,vec(indices_array))

                    elseif (dim == "y") && !("v" in velocity_names)
                        Xi_v = coord_dict["y_mesh"]

                        # Lose the halo regions 
                        Xi_v = _remove_halos(Xi_v,grid)
                        push!(Xi_list,vec(Xi_v))
                        indices_array = [I[2] for I in CartesianIndices(size(Xi_v))]
                        push!(Xi_list,vec(indices_array))

                    elseif (dim == "z") && ("w" in velocity_names)
                        if !Immersed
                            Xi_w = coord_dict["z_mesh"] .+ file["timeseries/xi_w"*label*"/$iter"]
                        else
                            Xi_w = coord_dict["z_mesh"] .+ _mask_immersed(file["timeseries/xi_w"*label*"/$iter"], grid, fill_value = 0.)
                        end

                        # Lose the halo regions 
                        Xi_w = _remove_halos(Xi_w,grid)

                        # move xi points outside of domain + halo regions into domain
                        if "z" in periodic_dimensions
                            Xi_w .-= floor.((Xi_w .- grid.z.cᵃᵃᶠ[1])./grid.Lz) .* grid.Lz
                        end
                        push!(Xi_list,vec(Xi_w))
                        indices_array = [I[3] for I in CartesianIndices(size(Xi_w))]
                        push!(Xi_list,vec(indices_array))

                    elseif (dim == "z") && !("w" in velocity_names)
                        Xi_w = coord_dict["z_mesh"] .+ zeros(coord_dict["original_size"])

                        # Lose the halo regions 
                        Xi_w = _remove_halos(Xi_w,grid)
                        push!(Xi_list,vec(Xi_w))
                        indices_array = [I[3] for I in CartesianIndices(size(Xi_w))]
                        push!(Xi_list,vec(indices_array))

                    else
                        error("Something's wrong with the dimensions of the regridding routine")
                    end
                end
            
            end
            # Now we do some padding on the periodic dimensions, introducing new elements to the list near the periodic boundaries
            # First construct a matrix that contains the coordinates and the fields to interpolate
            for var in var_names_to_regrid
                if !Immersed
                    var_data = file["timeseries/$var/$iter"]
                else
                    var_data = _mask_immersed(file["timeseries/$var/$iter"], grid, fill_value = NaN)
                end
                
                # Lose the halo regions 
                var_data = _remove_halos(var_data,grid)

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
                    column_number += 2 # We added an extra column of indices, so skip these

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
                    column_number += 2 # We added an extra column of indices, so skip these
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
            
            coords = data_array[:,1:2:2*n_true_dims] # The first columns are the coordinates and indices, just take the coordinates
            var_data = data_array[:,(2*n_true_dims+1):end] # The rest is the data
            
            # Normalisation to help with interpolation stability. Delaunay triangulation can struggle with very different scales in 
            # different dimensions (e.g. x - z slice of ocean).
            coord_mins = minimum(coords, dims=1)
            coord_maxs = maximum(coords, dims=1)
            coord_ranges = coord_maxs .- coord_mins
            coords_norm = (coords .- coord_mins) ./ coord_ranges
            
            regular_coord_mesh_norm = []
            for (i, mesh) in enumerate(regular_coord_mesh)
                norm_mesh = mesh .-coord_mins[i]
                norm_mesh = norm_mesh ./coord_ranges[i]
                push!(regular_coord_mesh_norm, norm_mesh)
            end

            for (ivar,var) in enumerate(var_names_to_regrid)    
                
                # This is the main interpolation
                values = var_data[:,ivar]
                interpolator = scipy_interpolate.LinearNDInterpolator(coords_norm, values)
                interp_data = pyconvert(Array,interpolator(regular_coord_mesh_norm...))

                # We already dealt with the periodic boundaries with padding, but we now make sure
                # that the interpolation is accurate at fixed boundaries by doing an interpolation 
                # at fixed coordinate normal to the boundary.
                # Reuse the data array using the indices columns to pull out the data we need. 
                for bounded_dim in bounded_dimensions # Just the non-singleton dimensions with fixed boundary
                    halo_size = getproperty(grid, Symbol("H$bounded_dim"))
                    true_dim_index = findfirst(isequal(bounded_dim), true_dims)
                    for edge_index in (1, getproperty(grid, Symbol("N$bounded_dim"))) # The index of the first and last real grid points in this dimension - halos have already been trimmed before indices assigned
                        edge_index_with_halos = edge_index + halo_size

                        # Filter to coordinates where the index in this dimension is edge_index
                        data_array_cut = data_array[Int.(data_array[:,2*true_dim_index]) .== Int(edge_index),:] 
                        
                        # Now remove the columns with the bounded coordinate
                        data_array_cut = hcat(data_array_cut[:, 1:2*true_dim_index-2], data_array_cut[:, 2*true_dim_index+1:end])
                        
                        # Then sort data_array_cut by the first remaining coordinate - this is needed for 1D interpolation
                        data_array_cut = sortslices(data_array_cut, dims=1, by=row -> row[1])

                        # data_array_cut can now be used for interpolation in one fewer dimension
                        coords_cut = data_array_cut[:,1:2:2*(n_true_dims-1)] # The first columns are the coordinates and indices, just take the coordinates
                        values_cut = data_array_cut[:,2*(n_true_dims-1)+ivar] # The rest is the data
                        
                        # Normalise coords_cut for interpolation stability
                        coord_mins_cut = minimum(coords_cut, dims=1)
                        coord_maxs_cut = maximum(coords_cut, dims=1)
                        coord_ranges_cut = coord_maxs_cut .- coord_mins_cut
                        coords_cut_norm = (coords_cut .- coord_mins_cut) ./ coord_ranges_cut
                        if bounded_dim == "x"
                            if n_true_dims == 2
                                # No need to normalise coords_cut as we are only doing 1D interpolation here
                                if "y" in true_dims
                                    mesh = coord_dict["y_mesh"][edge_index_with_halos,:,1]
                                    interp_data[edge_index_with_halos,:,1] = pyconvert(Array,numpy.interp(mesh, coords_cut[:,1], values_cut))
                                elseif "z" in true_dims
                                    mesh = coord_dict["z_mesh"][edge_index_with_halos,1,:]
                                    interp_data[edge_index_with_halos,1,:] = pyconvert(Array,numpy.interp(mesh, coords_cut[:,1], values_cut))
                                else
                                    @info "Dimension combo $true_dims is not implemented"
                                end
                                
                            elseif n_true_dims == 3
                                # Need to normalise for interpolation stability
                                y_mesh = coord_dict["y_mesh"][edge_index_with_halos,:,:]
                                z_mesh = coord_dict["z_mesh"][edge_index_with_halos,:,:]
                                y_mesh_norm = (y_mesh .- coord_mins_cut[2]) ./ coord_ranges_cut[2]
                                z_mesh_norm = (z_mesh .- coord_mins_cut[3]) ./ coord_ranges_cut[3]
                                meshes = (y_mesh_norm,z_mesh_norm)
                                interpolator = scipy_interpolate.LinearNDInterpolator(coords_cut_norm, values_cut)
                                interp_data[edge_index_with_halos,:,:] = pyconvert(Array,interpolator(meshes...))
                            else
                                @info "Number of dimensions $n_true_dims is not implemented"
                            end
            
                        elseif bounded_dim =="y"
                            if n_true_dims == 2
                                # No need to normalise coords_cut as we are only doing 1D interpolation here
                                if "x" in true_dims
                                    mesh = coord_dict["x_mesh"][:,edge_index_with_halos,1]
                                    interp_data[:,edge_index_with_halos,1] = pyconvert(Array,numpy.interp(mesh, coords_cut[:,1], values_cut))
                                elseif "z" in true_dims
                                    mesh = coord_dict["z_mesh"][1,edge_index_with_halos,:]
                                    interp_data[1,edge_index_with_halos,:] = pyconvert(Array,numpy.interp(mesh, coords_cut[:,1], values_cut))
                                else
                                    @info "Dimension combo $true_dims is not implemented"
                                end
                                
                                
                            elseif n_true_dims == 3
                                # Need to normalise for interpolation stability
                                x_mesh = coord_dict["x_mesh"][:,edge_index_with_halos,:]
                                z_mesh = coord_dict["z_mesh"][:,edge_index_with_halos,:]
                                x_mesh_norm = (x_mesh .- coord_mins_cut[1]) ./ coord_ranges_cut[1]
                                z_mesh_norm = (z_mesh .- coord_mins_cut[3]) ./ coord_ranges_cut[3]
                                meshes = (x_mesh_norm,z_mesh_norm)
                                interpolator = scipy_interpolate.LinearNDInterpolator(coords_cut_norm, values_cut)
                                interp_data[:,edge_index_with_halos,:] = pyconvert(Array,interpolator(meshes...))
                            else
                                @info "Number of dimensions $n_true_dims is not implemented"
                            end
                        elseif bounded_dim == "z"
                            if n_true_dims == 2
                                # No need to normalise coords_cut as we are only doing 1D interpolation here
                                if "x" in true_dims
                                    mesh = coord_dict["x_mesh"][:,1,edge_index_with_halos]
                                    interp_data[:,1,edge_index_with_halos] = pyconvert(Array,numpy.interp(mesh, coords_cut[:,1], values_cut))                           
                                elseif "y" in true_dims
                                    mesh = coord_dict["y_mesh"][1,:,edge_index_with_halos]
                                    interp_data[1,:,edge_index_with_halos] = pyconvert(Array,numpy.interp(mesh, coords_cut[:,1], values_cut))
                                else
                                    @info "Dimension combo $true_dims is not implemented"
                                end
                                
                            elseif n_true_dims == 3
                                # Need to normalise for interpolation stability
                                x_mesh = coord_dict["x_mesh"][:,:,edge_index_with_halos]
                                y_mesh = coord_dict["y_mesh"][:,:,edge_index_with_halos]
                                x_mesh_norm = (x_mesh .- coord_mins_cut[1]) ./ coord_ranges_cut[1]
                                y_mesh_norm = (y_mesh .- coord_mins_cut[2]) ./ coord_ranges_cut[2]
                                meshes = (x_mesh_norm,y_mesh_norm)
                                interpolator = scipy_interpolate.LinearNDInterpolator(coords_cut_norm, values_cut)
                                interp_data[:,:,edge_index_with_halos] = pyconvert(Array,interpolator(meshes...))
                            else
                                @info "Number of dimensions $n_true_dims is not implemented"
                            end
                        end
                    end 

                end
                
                # Now we make sure the halo regions are filled with zeros again
                _fill_halos!(interp_data, grid, 0.0)
                
                # Finally write the data to a new variable in the file
                new_var_loc = "timeseries/$var"*"_at_mean/$iter"
                if haskey(file, new_var_loc)
                    Base.delete!(file, new_var_loc) #incase we already tried to write this variable
                end
                file[new_var_loc] = interp_data
                
            end
        end
    end
    @info "Wrote regridded data to new variables with _at_mean suffix in file $output_filename"
    
end

# The terms of a weight function (see get_weight_function) as complex exponentials: at a time a lag tau before the
# output time, a term is Re[past * exp(-lambda * tau)], and at a time a lag tau after it, Re[future * exp(-lambda * tau)],
# with lambda = c - i d. The sine terms change sign after the output time if they're odd.
function exponential_terms(filter_params)
    if filter_params.N_coeffs == 0.5
        a, c = complex(filter_params.a1), complex(filter_params.c1)
        return [(; past = a, future = a, lambda = c)]
    end
    _, sine_parity = filter_output_options(filter_params)
    sign_after = sine_parity === :odd ? -1 : 1
    return map(1:filter_params.N_coeffs) do i
        a, b, c, d = (getproperty(filter_params, Symbol(coefficient, i)) for coefficient in ("a", "b", "c", "d"))
        (; past = complex(a, -b), future = complex(a, -sign_after * b), lambda = complex(c, -d))
    end
end

# With a frequency scale M that doesn't vary in time, the weight function at a fixed point is M G(M tau), so each term's
# amplitudes and decay rate are multiplied by M (an array of its values at the saved points)
scale_terms(terms, ::Nothing) = terms
scale_terms(terms, M) = [(; past = M .* term.past, future = M .* term.future, lambda = M .* term.lambda) for term in terms]

# The frequency scale at the points where var_name is saved in the output file (which may include halo points), or
# nothing for a uniform filter. It can only depend on the coordinates (not on frequency_scale_dependencies), so that it
# doesn't vary in time
function saved_frequency_scale(config::AbstractConfig, file, var_name, saved_size)
    isnothing(config.frequency_scale_func) && return nothing
    func, params = config.frequency_scale_func, config.frequency_scale_params
    grid = on_architecture(CPU(), load_serialized_grid(file, var_name))
    LX, LY, LZ = file["timeseries/$var_name/serialized/location"]
    M = Field{LX, LY, LZ}(grid; indices = file["timeseries/$var_name/serialized/indices"])
    for index in CartesianIndices(M.data)
        M.data[index] = func(node(Tuple(index)..., grid, LX(), LY(), LZ())..., params)
    end
    size(parent(M)) == saved_size && return Array(parent(M))
    size(interior(M)) == saved_size && return Array(interior(M))
    error("The saved $var_name, of size $saved_size, doesn't match its grid, of size $(size(parent(M))) with halos or " *
          "$(size(interior(M))) without")
end

# Decays the running sums of each term over the time interval dt, adds f * dt_weight to them, and returns the sum over
# the terms of Re[amplitude * running sum], where amplitude is the past or future amplitude of the term
function update_running_sums!(running_sums, terms, amplitude, f, dt_weight, dt)
    value = zeros(size(f))
    for (running_sum, term) in zip(running_sums, terms)
        running_sum .*= exp.(-term.lambda .* dt)
        running_sum .+= f .* dt_weight
        value .+= real.(getproperty(term, amplitude) .* running_sum)
    end
    return value
end

# The Eulerian filter as running sums, so that each field is read at most twice (rather than once for every output
# time). read(j) gives the input f_j for each weight function in terms_list (each a list of exponential terms), and for
# each output time t_i, output(i, values) is called with the sums over j of G(t_i - t_j) * f_j * dt_weights[j], one for
# each weight function. For direction = :both, the sums over
# later times are built up first, backwards in time, and kept in backward_store (anything indexed by a String, e.g. a
# Dict or a JLD2 file), and then the sums over earlier times, forwards in time, so that the outputs are written in order.
function Eulerian_running_sums(read, output, times, dt_weights, terms_list, direction, backward_store)
    n = length(times)
    new_running_sums(inputs) = [[zeros(ComplexF64, size(f)) for _ in terms] for (terms, f) in zip(terms_list, inputs)]

    if direction === :both
        running_sums = nothing
        for i in n:-1:1
            inputs = read(i)
            running_sums = something(running_sums, new_running_sums(inputs))
            dt = i < n ? times[i+1] - times[i] : 0.0
            for (w, (terms, sums, f)) in enumerate(zip(terms_list, running_sums, inputs))
                value = update_running_sums!(sums, terms, :future, f, dt_weights[i], dt)
                # The output time itself is counted in both the earlier and the later times, with G(0) = sum of a
                value .-= sum(real.(term.future) for term in terms) .* f .* dt_weights[i]
                backward_store["$i/$w"] = value
            end
        end
    end

    running_sums = nothing
    for i in 1:n
        inputs = read(i)
        running_sums = something(running_sums, new_running_sums(inputs))
        dt = i > 1 ? times[i] - times[i-1] : 0.0
        values = map(enumerate(zip(terms_list, running_sums, inputs))) do (w, (terms, sums, f))
            value = update_running_sums!(sums, terms, :past, f, dt_weights[i], dt)
            direction === :both ? value .+ backward_store["$i/$w"] : value
        end
        output(i, values)
    end
    return nothing
end

"""
    compute_Eulerian_filter!(config::AbstractConfig)

Computes the Eulerian filter for specified variables and writes the results to a
combined output file. At each grid point, the saved fields are filtered in time with
the filter's weight function. Each term of the weight function is a decaying complex
exponential, so the filtered fields at all the output times are built up as running
sums: forwards in time over the earlier times, and offline also backwards over the later
times. Each saved field is read once in each direction, so the cost grows linearly with
the number of saved times. Offline, the sums over the later times are kept in a temporary
file next to the output file. The results are saved to the same JLD2 file in a new group
with the `_Eulerian_filtered` suffix.

This method serves as a benchmark for comparison with the main Lagrangian filter.

If the filter outputs its terms separately (`outputs = :separate`), each term is Eulerian filtered with its own
weight function, and named like the Lagrangian outputs with `_Eulerian_filtered` in the stem, e.g.
`b_Eulerian_filtered_C1_scaled` and `b_Eulerian_filtered_S1_scaled`. The outputs are only divided by the integral
of the weight function (correcting for it being cut off at the ends of the data) if the whole filter is normalised,
in which case every term is divided by the same integral, so that the terms add up to the combined output.

With a frequency scale ``M`` that only depends on the coordinates (no `frequency_scale_dependencies`), the weight
function at each point is ``M G(M \\tau)``, with ``M`` evaluated where each variable is saved (e.g. at the `u` points
for `u`). A frequency scale with `frequency_scale_dependencies` varies in time, which isn't supported here: for the
Eulerian filter with that frequency scale, run the filter with `lagrangian = false` instead.

The method uses data saved to the filter output file - incase we decide to save this
at lower frequency than the original data, it should be rewritten to use the original
data file instead.

Arguments
=========
- `config`: An instance of `AbstractConfig` containing the file path, variable
  names, and filter parameters.
"""
function compute_Eulerian_filter!(config::AbstractConfig)
    isempty(config.frequency_scale_dependencies) ||
        error("compute_Eulerian_filter! doesn't support a frequency scale with frequency_scale_dependencies, which varies " *
              "in time. For the Eulerian filter with the same frequency scale, run the filter without advection " *
              "(lagrangian = false)")
    filter_params = config.filter_params
    outputs, _ = filter_output_options(filter_params)
    output_filename = config.output_filename
    var_names_to_filter = config.var_names_to_filter
    compute_mean_velocities = config.compute_mean_velocities
    velocity_names = config.velocity_names
    label = config.label

    var_names_to_Eulerian_filter = var_names_to_filter
    if compute_mean_velocities
        var_names_to_Eulerian_filter = (var_names_to_Eulerian_filter..., velocity_names...)
    end
    
    if (config isa AbstractOfflineConfig)
        direction = :both
    elseif (config isa AbstractOnlineConfig)
        direction = :forward
    end

    # Dividing by the discrete integral of the weight function corrects for it being cut off at the ends of the
    # data, which only makes sense for a normalised filter, whether its terms are output combined or separately
    # (for a band-pass filter, it can be close to zero)
    normalised = zero_frequency_gain(filter_params, direction) ≈ 1

    # The weight functions to filter with: the whole filter, or each of its terms if they are output separately
    weights = outputs === :combined ? [("", filter_params)] : separate_term_params(filter_params)

    # Open existing file
    jldopen(output_filename,"r+") do file
        iterations = parse.(Int, keys(file["timeseries/t"]))
        times = [file["timeseries/t/$iter"] for iter in iterations]
        n = length(times)

        # Each field is weighted by the time interval to the next one (the last by the one before it)
        dt_weights = [diff(times); times[n] - times[n-1]]


        # Loop over variables to filter
        for var_name in var_names_to_Eulerian_filter
            @info "Computing Eulerian filter for variable $var_name"

            # Create a group for each filtered output, copying over the serialized properties
            groups = map(weights) do (suffix, _)
                name = var_name * label * "_Eulerian_filtered" * suffix
                g_EF = Group(file, "timeseries/$name")
                g_EF_serialized = Group(file, "timeseries/$name/serialized")
                for property in keys(file["timeseries/$var_name/serialized"])
                    g_EF_serialized[property] = file["timeseries/$var_name/serialized/$property"]
                end
                g_EF
            end

            # Every output is normalised by the discrete integral of the whole filter's weight function, so the terms
            # add up to the whole. It's found with the same running sums, as the filter of one at each time except the
            # last (which it leaves out, unlike the filtered fields), and depends on position with a frequency scale
            M = saved_frequency_scale(config, file, var_name, size(file["timeseries/$(var_name)/$(iterations[1])"]))
            terms_list = [[scale_terms(exponential_terms(p), M) for (_, p) in weights]; [scale_terms(exponential_terms(filter_params), M)]]
            ones_shape = isnothing(M) ? () : size(M)
            function read_inputs(i)
                field = file["timeseries/$(var_name)/$(iterations[i])"]
                return [fill(field, length(weights)); [fill(i < n ? 1.0 : 0.0, ones_shape)]]
            end
            function write_outputs(i, values)
                normalisation = values[end]
                for (g_EF, value) in zip(groups, values[1:end-1])
                    g_EF["$(iterations[i])"] = normalised ? value ./ normalisation : value
                end
            end

            # The offline filter keeps the sums over later times in a temporary file next to the output
            if direction === :both
                backward_filename = output_filename * ".Eulerian_backward_sums.jld2"
                try
                    jldopen(backward_filename, "w") do backward_store
                        Eulerian_running_sums(read_inputs, write_outputs, times, dt_weights, terms_list, direction, backward_store)
                    end
                finally
                    rm(backward_filename, force = true)
                end
            else
                Eulerian_running_sums(read_inputs, write_outputs, times, dt_weights, terms_list, direction, nothing)
            end
        end
    end
end

"""
    compute_time_shift!(config::AbstractConfig)

Computes the time shift for a forward filter (either online or forward-only offline)
based on its coefficients and writes the shifted time series to the output file.

The time shift is computed as the time delay introduced by the filter's transfer function. 
This new time series is stored in a new group called `timeseries/t_shifted` within the output JLD2 file.

# Arguments
- `config`: A configuration object of type `AbstractFilterConfig` which contains
  the `output_filename` and `filter_params` (filter coefficients).

"""
function compute_time_shift!(config::AbstractConfig)
    isnothing(config.frequency_scale_func) ||
        error("The time shift is not defined for a filter with a frequency scale, since the centre of its weight function varies")
    if !(config isa AbstractOnlineConfig)
        @warn "Time shift computation is only relevant when filtering forward only. The offline filter uses both the past " *
              "and the future, so no time shift is needed. This function will compute a time shift regardless, but " *
              "it may not be meaningful for the offline filter."
    end
    output_filename = config.output_filename
    filter_params = config.filter_params
    direction = config isa AbstractOnlineConfig ? :forward : :both
    if !(zero_frequency_gain(filter_params, direction) ≈ 1)
        @warn "The filter is not normalised, so the time shift computed here (the centre of its weight function) " *
              "is not meaningful. For a spectral filter, changes in the amplitude of the signal at each frequency " *
              "are delayed by about 1/alpha (so t - 1/alpha is the time for the envelope or power), while its phase " *
              "is only shifted by about alpha/(2 omega_n) radians."
    end
    N_coeffs = filter_params.N_coeffs
    time_shift = 0.0
    if N_coeffs == 0.5 # exponential special case
        time_shift = 1/filter_params.c1
    else
        for i in 1:N_coeffs
            a = getproperty(filter_params, Symbol("a$i"))
            b = getproperty(filter_params, Symbol("b$i"))
            c = getproperty(filter_params, Symbol("c$i"))
            d = getproperty(filter_params, Symbol("d$i"))
            time_shift += (a*c^2 + 2*b*c*d - a*d^2)/(c^2 + d^2)^2
        end
    end
    jldopen(output_filename,"r+") do file
        iterations = parse.(Int, keys(file["timeseries/t"]))
        times = [file["timeseries/t/$iter"] for iter in iterations]
        t_shift_group = JLD2.Group(file, "timeseries/t_shifted")

        for (i, t) in enumerate(times)
            t_shift_group["$(iterations[i])"] = t - time_shift
        end

    end
    @info "Wrote time shift data to new group timeseries/t_shifted in file $output_filename"

end
