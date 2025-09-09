using Oceananigans.BoundaryConditions: PeriodicBoundaryCondition
using DataStructures: OrderedDict
using Oceananigans.Grids: AbstractGrid
using NCDatasets


# Python import to use the LinearNDInterpolator for regridding
using PythonCall
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

Arguments
=========
- `config`: An instance of `AbstractConfig` containing the file paths and
  variable names.
"""
function sum_forward_backward_contributions!(config::AbstractConfig)
    # Combine the forward and backward simulations by summing them into a single file

    output_filename = config.output_filename
    forward_output_filename = config.forward_output_filename
    backward_output_filename = config.backward_output_filename
    T = config.T
    velocity_names = config.velocity_names
    var_names_to_filter = config.var_names_to_filter

    # List the names of the fields that we will combine
    filtered_var_names = vcat(["xi_" * vel for vel in velocity_names], [var * "_Lagrangian_filtered" for var in var_names_to_filter])
    
    jldopen(output_filename,"w") do combined_file
        jldopen(forward_output_filename,"r") do forward_file

            # First copy the forward file metadata and file structure
            copy_file_metadata!(forward_file, combined_file, (var_names_to_filter..., filtered_var_names...))

            forward_iterations = parse.(Int, keys(forward_file["timeseries/t"]))

            # Copy over the unfiltered field data
            for var_name in (var_names_to_filter...,"t")
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
    - `"x"`, `"y"`, `"z"`: 1D `AbstractArray`s containing the coordinates along each axis,
      including halos.
    - `"x_mesh"`, `"y_mesh"`, `"z_mesh"`: 3D `AbstractArray`s representing the coordinate
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
"""
#TODO split this function into smaller functions for readability
function regrid_to_mean_position!(config::AbstractConfig)
    output_filename = config.output_filename
    var_names_to_filter = config.var_names_to_filter
    velocity_names = config.velocity_names
    npad = config.npad 
    
    jldopen(output_filename,"r+") do file
        iterations = parse.(Int, keys(file["timeseries/t"]))
        grid = file["serialized/grid"]
        coord_dict = _create_coords(grid)

        # Check if it is an immersed boundary grid:
        if isa(grid, ImmersedBoundaryGrid)
            @warn "Grid is an ImmersedBoundaryGrid, regridding to mean position may not be sensible"
            Immersed = true
        else
            Immersed = false
        end

        # Work out the periodic and bounded directions
        test_var = var_names_to_filter[1]
        BCs = file["timeseries/$test_var/serialized/boundary_conditions"]
        periodic_dimensions = []
        bounded_dimensions = []
        if BCs.west == PeriodicBoundaryCondition()
            push!(periodic_dimensions,"x")
        elseif Immersed == false && !isnothing(BCs.west)
            # We only do the fixed boundary correction if we know there's no immersed boundary and this is a non singleton dimension
            push!(bounded_dimensions,"x")
            @info "Assuming velocities normal to x boundaries are zero"
        end
        if BCs.south == PeriodicBoundaryCondition()
            push!(periodic_dimensions,"y")
        elseif Immersed == false && !isnothing(BCs.south)   
            push!(bounded_dimensions,"y")
            @info "Assuming velocities normal to y boundaries are zero"
        end
        if BCs.top == PeriodicBoundaryCondition()
            push!(periodic_dimensions,"z")
        elseif Immersed == false && !isnothing(BCs.top)
            push!(bounded_dimensions,"z")
            @info "Assuming velocities normal to z boundaries are zero"
        end
        
        # First add the necessary serialized entry for each new variable
        for var in var_names_to_filter 
            new_path = "timeseries/$var"*"_Lagrangian_filtered_at_mean/serialized"
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
                        Xi_u = coord_dict["x_mesh"] .+ file["timeseries/xi_u/$iter"]

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
                        Xi_v =  coord_dict["y_mesh"] .+ file["timeseries/xi_v/$iter"]

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
                        Xi_w = coord_dict["z_mesh"] .+ file["timeseries/xi_w/$iter"]

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
            for var in var_names_to_filter
                var_data = file["timeseries/$var"*"_Lagrangian_filtered/$iter"]

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
            
            for (ivar,var) in enumerate(var_names_to_filter)
                
                # This is the main interpolation
                values = var_data[:,ivar]
                interpolator = scipy_interpolate.LinearNDInterpolator(coords, values)
                interp_data = pyconvert(Array,interpolator(regular_coord_mesh...))

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
                        
                        if bounded_dim == "x"
                            if n_true_dims == 2
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
                                y_mesh = coord_dict["y_mesh"][edge_index_with_halos,:,:]
                                z_mesh = coord_dict["z_mesh"][edge_index_with_halos,:,:]
                                meshes = (y_mesh,z_mesh)
                                interpolator = scipy_interpolate.LinearNDInterpolator(coords_cut, values_cut)
                                interp_data[edge_index_with_halos,:,:] = pyconvert(Array,interpolator(meshes...))
                            else
                                @info "Number of dimensions $n_true_dims is not implemented"
                            end
            
                        elseif bounded_dim =="y"
                            if n_true_dims == 2
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
                                x_mesh = coord_dict["x_mesh"][:,edge_index_with_halos,:]
                                z_mesh = coord_dict["z_mesh"][:,edge_index_with_halos,:]
                                meshes = (x_mesh,z_mesh)
                                interpolator = scipy_interpolate.LinearNDInterpolator(coords_cut, values_cut)
                                interp_data[:,edge_index_with_halos,:] = pyconvert(Array,interpolator(meshes...))
                            else
                                @info "Number of dimensions $n_true_dims is not implemented"
                            end
                        elseif bounded_dim == "z"
                            if n_true_dims == 2
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
                                x_mesh = coord_dict["x_mesh"][:,:,edge_index_with_halos]
                                y_mesh = coord_dict["y_mesh"][:,:,edge_index_with_halos]
                                meshes = (x_mesh,y_mesh)
                                interpolator = scipy_interpolate.LinearNDInterpolator(coords_cut, values_cut)
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
                new_var_loc = "timeseries/$var"*"_Lagrangian_filtered_at_mean/$iter"
                if haskey(file, new_var_loc)
                    Base.delete!(file, new_var_loc) #incase we already tried to write this variable
                end
                file[new_var_loc] = interp_data
                
            end
        end
    end
    @info "Wrote regridded data to new variables with _at_mean suffix in file $output_filename"
    
end

"""
    jld2_to_netcdf(jld2_filename::String, nc_filename::String)

Converts a JLD2 output file generated by an Oceananigans simulation into a
standard NetCDF file. This function is useful for post-processing and for
sharing data with other tools that expect the NetCDF format.

The conversion process involves the following steps:
1.  **Read JLD2 data**: Opens the input JLD2 file and reads the grid, time,
    and all timeseries variables.
2.  **Create NetCDF file**: Creates a new NetCDF file with a `.nc` extension.
3.  **Define dimensions**: Defines NetCDF dimensions based on the grid sizes
    and staggered locations (e.g., `x_caa` for cell centers, `x_faa` for
    cell faces).
4.  **Define grid variables**: Writes the grid coordinates and metadata
    (e.g., `Lx`, `Ny`, `Hx`) as variables to the NetCDF file.
5.  **Write timeseries data**: Iterates through each variable in the JLD2
    file's timeseries, determines its location on the grid, and writes the
    data to a new variable in the NetCDF file.
6.  **Add metadata**: Adds attributes to each variable, including boundary
    conditions and units, for better documentation.

Arguments
=========
- `jld2_filename`: A `String` specifying the path to the input JLD2 file.
- `nc_filename`: A `String` specifying the path for the output NetCDF file.
"""
function jld2_to_netcdf(jld2_filename::String, nc_filename::String)
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

        # There might be a time shifted variable
        if haskey(file["timeseries"], "t_shifted")
            t_shifted = [file["timeseries/t_shifted/$iter"] for iter in iterations]
            nctime = defVar(ds,"time_shifted", Float64, ("time",), attrib = OrderedDict(
            "units"                     => "seconds",
            "long_name"                 => "Time shifted to mean time",
            ))
            nctime[:] = t_shifted
        end

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
            if varname ∉ ("t","t_shifted") 

                # Create a variable in the NetCDF file
                bc_string = sprint(show, file["timeseries/$varname/serialized/boundary_conditions"])
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
    @info "Wrote NetCDF file to $nc_filename"
end

"""
    get_weight_function(;t::AbstractArray, tref::Real, filter_params::NamedTuple, direction::String = "both")

Computes the weighting function for the offline filter. This function calculates
the filter's impulse response, which determines how much each point in the
timeseries `t` contributes to the filtered value at a reference time `tref`.
The weighting function is based on the provided `filter_params`, which contains
the coefficients for the filter's impulse response.

Keyword arguments
=========
- `t`: A collection of time points in the timeseries.
- `tref`: The reference time at which the filter is being evaluated.
- `filter_params`: A `NamedTuple` containing the coefficients (`a`, `b`, `c`,
  `d`) and the number of coefficient pairs (`N_coeffs`).
- `direction`: A `String` indicating the direction of the filter. It can be
  "both" (default), "forward", or "backward". This determines whether the
  filter is applied symmetrically around `tref`, only to past times, or only
  to future times.

Returns
=======
- A vector of weights `G`, with the same dimensions as `t`, representing the
  value of the filter's impulse response at each time point relative to `tref`.
"""
function get_weight_function(;t::AbstractArray, tref::Real, filter_params::NamedTuple, direction::String = "both")

    G = 0*t
    N_coeffs = filter_params.N_coeffs
    for i in 1:N_coeffs
        
        a = getproperty(filter_params, Symbol("a$i"))
        b = getproperty(filter_params, Symbol("b$i"))
        c = getproperty(filter_params, Symbol("c$i"))
        d = getproperty(filter_params, Symbol("d$i"))

        G += (a.*cos.(d.*abs.(t .- tref)) .+ b.*sin.(d.*abs.(t .- tref))).*exp.(-c.*abs.(t .- tref))
    end
    if direction == "forward"
        G[t .> tref] .= 0
    elseif direction == "backward"
        G[t .< tref] .= 0
    elseif direction != "both"
        error("Direction must be 'forward', 'backward' or 'both'")
    end
    return G
end

#TODO do this for forward and backward only too
"""
    get_frequency_response(freq::AbstractArray, filter_params::NamedTuple)

Calculates the frequency response of the offline filter. This function takes a set
of frequencies and the filter's coefficients to compute how the filter amplifies
or attenuates different frequency components of a signal.

The response is computed by summing the contributions of each coefficient pair
based on the filter's transfer function in the frequency domain. The result is
a measure of the filter's gain at each given frequency.

Arguments
=========
- `freq`: A vector of frequencies (in radians per unit time).
- `filter_params`: A `NamedTuple` containing the filter coefficients (`a`, `b`,
  `c`, `d`) and the number of coefficient pairs (`N_coeffs`).

Returns
=======
- A vector `Ghat` representing the filter's frequency response at each
  corresponding frequency in `freq`.
"""
function get_frequency_response(freq::AbstractArray, filter_params::NamedTuple)
    
    Ghat = 0*freq
    N_coeffs = filter_params.N_coeffs
    for i in 1:N_coeffs
        
        a = getproperty(filter_params, Symbol("a$i"))
        b = getproperty(filter_params, Symbol("b$i"))
        c = getproperty(filter_params, Symbol("c$i"))
        d = getproperty(filter_params, Symbol("d$i"))

        Ghat += (a*c .+ b.*(d .+ freq))./(c^2 .+ (d .+ freq).^2) .+ (a*c .+ b.*(d .- freq))./(c^2 .+ (d .- freq).^2)
    end

    return Ghat
end

"""
    compute_Eulerian_filter!(config::AbstractConfig)

Computes the Eulerian filter for specified variables and writes the results to a
combined output file. This function performs a direct, convolution-style
filtering of a time series by applying a weighting function to the data at each
time step.

The function iterates through each variable to be filtered:
1.  **Reads data**: The entire time series of the variable is read from the
    JLD2 file.
2.  **Applies weighting**: At each output time, a weighting function `G` is
    computed and applied to the entire time series. The weighted data is summed
    to produce the filtered field.
3.  **Writes output**: The resulting filtered field is saved back to the
    same JLD2 file in a new group with the `_Eulerian_filtered` suffix.

This method serves as a benchmark for comparison with the main Lagrangian filter.

Arguments
=========
- `config`: An instance of `AbstractConfig` containing the file path, variable
  names, and filter parameters.
"""
function compute_Eulerian_filter!(config::AbstractConfig)
    filter_params = config.filter_params
    output_filename = config.output_filename
    var_names_to_filter = config.var_names_to_filter
    filter_mode = config.filter_mode
    if filter_mode == "offline"
        direction = "both"
    elseif filter_mode == "online"
        direction = "forward"
    end

    # Open existing file
    jldopen(output_filename,"r+") do file
        iterations = parse.(Int, keys(file["timeseries/t"]))
        times = [file["timeseries/t/$iter"] for iter in iterations]
        dt = times[2] - times[1] # Initialise time interval

        # Loop over variables to filter
        for var_name in var_names_to_filter
            @info "Computing Eulerian filter for variable $var_name"

            # Create group for filtered data
            g_EF = Group(file, "timeseries/$(var_name * "_Eulerian_filtered")")

            # Copy over serialized properties
            g_EF_serialized = Group(file, "timeseries/$(var_name * "_Eulerian_filtered")/serialized") 
            for property in keys(file["timeseries/$var_name/serialized"])
                g_EF_serialized[property] = file["timeseries/$var_name/serialized/$property"]
            end

            # Loop over times to compute filtered field at each time
            for (i, t) in enumerate(times)
                @info "Computing Eulerian filter at time $t =  (index $i of $(length(times)))"
                G = get_weight_function(t = times, tref = t, filter_params = filter_params, direction = direction)
                
                # Initialise with zeros
                mean_field = file["timeseries/$(var_name)/$(iterations[1])"]*0.0

                normalisation = sum(G[2:end] .* diff(times))
                # Construct mean sequentially
                for j in 1:length(times)
                    field = file["timeseries/$(var_name)/$(iterations[j])"]
                    dt = j < length(times) ? times[j+1] - times[j] : dt
                    mean_field .+= G[j] .* field .* dt 
                end
                g_EF["$(iterations[i])"] = mean_field./normalisation
            end
        end
    end
end

"""
    compute_time_shift!(config::AbstractConfig)

Computes the time shift for an online filter based on its coefficients and writes
the shifted time series to the output file.

This function is only applicable for online filtering. The time shift is computed
as the time delay introduced by the filter's transfer function. This new time series
is stored in a new group called `timeseries/t_shifted` within the output JLD2 file.

# Arguments
- `config`: A configuration object of type `OfflineFilterConfig` which contains
  the `filter_mode`, `output_filename`, and `filter_params` (filter coefficients).

# Throws
- `error`: If `config.filter_mode` is not `"online"`. The time shift for offline
  (forward-backward) filtering is zero by definition due to an even weight function.

"""
function compute_time_shift!(config::AbstractConfig)
    if config.filter_mode != "online"
        error("Time shift computation is only relevant for online filtering. Offline forward-backward filtering
        has an even weight function, so time shift is zero .")
    end
    output_filename = config.output_filename
    filter_params = config.filter_params    
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