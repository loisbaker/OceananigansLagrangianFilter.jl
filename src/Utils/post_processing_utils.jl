using Oceananigans.BoundaryConditions: PeriodicBoundaryCondition
using DataStructures: OrderedDict
using NCDatasets
using JLD2
using JLD2: Group

# Python import to use the LinearNDInterpolator for regridding
using PythonCall
const scipy_interpolate = PythonCall.pynew()
function __init__()
    PythonCall.pycopy!(scipy_interpolate, pyimport("scipy.interpolate"))
end

"""
    sum_forward_backward_contributions!(combined_output_filename::String, 
                                       forward_output_filename::String, 
                                       backward_output_filename::String, 
                                       T::Real, 
                                       velocity_names::Tuple{Vararg{String}}, 
                                       var_names_to_filter::Tuple{Vararg{String}})

Combines the results of a forward-in-time and a backward-in-time filter simulation by summing
their contributions at the correct times and saving in a new file.


# Arguments
- `combined_output_filename::String`: The path to the output file for the combined data.
- `forward_output_filename::String`: The path to the JLD2 file containing the forward simulation output.
- `backward_output_filename::String`: The path to the JLD2 file containing the backward simulation output.
- `T::Real`: The total duration of the filter simulation (end time).
- `velocity_names::Tuple{Vararg{String}}`: A tuple of names for velocity variables to be filtered (e.g., ("u", "v")).
- `var_names_to_filter::Tuple{Vararg{String}}`: A tuple of names for other variables to be filtered (e.g., ("T", "S")).
"""

function sum_forward_backward_contributions!(config)
# Combine the forward and backward simulations by summing them into a single file

combined_output_filename = config.combined_output_filename
forward_output_filename = config.forward_output_filename
backward_output_filename = config.backward_output_filename
T = config.T
velocity_names = config.velocity_names
var_names_to_filter = config.var_names_to_filter

    # List the names of the fields that we will combine
    filtered_var_names = vcat(["xi_" * vel for vel in velocity_names], [var * "_filtered" for var in var_names_to_filter])
    
    jldopen(combined_output_filename,"w") do combined_file
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
    println("Combined forward and backward contributions into $combined_output_filename")

end

function _remove_halos(data,grid)
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

function regrid_to_mean_position!(config)
    combined_output_filename = config.combined_output_filename
    var_names_to_filter = config.var_names_to_filter
    velocity_names = config.velocity_names
    npad = config.npad 

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
        for iter in iterations
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
                        Xi_u = _remove_halos(Xi_u,grid)
                        # move xi points outside of domain into domain
                        if "x" in periodic_dimensions                           
                            Xi_u .-= floor.((Xi_u .- grid.xᶠᵃᵃ[1])./grid.Lx) .* grid.Lx
                        end
                        push!(Xi_list,vec(Xi_u))


                    elseif (dim == "x") && !("u" in velocity_names)
                        Xi_u = coord_dict["x_mesh"]
                        # Lose the halo regions 
                        Xi_u = _remove_halos(Xi_u,grid)
                        push!(Xi_list,vec(Xi_u))

                    elseif (dim == "y") && ("v" in velocity_names)
                        Xi_v =  coord_dict["y_mesh"] .+ file["timeseries/xi_v/$iter"]
                        # Lose the halo regions 
                        Xi_v = _remove_halos(Xi_v,grid)
                        # move xi points outside of domain + halo regions into domain
                        if "y" in periodic_dimensions
                            Xi_v .-= floor.((Xi_v .- grid.yᵃᶠᵃ[1])./grid.Ly) .* grid.Ly
                        end
                        push!(Xi_list,vec(Xi_v))

                    elseif (dim == "y") && !("v" in velocity_names)
                        Xi_v = coord_dict["y_mesh"]
                        # Lose the halo regions 
                        Xi_v = _remove_halos(Xi_v,grid)
                        push!(Xi_list,vec(Xi_v))

                    elseif (dim == "z") && ("w" in velocity_names)
                        Xi_w = coord_dict["z_mesh"] .+ file["timeseries/xi_w/$iter"]
                        # Lose the halo regions 
                        Xi_w = _remove_halos(Xi_w,grid)
                        # move xi points outside of domain + halo regions into domain
                        if "z" in periodic_dimensions
                            Xi_w .-= floor.((Xi_w .- grid.z.cᵃᵃᶠ[1])./grid.Lz) .* grid.Lz
                        end
                        push!(Xi_list,vec(Xi_w))

                    elseif (dim == "z") && !("w" in velocity_names)
                        Xi_w = coord_dict["z_mesh"] .+ zeros(coord_dict["original_size"])
                        # Lose the halo regions 
                        Xi_w = _remove_halos(Xi_w,grid)
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
    println("Wrote regridded data to new variables with _regrid suffix in file $combined_output_filename")
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
            if varname ∉ ("t","t_simulation") 
                    
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
    println("Wrote NetCDF file to $nc_filename")
end

function get_weight_function(t, tref, filter_params)
    
    G = 0*t
    N_coeffs = filter_params.N_coeffs
    for i in 1:N_coeffs
        
        a = getproperty(filter_params, Symbol("a$i"))
        b = getproperty(filter_params, Symbol("b$i"))
        c = getproperty(filter_params, Symbol("c$i"))
        d = getproperty(filter_params, Symbol("d$i"))

        G += (a.*cos.(d.*abs.(t .- tref)) .+ b.*sin.(d.*abs.(t .- tref))).*exp.(-c.*abs.(t .- tref))
    end

    return G
end

function get_frequency_response(freq, filter_params)
    
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