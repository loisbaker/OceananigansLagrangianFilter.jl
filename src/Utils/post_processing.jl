# Post-processing of the filter output: combining the forward and backward passes, regridding to the mean
# position, and the Eulerian filter and time shift for comparison

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
