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

# The integral of exp(-lambda u) f over an interval 0 <= u <= dt, where f varies linearly from f_near (at u = 0) to
# f_far (at u = dt), is near_weight * f_near + far_weight * f_far, with
#   far_weight = (1 - exp(-lambda dt) (1 + lambda dt)) / (lambda^2 dt)
#   near_weight = (1 - exp(-lambda dt)) / lambda - far_weight
# Where lambda dt is small, far_weight comes from its series, as the formula loses precision.
@inline function far_weight(lambda, dt)
    z = lambda * dt
    return abs(z) < 1e-3 ? dt * (1/2 - z/3 + z^2/8) : (-expm1(-z) - z * exp(-z)) / (lambda^2 * dt)
end
@inline near_weight(lambda, dt) = -expm1(-lambda * dt) / lambda - far_weight(lambda, dt)

# Moves the running sum of each term on by dt, to the next saved time: it decays by exp(-lambda dt) and gains the
# integral over the interval, with the field taken to vary linearly from f_far (the previous saved time) to f_near
function step_running_sums!(running_sums, terms, f_near, f_far, dt)
    for (running_sum, term) in zip(running_sums, terms)
        running_sum .= exp.(-term.lambda .* dt) .* running_sum .+
                       near_weight.(term.lambda, dt) .* f_near .+ far_weight.(term.lambda, dt) .* f_far
    end
    return running_sums
end

# The sum over the terms of Re[amplitude * running sum], where amplitude is the past or future amplitude of the term
running_sums_value(running_sums, terms, amplitude) =
    sum(real.(getproperty(term, amplitude) .* running_sum) for (running_sum, term) in zip(running_sums, terms))

# The Eulerian filter as running sums, so that each field is read at most twice (rather than once for every output
# time). read(j) gives the input f_j for each weight function in terms_list (each a list of exponential terms), and for
# each output time t_i, output(i, values) is called with the integrals of G(t_i - t) f(t) over t, one for each weight
# function. Between saved times each field is taken to vary linearly, and the exponentials are integrated exactly. The
# ends are treated as by the Lagrangian filter: the sums over earlier times start from their steady values f / lambda,
# as if each field were constant before the first saved time; and offline (direction = :both), the sums over later
# times start at the last saved time from the sums over earlier times, as the backward pass of the offline filter
# starts from the end of the forward pass. The sums are kept in `store` (anything indexed by a String, e.g. a JLD2
# file) until the outputs are written, in order.
function Eulerian_running_sums(read, output, times, terms_list, direction, store)
    n = length(times)

    # Sweep through the saved times in `order` from the running sums `start` (the steady values if `nothing`), calling
    # emit(i, values) at each; returns the running sums at the last time of the sweep
    function sweep(order, amplitude, start, emit)
        running_sums = deepcopy(start)
        previous_inputs, previous_i = nothing, 0
        for i in order
            inputs = read(i)
            if isnothing(running_sums)
                running_sums = [[complex.(f ./ term.lambda) for term in terms]
                                for (terms, f) in zip(terms_list, inputs)]
            elseif previous_i != 0
                dt = abs(times[i] - times[previous_i])
                for (terms, sums, f, f_previous) in zip(terms_list, running_sums, inputs, previous_inputs)
                    step_running_sums!(sums, terms, f, f_previous, dt)
                end
            end
            emit(i, [running_sums_value(sums, terms, amplitude) for (terms, sums) in zip(terms_list, running_sums)])
            previous_inputs, previous_i = inputs, i
        end
        return running_sums
    end

    if direction === :forward
        sweep(1:n, :past, nothing, output)
    else
        keep(key) = (i, values) -> foreach(((w, value),) -> store["$key/$i/$w"] = value, enumerate(values))
        sums_at_end = sweep(1:n, :past, nothing, keep("past"))
        sweep(n:-1:1, :future, sums_at_end, keep("future"))
        for i in 1:n
            output(i, [store["past/$i/$w"] .+ store["future/$i/$w"] for w in eachindex(terms_list)])
        end
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
the number of saved times. Offline, the sums are kept in a temporary file next to the output
file until they are written. The results are saved to the same JLD2 file in a new group
with the `_Eulerian_filtered` suffix.

This method serves as a benchmark for comparison with the main Lagrangian filter.

If the filter outputs its terms separately (`outputs = :separate`), each term is Eulerian filtered with its own
weight function, and named like the Lagrangian outputs with `_Eulerian_filtered` in the stem, e.g.
`b_Eulerian_filtered_C1_scaled` and `b_Eulerian_filtered_S1_scaled`.

The ends of the data are treated as by the Lagrangian filter, so that the two can be compared: each field is taken
to be constant before the first saved time (see [`initialise_filtered_vars_from_model`](@ref)), and offline, the
sums over later times start at the last saved time from the sums over earlier times, as the backward pass of the
offline filter starts from the end of the forward pass. Between saved times each field is taken to vary linearly, and
the weight function is integrated exactly, so a field that is constant in time is unchanged by a normalised filter.

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

    # The weight functions to filter with: the whole filter, or each of its terms if they are output separately
    weights = outputs === :combined ? [("", filter_params)] : separate_term_params(filter_params)

    # Open existing file
    jldopen(output_filename,"r+") do file
        iterations = parse.(Int, keys(file["timeseries/t"]))
        times = [file["timeseries/t/$iter"] for iter in iterations]

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

            # The terms of each weight function, scaled by the frequency scale where there is one
            M = saved_frequency_scale(config, file, var_name, size(file["timeseries/$(var_name)/$(iterations[1])"]))
            terms_list = [scale_terms(exponential_terms(p), M) for (_, p) in weights]
            read_inputs(i) = fill(file["timeseries/$(var_name)/$(iterations[i])"], length(weights))
            function write_outputs(i, values)
                for (g_EF, value) in zip(groups, values)
                    g_EF["$(iterations[i])"] = value
                end
            end

            # The offline filter keeps the running sums in a temporary file next to the output
            if direction === :both
                sums_filename = output_filename * ".Eulerian_sums.jld2"
                try
                    jldopen(sums_filename, "w") do store
                        Eulerian_running_sums(read_inputs, write_outputs, times, terms_list, direction, store)
                    end
                finally
                    rm(sums_filename, force = true)
                end
            else
                Eulerian_running_sums(read_inputs, write_outputs, times, terms_list, direction, nothing)
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
