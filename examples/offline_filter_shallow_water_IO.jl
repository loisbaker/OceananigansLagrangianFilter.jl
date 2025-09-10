using OceananigansLagrangianFilter
using Oceananigans.Units
using CUDA # If GPU available
using CairoMakie

# Set up the filter configuration
filter_config = OfflineFilterConfig(original_data_filename="SW_IO_with_tracer.jld2", # Where the original simulation output is
                                    output_filename = "SW_IO_with_tracer_filtered.jld2", # Where to save the filtered output
                                    var_names_to_filter = ("T",), # Which variables to filter
                                    velocity_names = ("u","v"), # Velocities to use for Lagrangian filtering
                                    architecture = GPU(), # CPU() or GPU() if available
                                    Δt = 1e-2, # Time step of filtering simulation
                                    T_out=0.1, # How often to output filtered data
                                    N=4, # Order of Butterworth filter
                                    freq_c = 0.5, # Cut-off frequency of Butterworth filter
                                    output_netcdf = true, # Whether to output filtered data to a netcdf file in addition to .jld2
                                    delete_intermediate_files = true, # Delete the individual output of the forward and backward passes
                                    compute_Eulerian_filter = true) # Whether to compute the Eulerian filter for comparison


# Run the offline Lagrangian filter
run_offline_Lagrangian_filter(filter_config)

# Animate the results
timeseries1 = FieldTimeSeries(filter_config.output_filename, "T")
timeseries2 = FieldTimeSeries(filter_config.output_filename, "T_Eulerian_filtered")
timeseries3 = FieldTimeSeries(filter_config.output_filename, "T_Lagrangian_filtered")
timeseries4 = FieldTimeSeries(filter_config.output_filename, "T_Lagrangian_filtered_at_mean")

times = timeseries1.times

set_theme!(Theme(fontsize = 20))
fig = Figure(size = (1300, 500))

axis_kwargs = (xlabel = "x",
               ylabel = "y",
               limits = ((0, 2*pi), (0, 2*pi)),
               aspect = AxisAspect(1))

ax1 = Axis(fig[2, 1]; title = "Raw", axis_kwargs...)
ax2 = Axis(fig[2, 2]; title = "Eulerian filtered", axis_kwargs...)
ax3 = Axis(fig[2, 3]; title = "Lagrangian filtered", axis_kwargs...)
ax4 = Axis(fig[2, 4]; title = "Lagrangian filtered at mean", axis_kwargs...)


n = Observable(1)
Observable(1)

var1 = @lift timeseries1[$n]
var2 = @lift timeseries2[$n]
var3 = @lift timeseries3[$n]
var4 = @lift timeseries4[$n]

heatmap!(ax1, var1; colormap = :Spectral, colorrange = (0, 1))
heatmap!(ax2, var2; colormap = :Spectral, colorrange = (0, 1))
heatmap!(ax3, var3; colormap = :Spectral, colorrange = (0, 1))
heatmap!(ax4, var4; colormap = :Spectral, colorrange = (0, 1))


title = @lift "Tracer T at time = " * string(round(times[$n], digits=2))
Label(fig[1, 1:4], title, fontsize=24, tellwidth=false)

fig

frames = 1:length(times)

@info "Making an animation"

CairoMakie.record(fig, "IO_filtered_tracer_movie_offline.mp4", frames, framerate=24) do i
    n[] = i
end