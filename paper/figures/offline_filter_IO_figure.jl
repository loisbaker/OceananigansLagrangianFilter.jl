using CairoMakie 
using Oceananigans
using Oceananigans.Units
using Printf
using Statistics
using MathTeXEngine

output_filename = "data/SW_IO_with_tracer_filtered.jld2"

# Visualisation
timeseries1 = FieldTimeSeries(output_filename, "T")
timeseries2 = FieldTimeSeries(output_filename, "T_Eulerian_filtered")
timeseries3 = FieldTimeSeries(output_filename, "T_Lagrangian_filtered")
timeseries4 = FieldTimeSeries(output_filename, "T_Lagrangian_filtered_at_mean")

times = timeseries1.times
Nx = timeseries1.grid.Nx
Ny = timeseries1.grid.Ny
x = timeseries1.grid.xᶜᵃᵃ[1:Nx]
y = timeseries1.grid.yᵃᶜᵃ[1:Ny]

Nt = length(times)
it = ceil(Int, Nt/2)

# Times from end when filter is 95% or 99% converged
t_99 = 15.17
t_95 = 7.79

# Extract Hovmuller slices
# Find index closest to y = 2.5
y_target = 2.5
j_target = argmin(abs.(y .- y_target))

# Extract Hovmuller data: T(x, y_target, t)
# We convert to Array to ensure plotting compatibility
hov1 = Array(timeseries1[1:Nx, j_target, 1, :])
hov2 = Array(timeseries2[1:Nx, j_target, 1, :])
hov3 = Array(timeseries3[1:Nx, j_target, 1, :])
hov4 = Array(timeseries4[1:Nx, j_target, 1, :])

# Pre-calculate the trajectories of the max tracer concentration
it_spinup = round(Int, t_95 / (times[2] - times[1]))
function get_max_trajectory(ts, x_grid, y_grid)
    traj = Point2f[]
    for t_idx in it_spinup:length(ts.times)-it_spinup
        data = Array(ts[t_idx])
        val, idx = findmax(data)
        push!(traj, Point2f(x_grid[idx[1]], y_grid[idx[2]]))
    end
    return traj
end

traj1 = get_max_trajectory(timeseries1, x, y)
traj2 = get_max_trajectory(timeseries2, x, y)
traj3 = get_max_trajectory(timeseries3, x, y)
traj4 = get_max_trajectory(timeseries4, x, y)

# convenience function : Point2f vector  →  (xs, ys)
xy(v::Vector{<:Point2f}) = (first.(v), last.(v))

xs1, ys1 = xy(traj1)
xs2, ys2 = xy(traj2)
xs3, ys3 = xy(traj3)
xs4, ys4 = xy(traj4)

centre_x = mean(xs3)
centre_y = mean(ys3)

range = 2

# Set theme to default 
set_texfont_family!(FontFamily("TeXGyreHeros"))

set_theme!(Theme(
    fontsize = 32,
    Axis = (
        xtickalign = 1, ytickalign = 1,
        xticksize = 10, yticksize = 10,
        xlabelsize = 32, ylabelsize = 32,
        titlesize = 34,
    ),
    Legend = (
        framevisible = true,
        backgroundcolor = :white,
        labelsize = 24,
        patchsize = (40, 20)
    ),
))

fig = Figure(size = (2000, 1300))

# Use L"..." for all axis labels
axis_top_kwargs = (xlabel = L"x",
               ylabel = L"y",
                limits = ((centre_x - range, centre_x + range), (centre_y - range, centre_y + range)),
               aspect = AxisAspect(1), alignmode = Inside())

axis_bottom_kwargs = (xlabel = L"x",
               ylabel = L"t",
               limits = ((centre_x - range, centre_x + range), (0, times[end])),
               alignmode = Inside())

# x-y plots at time index it
# Using L strings for titles to ensure they use the LaTeX engine
ax1 = Axis(fig[1, 1]; title = L"\text{Raw}", xlabelvisible = false, xticklabelsvisible = false, axis_top_kwargs...)
ax2 = Axis(fig[1, 2]; title = L"\text{Eulerian filtered}", xlabelvisible = false, xticklabelsvisible = false, ylabelvisible = false, yticklabelsvisible = false, axis_top_kwargs...)
ax3 = Axis(fig[1, 3]; title = L"\text{Lagrangian filtered}", xlabelvisible = false, xticklabelsvisible = false, ylabelvisible = false, yticklabelsvisible = false, axis_top_kwargs...)
ax4 = Axis(fig[1, 4]; title = L"\text{Lagrangian filtered (remapped)}", xlabelvisible = false, ylabelvisible = false, yticklabelsvisible = false, xticklabelsvisible = false, axis_top_kwargs...)

n = Observable(it)

var1 = @lift timeseries1[$n]
var2 = @lift timeseries2[$n]
var3 = @lift timeseries3[$n]
var4 = @lift timeseries4[$n]

hm1 = heatmap!(ax1, var1; colormap = :Spectral, colorrange = (0, 1))
heatmap!(ax2, var2; colormap = :Spectral, colorrange = (0, 1))
heatmap!(ax3, var3; colormap = :Spectral, colorrange = (0, 1))
heatmap!(ax4, var4; colormap = :Spectral, colorrange = (0, 1))

# Overlay max tracer trajectories
lines!(ax1, xs1, ys1 ; color = :black, linewidth = 2)
lines!(ax2, xs2, ys2 ; color = :black, linewidth = 2)
lines!(ax3, xs3, ys3 ; color = :black, linewidth = 2)
lines!(ax4, xs4, ys4 ; color = :black, linewidth = 2)

scatter!(ax1, [xs1[it-it_spinup]], [ys1[it-it_spinup]] ; color = :white, marker = :circle,
         markersize = 14, strokecolor = :black, strokewidth = 1.5)
scatter!(ax2, [xs2[it-it_spinup]], [ys2[it-it_spinup]] ; color = :white, marker = :circle,
         markersize = 14, strokecolor = :black, strokewidth = 1.5)
scatter!(ax3, [xs3[it-it_spinup]], [ys3[it-it_spinup]] ; color = :white, marker = :circle,
         markersize = 14, strokecolor = :black, strokewidth = 1.5)
scatter!(ax4, [xs4[it-it_spinup]], [ys4[it-it_spinup]] ; color = :white, marker = :circle,
         markersize = 14, strokecolor = :black, strokewidth = 1.5)

# Hovmuller plots
ax5 = Axis(fig[2, 1]; axis_bottom_kwargs...)
ax6 = Axis(fig[2, 2]; ylabelvisible = false, yticklabelsvisible = false, axis_bottom_kwargs...)
ax7 = Axis(fig[2, 3]; ylabelvisible = false, yticklabelsvisible = false, axis_bottom_kwargs...)
ax8 = Axis(fig[2, 4]; ylabelvisible = false, yticklabelsvisible = false, axis_bottom_kwargs...)

heatmap!(ax5, x, times, hov1; colormap = :Spectral, colorrange = (0, 1))
heatmap!(ax6, x, times, hov2; colormap = :Spectral, colorrange = (0, 1))
heatmap!(ax7, x, times, hov3; colormap = :Spectral, colorrange = (0, 1))
heatmap!(ax8, x, times, hov4; colormap = :Spectral, colorrange = (0, 1))

# Add horizontal lines to show 95 and 99% accuracy times

hlines!(ax6, [t_99, 40-t_99]; color = :black, linestyle = :dash, linewidth = 2)
hlines!(ax6, [t_95, 40-t_95]; color = :black, linestyle = :dot, linewidth = 2)

hlines!(ax7, [t_99, 40-t_99]; color = :black, linestyle = :dash, linewidth = 2)
hlines!(ax7, [t_95, 40-t_95]; color = :black, linestyle = :dot, linewidth = 2)

hlines!(ax8, [t_99, 40-t_99]; color = :black, linestyle = :dash, linewidth = 2)
hlines!(ax8, [t_95, 40-t_95]; color = :black, linestyle = :dot, linewidth = 2)

# Add horizontal lines on upper plot to show y = 2.5 (where the Hovmuller slices are taken)
hlines!(ax1, [2.5]; color = :black, linestyle = :dash, linewidth = 2)
hlines!(ax2, [2.5]; color = :black, linestyle = :dash, linewidth = 2)
hlines!(ax3, [2.5]; color = :black, linestyle = :dash, linewidth = 2)
hlines!(ax4, [2.5]; color = :black, linestyle = :dash, linewidth = 2)


Colorbar(fig[1:2, 5], hm1, label = L"\text{Tracer concentration, } T", valign = :bottom, height = Relative(0.935), labelsize = 32) 

rowgap!(fig.layout, -30)
fig

save("offline_IO_lagrangian_filtering.png", fig,px_per_unit = 4)