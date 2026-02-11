using CairoMakie
using Printf
using Oceananigans 
using Oceananigans.Units
using MathTeXEngine

# Load data
output_filename = "data/lee_wave_offline_filtered.jld2"
it = 121 # 5 days (halfway through filtering period)

timeseries1 = FieldTimeSeries(output_filename, "u")
timeseries2 = FieldTimeSeries(output_filename, "u_Eulerian_filtered")
timeseries3 = FieldTimeSeries(output_filename, "u_Lagrangian_filtered")
timeseries4 = FieldTimeSeries(output_filename, "u_Lagrangian_filtered_at_mean")

b_timeseries1 = FieldTimeSeries(output_filename, "b")
b_timeseries2 = FieldTimeSeries(output_filename, "b_Eulerian_filtered")
b_timeseries3 = FieldTimeSeries(output_filename, "b_Lagrangian_filtered")
b_timeseries4 = FieldTimeSeries(output_filename, "b_Lagrangian_filtered_at_mean")

times = timeseries1.times
grid = timeseries1.grid
bottom_height = vec(grid.immersed_boundary.bottom_height)
Nx = grid.underlying_grid.Nx
Nz = grid.underlying_grid.Nz

# Set up coords
x_km = Array(grid.underlying_grid.xᶜᵃᵃ[1:Nx]) ./ 1000
z = Array(grid.underlying_grid.z.cᵃᵃᶜ[1:Nz])

# Set up theme
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
        labelsize = 30,
        patchsize = (40, 20)
    ),
    Colorbar = (
        labelsize = 32,
        ticksize = 10
    )
))

fig = Figure(size = (2000, 980))

# Set up axes
axis_kwargs = (xlabel = L"x \text{ [km]}",
               ylabel = L"z \text{ [m]}",
               limits = ((-20, 10), (-2000, 0)),
               aspect = AxisAspect(2.2))

ax1 = Axis(fig[1, 1]; title = L"\text{Raw } u", xlabelvisible = false, xticklabelsvisible = false, axis_kwargs...)
ax2 = Axis(fig[1, 2]; title = L"\text{Eulerian filtered } \overline{u}", xlabelvisible = false, xticklabelsvisible = false, ylabelvisible = false, yticklabelsvisible = false, axis_kwargs...)
ax3 = Axis(fig[2, 1]; title = L"\text{Lagrangian filtered } u*", axis_kwargs...)
ax4 = Axis(fig[2, 2]; title = L"\text{Lagrangian filtered } \overline{u}^\mathrm{L}", ylabelvisible = false, yticklabelsvisible = false, axis_kwargs...)

# Observables
n = Observable(it)
var1 = @lift timeseries1[$n]
var2 = @lift timeseries2[$n]
var3 = @lift timeseries3[$n]
var4 = @lift timeseries4[$n]

b_var1 = @lift b_timeseries1[$n]
b_var2 = @lift b_timeseries2[$n]
b_var3 = @lift b_timeseries3[$n]
b_var4 = @lift b_timeseries4[$n]

# Heatmaps and contours 
eps_val = 1e-5

# Top row ranges (Raw/Eulerian)
hm_1 = heatmap!(ax1, x_km, z, var1; colormap = :balance, colorrange = (0, 0.4))
hm_2 = heatmap!(ax2, x_km, z, var2; colormap = :balance, colorrange = (0, 0.4))

# Bottom row ranges (Lagrangian - much narrower)
hm_3 = heatmap!(ax3, x_km, z, var3; colormap = :balance, colorrange = (0.18-eps_val, 0.22+eps_val))
hm_4 = heatmap!(ax4, x_km, z, var4; colormap = :balance, colorrange = (0.18-eps_val, 0.22+eps_val))

# Buoyancy contours
for (ax, b_v) in zip([ax1, ax2, ax3, ax4], [b_var1, b_var2, b_var3, b_var4])
    contour!(ax, x_km, z, b_v; levels = 40, color = :black, linewidth = 0.8, alpha = 0.5)
end

# Colorbars
cb1 = Colorbar(fig[1, 3], hm_1, label = L"u \text{ [m s}^{-1}\text{]}") 
cb2 = Colorbar(fig[2, 3], hm_3, label = L"u \text{ [m s}^{-1}\text{]}")

# Topography
floor_level = fill(-2000, length(x_km)) 
for ax in [ax1, ax2, ax3, ax4]
    band!(ax, x_km, floor_level, bottom_height, color=:black)
end

# Redo x-axis ticks to avoid overlap
ax3.xticks = ([-20, -15, -10, -5, 0, 5])
ax4.xticks = ([-20, -15, -10, -5, 0, 5])

text!(ax1, 0.01, 1, text = "a", space = :relative, fontsize = 40, font = :bold, align = (:left, :top))
text!(ax2, 0.01, 1, text = "b", space = :relative, fontsize = 40, font = :bold, align = (:left, :top))
text!(ax3, 0.01, 1, text = "c", space = :relative, fontsize = 40, font = :bold, align = (:left, :top))
text!(ax4, 0.01, 1, text = "d", space = :relative, fontsize = 40, font = :bold, align = (:left, :top))

save("offline_lee_wave_lagrangian_filtering.png", fig,px_per_unit = 4)
fig