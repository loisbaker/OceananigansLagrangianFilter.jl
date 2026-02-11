using CairoMakie, Colors
using MathTeXEngine
using OceananigansLagrangianFilter.Utils: set_offline_BW2_filter_params, set_online_BW_filter_params, get_weight_function, get_offline_frequency_response, get_online_frequency_response

# Compute time shift for online filter
function compute_time_shift(filter_params::NamedTuple)
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
    return time_shift
end

function compute_offline_threshold_time(;freq_c::Real, threshold::Real = 0.95, N::Int = 1, Tmax::Real = 10, Nt::Int = 100)
    
    times = Array(LinRange(0, Tmax, Nt))
    dt = times[2] - times[1]
    filter_params = set_offline_BW2_filter_params(N = N, freq_c = freq_c)
    G = get_weight_function(t = times, tref = 0, filter_params = filter_params, direction = "both")

    half_int = dt * cumsum( (G[1:end-1] .+ G[2:end]) / 2 )
    half_int = [0.0; half_int] # align with times
    full_int = 0.5 .+ half_int
    reverse_int = 1 .- full_int
    # Want earliest time after which integral never dips below threshold
    idx = findfirst(>=(1 - threshold), abs.(reverse_int[end:-1:1]))
    return isnothing(idx) ? NaN : times[end:-1:1][idx]
end

function compute_online_threshold_time(;freq_c::Real, threshold::Real = 0.95, N::Int = 1, Tmax::Real = 10, Nt::Int = 100)
    
    times = Array(LinRange(0, Tmax, Nt))
    dt = times[2] - times[1]
    filter_params = set_online_BW_filter_params(N = N, freq_c = freq_c)
    G = get_weight_function(t = times, tref = 0, filter_params = filter_params, direction = "backward")

    int = dt * cumsum( (G[1:end-1] .+ G[2:end]) / 2 )
    int = [0.0; int] # align with times
    reverse_int = 1 .- int
    # Want earliest time after which integral never dips below threshold

    idx = findfirst(>=(1 - threshold), abs.(reverse_int[end:-1:1]))
    return isnothing(idx) ? NaN : times[end:-1:1][idx]
end

Tmax = 15
freq_max = 4

# time and frequency arrays
t = Array(LinRange(-Tmax, Tmax, 1000))
freqs = Array(LinRange(-freq_max,  freq_max, 1000))
freq_c = 1

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
        labelsize = 30, # Increased from 24
        patchsize = (40, 20)
    ),
    Colorbar = (
       labelsize = 32
    )
))

fig = Figure(size = (2000, 1000))

c1 = RGB(0/255, 114/255, 178/255) # blue
c2 = RGB(0/255, 158/255, 115/255) # teal 
c3 = RGB(204/255, 121/255, 167/255) # purple
c4 = :black

# Set up axes
title1 = L"\text{Impulse response (weight function) } G(t)"
title2 = L"\text{Frequency response } \hat{G}(\omega)"

ax1 = Axis(fig[1, 1]; title = title1, xlabelvisible = false, xticklabelsvisible = false, yticklabelsvisible = true, ylabel = L"\text{Online filter}", limits = ((-Tmax, Tmax),(-0.1, 1.1)))
ax2 = Axis(fig[1, 2]; title = title2, xlabelvisible = false, xticklabelsvisible = false, ylabelvisible = false, yticklabelsvisible = true, limits = ((-freq_max, freq_max), (-1.1, 1.1)))
ax3 = Axis(fig[2, 1]; xlabel = L"\text{Time}", ylabel = L"\text{Offline filter}", xticklabelsvisible = true, yticklabelsvisible = true, limits = ((-Tmax, Tmax),(-0.1, 1.1)))
ax4 = Axis(fig[2, 2]; ylabelvisible = false, yticklabelsvisible = true, xticklabelsvisible = true, xlabel = L"\text{Frequency}", limits = ((-freq_max, freq_max), (-1.1, 1.1)))

# ax1 (impulse response online)
lines!(ax1, t, get_weight_function(t = t, tref = 0, filter_params = set_online_BW_filter_params(N = 1, freq_c = freq_c), direction = "backward"), linewidth = 5, label = L"N=1", color = c1)
lines!(ax1, t, get_weight_function(t = t, tref = 0, filter_params = set_online_BW_filter_params(N = 2, freq_c = freq_c), direction = "backward"), linewidth = 5, label = L"N=2", color = c2)
lines!(ax1, t, get_weight_function(t = t, tref = 0, filter_params = set_online_BW_filter_params(N = 4, freq_c = freq_c), direction = "backward"), linewidth = 5, label = L"N=4", color = c3)
vlines!(ax1, 2*pi/freq_c; color = c4, linewidth = 4, linestyle = :dash, label = L"\text{Cut-off period}")
vlines!(ax1, -2*pi/freq_c; color = c4, linewidth = 4, linestyle = :dash)

t_95_N1 = compute_online_threshold_time(freq_c = freq_c, threshold = 0.99, N = 1, Tmax = Tmax, Nt = 2000)
vlines!(ax1, t_95_N1; color = c1, linewidth = 4, linestyle = :dashdot, label = L"99\% \text{ threshold } N=1")

t_95_N2 = compute_online_threshold_time(freq_c = freq_c, threshold = 0.99, N = 2, Tmax = Tmax, Nt = 2000)
vlines!(ax1, t_95_N2; color = c2, linewidth = 4, linestyle = :dashdot, label = L"99\% \text{ threshold } N=2")

t_95_N4 = compute_online_threshold_time(freq_c = freq_c, threshold = 0.99, N = 4, Tmax = Tmax, Nt = 2000)
vlines!(ax1, t_95_N4; color = c3, linewidth = 4, linestyle = :dashdot, label = L"99\% \text{ threshold } N=4")

# Scatter points for time shifts
tshift_N1 = compute_time_shift(set_online_BW_filter_params(N = 1, freq_c = freq_c))
scatter!(ax1, tshift_N1, get_weight_function(t = Array([tshift_N1]), tref = 0, filter_params = set_online_BW_filter_params(N = 1, freq_c = freq_c), direction = "backward")[1], markersize = 20, color = c1)

tshift_N2 = compute_time_shift(set_online_BW_filter_params(N = 2, freq_c = freq_c))
scatter!(ax1, tshift_N2, get_weight_function(t = Array([tshift_N2]), tref = 0, filter_params = set_online_BW_filter_params(N = 2, freq_c = freq_c), direction = "backward")[1], markersize = 20, color = c2)

tshift_N4 = compute_time_shift(set_online_BW_filter_params(N = 4, freq_c = freq_c))
scatter!(ax1, tshift_N4, get_weight_function(t = Array([tshift_N4]), tref = 0, filter_params = set_online_BW_filter_params(N = 4, freq_c = freq_c), direction = "backward")[1], markersize = 20, color = c3)

# ax2 (frequency response online) ---
lines!(ax2, freqs, abs.(get_online_frequency_response(freq = freqs, filter_params = set_online_BW_filter_params(N = 1, freq_c = freq_c))), linewidth = 5, label = L"N=1", color = c1)
lines!(ax2, freqs, real.(get_online_frequency_response(freq = freqs, filter_params = set_online_BW_filter_params(N = 1, freq_c = freq_c))), linewidth = 5, color = c1, linestyle = :dash)
lines!(ax2, freqs, imag.(get_online_frequency_response(freq = freqs, filter_params = set_online_BW_filter_params(N = 1, freq_c = freq_c))), linewidth = 5, color = c1, linestyle = :dot)

lines!(ax2, freqs, abs.(get_online_frequency_response(freq = freqs, filter_params = set_online_BW_filter_params(N = 2, freq_c = freq_c))), linewidth = 5, label = L"N=2", color = c2)
lines!(ax2, freqs, real.(get_online_frequency_response(freq = freqs, filter_params = set_online_BW_filter_params(N = 2, freq_c = freq_c))), linewidth = 5, color = c2, linestyle = :dash)
lines!(ax2, freqs, imag.(get_online_frequency_response(freq = freqs, filter_params = set_online_BW_filter_params(N = 2, freq_c = freq_c))), linewidth = 5, color = c2, linestyle = :dot)

lines!(ax2, freqs, abs.(get_online_frequency_response(freq = freqs, filter_params = set_online_BW_filter_params(N = 4, freq_c = freq_c))), linewidth = 5, label = L"N=4", color = c3)
lines!(ax2, freqs, real.(get_online_frequency_response(freq = freqs, filter_params = set_online_BW_filter_params(N = 4, freq_c = freq_c))), linewidth = 5, color = c3, linestyle = :dash)
lines!(ax2, freqs, imag.(get_online_frequency_response(freq = freqs, filter_params = set_online_BW_filter_params(N = 4, freq_c = freq_c))), linewidth = 5, color = c3, linestyle = :dot)

vlines!(ax2, freq_c; color = c4, linewidth = 4, linestyle = :dash, label = L"\text{Cut-off frequency}")
vlines!(ax2, -freq_c; color = c4, linewidth = 4, linestyle = :dash)

# ax3 (impulse response offline)
lines!(ax3, t, get_weight_function(t = t, tref = 0, filter_params = set_offline_BW2_filter_params(N = 1, freq_c = freq_c), direction = "both"), linewidth = 5, label = L"N=1", color = c1)
lines!(ax3, t, get_weight_function(t = t, tref = 0, filter_params = set_offline_BW2_filter_params(N = 2, freq_c = freq_c), direction = "both"), linewidth = 5, label = L"N=2", color = c2)
lines!(ax3, t, get_weight_function(t = t, tref = 0, filter_params = set_offline_BW2_filter_params(N = 4, freq_c = freq_c), direction = "both"), linewidth = 5, label = L"N=4", color = c3)
vlines!(ax3, 2*pi/freq_c; color = c4, linewidth = 4, linestyle = :dash, label = L"\text{Cut-off period}") 
vlines!(ax3, -2*pi/freq_c; color = c4, linewidth = 4, linestyle = :dash) 

t_95_N1_off = compute_offline_threshold_time(freq_c = freq_c, threshold = 0.99, N = 1, Tmax = Tmax, Nt = 2000)
vlines!(ax3, t_95_N1_off; color = c1, linewidth = 4, linestyle = :dashdot, label = L"99\% \text{ threshold } N=1")

t_95_N2_off = compute_offline_threshold_time(freq_c = freq_c, threshold = 0.99, N = 2, Tmax = Tmax, Nt = 2000)
vlines!(ax3, t_95_N2_off; color = c2, linewidth = 4, linestyle = :dashdot, label = L"99\% \text{ threshold } N=2")

t_95_N4_off = compute_offline_threshold_time(freq_c = freq_c, threshold = 0.99, N = 4, Tmax = Tmax, Nt = 2000)
vlines!(ax3, t_95_N4_off; color = c3, linewidth = 4, linestyle = :dashdot, label = L"99\% \text{ threshold } N=4")

# ax4 (frequency response offline)
lines!(ax4, freqs, get_offline_frequency_response(freq = freqs, filter_params = set_offline_BW2_filter_params(N = 1, freq_c = freq_c)), linewidth = 5, label = L"N=1", color = c1)
lines!(ax4, freqs, get_offline_frequency_response(freq = freqs, filter_params = set_offline_BW2_filter_params(N = 2, freq_c = freq_c)), linewidth = 5, label = L"N=2", color = c2)
lines!(ax4, freqs, get_offline_frequency_response(freq = freqs, filter_params = set_offline_BW2_filter_params(N = 4, freq_c = freq_c)), linewidth = 5, label = L"N=4", color = c3)
vlines!(ax4, freq_c; color = c4, linewidth = 4, linestyle = :dash, label = L"\text{Cut-off frequency}")
vlines!(ax4, -freq_c; color = c4, linewidth = 4, linestyle = :dash)

# Final ticks and legends 
ωc = freq_c
ax3.xticks = ([-2π/ωc, 0.0, 2π/ωc], [L"-2\pi/\omega_c", L"0", L"2\pi/\omega_c"])
ax4.xticks = ([-ωc, 0.0, ωc], [L"-\omega_c", L"0", L"\omega_c"])

ax1.yticks = ([0, ωc/2, ωc], [L"0", L"\omega_c/ 2", L"\omega_c"])
ax3.yticks = ([0, ωc/2, ωc], [L"0", L"\omega_c/ 2", L"\omega_c"])

ax2.yticks = ([-1.0, 0.0, 1.0], [L"-1.0", L"0.0", L"1.0"])
ax4.yticks = ([-1.0, 0.0, 1.0], [L"-1.0", L"0.0", L"1.0"])

axislegend(ax1, position = :lb)
axislegend(ax2, position = :lb)

# Custom legend for common line styles in ax2
elem_abs = [LineElement(color = :grey50, linewidth = 5, linestyle = nothing)]
elem_real = [LineElement(color = :grey50, linewidth = 5, linestyle = :dash)]
elem_imag = [LineElement(color = :grey50, linewidth = 5, linestyle = :dot)]

Legend(fig[1, 2], 
    [elem_abs, elem_real, elem_imag], 
    [L"\text{Absolute value}", L"\text{Real part}", L"\text{Imaginary part}"],
    labelsize = 30, # Explicitly setting size for consistency
    tellwidth = false, tellheight = false,
    halign = :right, valign = :bottom,
    margin = (10, 10, 10, 10)
)

# Add simple text labels (a, b, c, d)
# We place them in the TopLeft() of the grid position with some padding
text!(ax1, 0.02, 0.98, text = "a", space = :relative, fontsize = 40, font = :bold, align = (:left, :top))
text!(ax2, 0.02, 0.98, text = "b", space = :relative, fontsize = 40, font = :bold, align = (:left, :top))
text!(ax3, 0.02, 0.98, text = "c", space = :relative, fontsize = 40, font = :bold, align = (:left, :top))
text!(ax4, 0.02, 0.98, text = "d", space = :relative, fontsize = 40, font = :bold, align = (:left, :top))

save("filter_responses.png", fig,px_per_unit = 4)
fig