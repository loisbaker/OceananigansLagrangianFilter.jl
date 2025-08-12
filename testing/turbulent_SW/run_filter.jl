using OceananigansLagrangianFilter
using Printf
using NCDatasets
using JLD2
using Oceananigans.Units: Time
using CairoMakie
using CUDA
using Oceananigans.TimeSteppers: reset!
using GeoStats: IDW

# User defined options

original_data_filename = joinpath(@__DIR__, "SW_vort.jld2")
T_start = 0
T_end = 10

arch = GPU()

# Set the output period
T_out = 0.1

# Set filter order and cut-off frequency
# Amplitude of frequency response of filter will be squared Butterworth order 2^N
N = 2 
freq_c = 2

# Define variables to filter
original_var_names = ("ω",)

# Define velocities to use for filtering
velocity_names = ("u","v")

# Set filtering parameters (this is for Butterworth-type, could define others here)
filter_params = set_BW_filter_params(N=N,freq_c=freq_c)

# Set the time step for the simulation
Δt = 1e-3

# Decide whether to solve for and output maps to generalised Lagrangian mean
map_to_mean = true

# Name output files
forward_output_filename = joinpath(@__DIR__, "forward_LF.jld2")
backward_output_filename = joinpath(@__DIR__, "backward_LF.jld2")
combined_output_filename = joinpath(@__DIR__, "combined_LF.jld2")

# Manipulate data on disk to have correct order 
T = set_data_on_disk!(original_data_filename, direction="forward", T_start = T_start, T_end = T_end)

