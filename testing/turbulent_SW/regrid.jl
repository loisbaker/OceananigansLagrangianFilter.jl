using OceananigansLagrangianFilter
using Printf
using JLD2
using Oceananigans.Units: Time
using CairoMakie
using CUDA
using Oceananigans.TimeSteppers: reset!
using PythonCall


# Define variables to filter
original_var_names = ("ω","T")

# Define velocities to use for filtering
velocity_names = ("u","v")

npad = 0
regrid_to_mean_position!(combined_output_filename, original_var_names, velocity_names, npad)
