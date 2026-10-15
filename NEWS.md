# OceananigansLagrangianFilter.jl release notes

## v0.2.0

### Breaking changes

- Oceananigans 0.113 and JLD2 0.6 are now required. Data saved by older Oceananigans versions can still be filtered, by passing the grid of the original simulation to the config.
- `map_to_mean` is replaced by `compute_maps` (output the maps) and `regrid_to_mean` (regrid to the mean position, which needs the maps). With Julie Meunier.
- The relaxation options are renamed: `relaxation_mask_func` and `relaxation_params` replace `boundary_relaxation`, `relax_timescale`, `mask_func` and `mask_params`.
- Offline, `advection = nothing` is replaced by `lagrangian = false` (see below).
- The `npad` and `backend` options are removed.
- Each pass of the offline filter now outputs its contribution to the filtered fields (`create_output_fields` has a `direction` keyword, and negates the mean velocities on the backward pass), and `sum_forward_backward_contributions!` sums them all, so it no longer takes `extra_filtered_velocity_names`.
- Input data is now always read by a `BufferedDataReader` from the original data file, so the intermediate input file and `create_input_data_on_disk` and `load_data` are removed. `initialise_filtered_vars_from_data` and `update_input_data!` now take a `BufferedDataReader` (made with `create_buffered_reader`).
- `get_frequency_response` replaces `get_offline_frequency_response` and `get_online_frequency_response`, and the `direction` of `get_weight_function` is now a `Symbol` (`:both`, `:forward` or `:backward`).
- Python is no longer needed: PythonCall, CondaPkg and SciPy are no longer dependencies.

### New features

- Regridding to the mean position now inverts the map with Newton's method, on the CPU or GPU, instead of using SciPy's Delaunay interpolation. It works in 3D, on stretched grids, with immersed boundaries and on a `LatitudeLongitudeGrid`, and returns `NaN` where no particle in the fluid has its mean position. `regrid_to_mean_position!` takes the `architecture` to run on and the maximum number of iterations `maxiter`. See the new theory page, "Regridding to the mean position".
- A frequency scale (`frequency_scale_func`, `frequency_scale_params` and `frequency_scale_dependencies`) lets the cutoff frequency vary in space and time. With Tatsu D Monkman.
- `lagrangian = false` computes the Eulerian filter with the same equations, online or offline. Online, `filtered_vars_advection` gives the model's advection schemes for it.
- Spectral filters, with `set_offline_spectrum_filter_params` and `set_online_spectrum_filter_params`. Filters can have sine terms that are odd in time (`sine_parity = :odd`), and their terms can be output separately (`outputs = :separate`). `get_weight_function` and `get_frequency_response` can give the response of a single term. With Julie Meunier.
- `extra_auxiliary_var_names` reads unfiltered variables into the offline filter, e.g. for the frequency scale.
- Offline input data can be read from NetCDF files.

### Changes and fixes

- The Eulerian filter in post-processing (`compute_Eulerian_filter!`) is computed with running sums, so its cost grows linearly with the number of saved times, and it treats the ends of the data as the Lagrangian filter does. A field that is constant in time is now unchanged by it. It also supports separate outputs, and a frequency scale that depends only on position. The reference test data was regenerated for this and for the new regridding.
- The offline filter compiles its output writer once, rather than again for the backward pass.
- The code is reorganised: the utilities are split into files by purpose (filter parameters, variables, forcing, outputs, initialisation, post-processing and regridding), the forcing of each filtered variable is a single `DiscreteForcing`, and input data is read in a new `DataIO` module, which now holds `update_input_data!` (still exported).
- The integration tests also run on a GPU, if there is one, and CI tests on the LTS and latest Julia.

### Bug fixes

- The backward pass could use the forward pass's input data when there were few saved times in the filter window.
- High-order advection schemes (e.g. `WENO(order=9)`) no longer fail with a `DimensionMismatch`.
- Offline runs without maps (no regridding or mean velocities) no longer fail at the start of the backward pass.
- Passing consistent `T_start`, `T` and `T_end` no longer throws an `UndefVarError`, and `T_start + T` is compared with `T_end` to rounding error.
- When `T_start` or `T_end` is between saved times, the input data is interpolated from the saved times either side.
- The default `output_filename` of `OnlineFilterConfig` now ends in `.jld2`, the name of the file that `JLD2Writer` writes, so that post-processing finds it.
