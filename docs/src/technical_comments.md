# Technical comments

- The online and offline filtering both work on GPU, just make sure CUDA is installed in the environment and loaded, and set `architecture = GPU()`
- The post-processing interpolation step using [`regrid_to_mean_position!`](@ref "regrid_to_mean_position!") isn't currently set up to work in domains with immersed boundaries or with more than one bounded dimension. Interpolation results might be junk near the boundary!
- For the offline filter, saved tracers to be filtered need to be on `Center` grid. Velocities should be on their standard grid. 
- As with any moving average, there are endpoint effects. For the offline filter, a time window at each end of the filtered timeseries of the order of the inverse of the cutoff frequency should be excluded from any further analysis. For the online filter, this is only necessary at the beginning of the timeseries. 