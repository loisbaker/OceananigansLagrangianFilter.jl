# Regridding to the mean position

The filter finds ``f^*(\vb*{x},t)``, the average along the trajectory of the particle that is at ``\vb*{x}`` at time ``t`` (see [Eulerian and Lagrangian averaging](@ref)). The Lagrangian mean ``\bar{f}^{\mathrm{L}}`` belongs instead at that particle's mean position, ``\vb*{\Xi}(\vb*{x}) = \vb*{x} + \vb*{\xi}(\vb*{x})``, where ``\vb*{\xi}`` are the maps computed by the filter. [`regrid_to_mean_position!`](@ref) moves ``f^*`` from the grid points, which label the particles, to their mean positions, and interpolates it back onto the grid. This is done separately at each time, so ``t`` is dropped below:
```math
\begin{equation}
    \bar{f}^{\mathrm{L}}(\vb*{\Xi}(\vb*{x})) = f^*(\vb*{x})\,, \qquad \text{so} \qquad \bar{f}^{\mathrm{L}}(\vb*{y}) = f^*(\vb*{\Xi}^{-1}(\vb*{y}))\,.
\end{equation}
```

## How the regridding works

Carrying every grid point to its mean position gives a deformed mesh with the same connections as the grid, so no triangulation is needed. For each grid point ``\vb*{y}``:

1. the position ``\vb*{x} = \vb*{\Xi}^{-1}(\vb*{y})``, whose mean position is ``\vb*{y}``, is found by solving ``\vb*{x} + \vb*{\xi}(\vb*{x}) = \vb*{y}`` with Newton's method, with the maps ``\vb*{\xi}`` interpolated linearly between cell centres;
2. if ``\vb*{x}`` is in the fluid (inside the domain, in a fluid cell, and above the bottom of an immersed boundary), ``\bar{f}^{\mathrm{L}}(\vb*{y}) = f^*(\vb*{x})``, with ``f^*`` interpolated in the same way;
3. otherwise, no particle in the fluid has its mean position at ``\vb*{y}``, and ``\bar{f}^{\mathrm{L}}(\vb*{y})`` is `NaN`.

Beyond a bounded boundary, the maps and fields are extrapolated linearly into one ghost cell, so nothing is assumed about the flow at the boundary; periodic dimensions wrap around. Both steps are second order in the grid spacing, and the same kernel runs on CPU and GPU.

## When the result is valid

**The maps must be resolved by the grid.** Linear interpolation of ``\vb*{\xi}`` assumes it varies smoothly from cell to cell. Where neighbouring particles' mean positions cross over, a deformed cell is turned inside out and the mesh *folds*: some grid points are then the mean position of more than one particle, and the regridded value there is not unique. This is detected (the Jacobian determinant of the map changes sign) and reported as a warning. Folds typically come from strong shear that is resolved by only a few cells, for example next to a boundary.

**Away from the ends of the filter window.** Near the start and end of the data, the filter only sees one side in time, so the maps there do not represent the mean position of the full trajectory, and the regridded fields are not the Lagrangian mean.

**`NaN` where no particle has its mean position.** This is not a failure: it marks points that are not the mean position of any particle in the fluid. It typically happens next to boundaries, where particles' mean positions are pulled away from the boundary. For example, particles near a flat bottom that ride over a hill within the filter window have mean positions above the bottom, so the bottom row is `NaN`. Points where Newton's method does not converge, usually at the edge of a fold, are also `NaN`, and counted in a warning.

**Next to an immersed boundary the result is approximate.** Solid cells next to the fluid are filled with the mean of their fluid neighbours, so the first fluid cells use maps and fields extrapolated from the fluid. Labels, and grid points, are in the solid if they are in a solid cell, or below the bottom height where the immersed boundary has one (`GridFittedBottom`, `PartialCellBottom`), so grid points below the bottom in a partial cell are `NaN` too.

## On a `LatitudeLongitudeGrid`

The mean flow map is the filtered position of each particle, in whatever coordinates the grid uses (see [Eulerian and Lagrangian averaging](@ref)). The filter finds the maps by filtering the velocity along trajectories, so it finds the mean flow map exactly only if the rates of change of the coordinates are the velocity components. That holds in Cartesian coordinates. On a `LatitudeLongitudeGrid`, the maps ``\xi_u`` and ``\xi_v`` are instead the eastward and northward distances travelled (in metres), and they are converted to displacements in longitude and latitude at the latitude ``\phi`` of each particle's current position,
```math
\begin{equation}
    \Delta\lambda = \frac{\xi_u}{R\cos\phi}\,, \qquad \Delta\phi = \frac{\xi_v}{R}\,,
\end{equation}
```
where ``R`` is the radius of the sphere. This assumes the grid is locally Cartesian over the trajectory excursions, that is, that the excursions are small compared with the coordinate system's radius of curvature.

The mean position on a sphere is itself not uniquely defined. The generalised Lagrangian mean of [andrewsExactTheoryNonlinear1978](@citet) assumes a Euclidean (Cartesian) structure, and on a curved surface it depends on the choice of coordinates; coordinate-free definitions of the mean flow map are possible, but involve a choice [gilbertGeometricGeneralisedLagrangianmean2018](@citep).
