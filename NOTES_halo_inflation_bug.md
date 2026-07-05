# WENO9 halo-inflation DimensionMismatch — investigation notes

Triggered by: `JAMESpaper/code/data_creation/offline_filter_modified_shallow_water_wave_turbulence_WENO9.jl`,
using `advection = WENO(order=9)` with the offline filter on `main` (v0.1.3, what JAMESpaper's
Manifest currently pins).

## Root cause

`OfflineFilterConfig` infers `config.grid` from the original data file, inheriting whatever halo
the source simulation was saved with. Here that's `RectilinearGrid{Periodic, Periodic, Flat}`,
`Hx=Hy=3, Hz=0` (confirmed by reading `code/data/MSW_wave_turbulence_IO/MSW_wave_turbulence_IO.jld2`
directly — grid and BCs on `ω`/`u` are periodic in x and y).

When `LagrangianFilter` is constructed with `advection=WENO(order=9)`, `inflate_grid_halo_size`
(`src/OfflineLagrangianFilter/lagrangian_filter.jl:185`) detects WENO9 needs halo ≥5 and replaces
the grid via `with_halo((5,5,0), grid)`. Model tracer/velocity/auxiliary fields now live on a
halo-(5,5,0) grid.

But `load_data` (`src/Utils/lagrangian_filter_utils.jl:201`) builds `FieldTimeSeries` straight from
the `_filter_input.jld2` file, which keeps the *original* halo (3,3,0) — never inflated to match.

Two call sites then do `parent(field) .= parent(data_field)`, broadcasting a (5,5,0)-halo array
against a (3,3,0)-halo array → `DimensionMismatch`:

1. **`initialise_filtered_vars_from_data`** (`lagrangian_filter_utils.jl:1086`, tracer C/S init +
   map init) — this is where the reported error actually happens first.
2. **`update_input_data!`** (`lagrangian_filter_utils.jl:1036`) — same pattern for
   `model.velocities` and `model.auxiliary_fields`, run every timestep via an
   `UpdateStateCallsite` callback. Not yet hit by the reported traceback (crashes before this
   runs), but will hit the identical error next once (1) is fixed.

## Fix for (1): safe in general, not just for this periodic case

`update_state!(model::LagrangianFilter)` (`update_lagrangian_filter_state.jl:32-34`) already calls
`fill_halo_regions!` on `model.tracers` every step, unconditionally, using each tracer's own BCs.
So the halo values written directly during initialisation get overwritten by the standard fill
almost immediately anyway. That means swapping `parent(...) .= parent(...)` for
`interior(...) .= interior(...)` in `initialise_filtered_vars_from_data` is safe regardless of
domain topology — sizes always match (interior dims are `Nx,Ny,Nz`, unaffected by halo), and the
halo gets correctly filled a moment later by the existing mechanism either way.

**Not yet applied on this branch.** (A version of this edit was prototyped on `buffered-data-reader`,
in its stash, but that branch's `NamedTuple`-based functions are dead code there — that branch
already fully switched to `BufferedDataReader`. The live code for JAMESpaper is here on `main`.)

## Fix for (2): genuine design tension — NOT just a periodic-vs-general question

`update_state!` deliberately does **not** call `fill_halo_regions!` for velocities/auxiliary
fields (explicit comment at `update_lagrangian_filter_state.jl:32`: "we're only applying the BCs
to the tracers - not the velocities or the auxiliary fields"). `update_input_data!` fills their
halos by copying real saved data instead (comment at `lagrangian_filter_utils.jl:1042-1043`,
`1053-1054`: "This also fills the halo regions, which we'll need to help with the filtered field
boundaries").

Why this matters: for THIS case (doubly periodic), a periodic BC-based `fill_halo_regions!` would
be mathematically identical to copying real data — periodic data repeats exactly. But the package
is also used for non-periodic/open-boundary domains elsewhere in JAMESpaper (e.g.
`offline_filter_lee_wave.jl`). For those, a generic BC-based halo fill is NOT equivalent to the
real off-domain snapshot — it would silently substitute a boundary condition for actual data,
which could quietly change results rather than crash loudly.

So simply switching `update_input_data!` to `interior(...)` + generic `fill_halo_regions!` is only
safe when the (possibly-inflated) grid is fully periodic in the directions where halo got inflated.
A general-purpose fix needs to either:
- check the grid topology and branch (periodic → BC-fill is fine; anything else → error clearly
  rather than corrupt silently), or
- avoid the situation entirely for non-periodic cases by ensuring source data is saved with a
  halo ≥ whatever advection scheme will ultimately be used offline, or
- something else — not decided yet.

Three options were on the table when this got parked (none applied yet):
1. **Periodic-only fix**: branch on grid topology; BC-fill only when periodic in the inflated
   directions, error otherwise.
2. **Always BC-fill**: simplest, correct here, silently wrong for non-periodic cases.
3. **Don't touch `update_input_data!`**: instead regenerate `MSW_wave_turbulence_IO.jld2` (and any
   other source data destined for high-order offline advection) with halo ≥5 from the start, so
   inflation past the data's native halo never happens.

## Status

Nothing on this branch (`fix/weno-halo-inflation-mismatch`, based on `main`/v0.1.3) has been
changed yet — this file is the only commit. Revisit and pick one of the above before touching
`update_input_data!`. The `initialise_filtered_vars_from_data` fix (interior instead of parent) is
low-risk and can likely be applied first independently.
