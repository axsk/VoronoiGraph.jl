# Release Notes

## v0.2.3

### Breaking

- Julia compat raised to 1.9 (required for package extensions).

### Changes

**Polyhedra moved to a package extension**

- `Polyhedra.jl` is no longer a hard dependency. It is now a weak dependency loaded via the `PolyhedraExt` extension when `Polyhedra` is explicitly imported by the user.
- `Polyhedra` compat extended to include 0.8.
- `boundary_area()` now throws an informative error if called without `Polyhedra` loaded; the real implementation lives in the extension.
- `boundary_area_vrep()` removed (was also Polyhedra-dependent).

**Runtime dispatch fixes (`voronoi.jl`, `raycast.jl`)**

- `u_default = u_qr` made `const` — fixes a runtime dispatch that prevented the compiler from specializing calls to `u_default(...)`.
- New module-level globals (`warn_degenerate`, `USE_HEURISTIC`, `cache`, `debug`, `hit`, `miss`, `new`) are declared with type annotations (`::Bool`, `::Int`), avoiding runtime dispatch when they are accessed inside hot loops.

**`raycast.jl`**

- Degenerate vertex warnings are now opt-in via `VoronoiGraph.warn_degenerate = true` (previously always emitted).
- `raycast_start_heuristic()` extracted into its own function; togglable via `VoronoiGraph.use_heuristic(false)`.
- Reverted to the simpler (and correct) start heuristic after benchmarking showed the "more correct" variant was slower without meaningful accuracy gain.
- Fixed the degenerate warning message: it now reports the new candidate generator `[j]` instead of `[i]`.
- `local t` declaration added before the convergence loop to fix variable scoping.

**`voronoi.jl`**

- `yield()` removed from the BFS loop in `explore()` — recovers ~10% runtime (it was only present to allow Ctrl+C interrupts).
- `Progress(...)` updated to keyword-argument form (`dt=1, desc=`) to fix a ProgressMeter deprecation warning.
- Added `report_memory()` helper and a `debug` flag; when `VoronoiGraph.debug = true`, `explore()` prints entry counts and memory usage of internal dictionaries after the BFS completes.
