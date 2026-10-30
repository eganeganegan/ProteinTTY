# Workflow validation

Run `cmake -S . -B build && cmake --build build && ctest --test-dir build --output-on-failure`.
Both test executables run offline and are included in the existing Linux GCC/Clang, macOS, and Windows CI matrix.

## Numerical references

`tests/workflow_tests.cpp` checks quaternion least-squares fitting against exact rigid transforms,
including a 180-degree rotation, a one-atom selection, and a chiral reflection that a proper rotation cannot fit.
The selected anchor and reference fit radius remain stable when an unselected atom moves.
The tolerance is 1e-8 Angstrom for these small double-precision fixtures.

The RMSD fixture uses positions (0,0,0)/(1,0,0), (0,0,0)/(3,0,0), and (0,0,0)/(5,0,0).
The second frame's unaligned RMSD is sqrt(2) Angstrom over both atoms.
The fixed atom's RMSF is zero and the moving atom's population RMSF is sqrt(8/3) Angstrom.
Sampling first and last frames gives an RMSF of 2 Angstrom. Rigidly transformed frames have
aligned RMSD and RMSF below 1e-8 Angstrom. Results use equal atom/frame weights, without mass weighting,
periodic unwrapping, or unit-cell image correction. Apply unwrapping before loading wrapped simulations.

XYZ/DCD coordinates remain disk-indexed. Analysis retains a working frame and alignment reference;
RMSF uses per-selected-atom Welford accumulators and RMSD stores scalar results only.
Loading, stepping, loop boundaries, playback ranges, stride, and malformed inputs have separate tests.

## Molecular and map fixtures

PDB, mmCIF, and BinaryCIF fixtures test model selection, insertion codes, differing author/label chain IDs,
explicit long-distance connections, and assembly transformations. CIF tests include singleton categories,
operator ranges and products. PDB/mmCIF export/reload checks atom identifiers, selected alternate markers,
coordinates, assembly copies, and connections; selection export drops bonds to excluded atoms.
PDB coordinates use a 0.001 Angstrom round-trip tolerance; mmCIF writes six decimals.

Generated MRC fixtures cover little/big endian modes 0/1/2/6, column/row axis permutation, extended-header
skipping, start indices, explicit origin, truncation, and short headers. Axis/origin checks use a 1e-8
Angstrom tolerance. Complex, packed, and half-float modes are rejected. Density maps must carry a MAP signature.
Map origin and handedness conventions vary by producer; use maps in the same deposited coordinate system
as the loaded structure. Crystal-mate struct_conn links are excluded. Only deposited biological assembly operators are expanded;
arbitrary crystallographic mate generation is not implemented.

Sessions restore identical plain-text frames including labels and map slices, trajectory alignment, colors,
selections, and camera state. Invalid sessions preserve the current view. Referenced inputs must remain
available and unchanged. Macro tests exercise command dispatch and bounded recursion.

## Reproducible performance check

```
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release
cmake --build build-release --target workflow_benchmark
./build-release/workflow_benchmark
```

The benchmark finalizes and indexes 100,000 atoms, then renders ten rotating 120x36 Braille frames, with
and without a slice of a 256³ float density map (64 MiB). It reports setup and mean per-frame times;
there is no machine-dependent timing assertion in unit tests.

On the local macOS host, the initial implementation measured approximately 54 ms/frame for the structure
and 59 ms/frame with the map slice. Cached atom silhouette geometry reduced these to approximately
31 ms/frame for both workloads. These are local measurements, not cross-platform guarantees.
The existing 1crn snapshot was byte-identical before and after the sprite optimization.

The decoded MRC value limit is 512 MiB; slice rendering uses a uniform stride above approximately
200,000 samples. Assembly expansion is capped at 10,000 operator products and 5,000,000 atoms.
