# ProteinTTY Roadmap

ProteinTTY v0.11 covers the core viewer: PDB, mmCIF, and BinaryCIF structures; RCSB loading; interactive selections and measurements; multiple terminal renderers; secondary-structure assignment; XYZ trajectories; native POSIX and Windows consoles; hardened input handling; and automated cross-platform releases.

The next releases should deepen correctness and scientific workflows rather than add more display modes. Versions indicate intended order, not dates.

## Completed Milestones

### v0.1: Initial Prototype

- Established the C++20/CMake application, molecular data model, PDB ingestion, and interactive terminal loop.
- Introduced the first software-rendered atom and bond views with camera controls and ANSI output.

This prototype predates the repository's retained release history; v0.2 is the first tagged commit.

### v0.2: Complete Core Viewer

- Added alternate-conformer handling, secondary-structure records, automatic fitting, depth buffering, lighting, and multiple color and representation modes.
- Added ASCII, Unicode, and Braille renderers; metadata and help panels; inspection, sequence, snapshots, autorotation, RCSB fetching, and safe terminal restoration.

### v0.3: mmCIF Support

- Added tag-indexed mmCIF parsing with quoted and multiline values, flexible column ordering, and deposited secondary-structure ranges.
- Added extension and content-based format detection across PDB and mmCIF inputs.

### v0.4: Large-Structure Scaling

- Replaced all-pairs bond candidate searching with a spatial index.
- Added full string chain identifiers for modern mmCIF assemblies and improved large-structure behavior.

### v0.5: Command Console and Automation

- Added the interactive `:` command console with history, quoted arguments, chained commands, help, and query output.
- Added repeatable `--command` execution and reusable `.ptty` scripts through the same dispatcher.

### v0.6: Molecular Selections and Analysis

- Added boolean atom selections, named selections, visibility masks, per-selection coloring, and selection-aware centering, zooming, and inspection.
- Added structure queries and geometric measurements, command completion and positioned errors, plus runtime load, fetch, and snapshot commands.

### v0.7: Secondary-Structure Fallback

- Added chain-local C-alpha geometry inference for unannotated protein helices and sheets.
- Kept deposited PDB/mmCIF annotations authoritative and corrected selection/rendering edge cases.

### v0.8: Trajectory Playback

- Added strict multi-frame XYZ coordinate sources with topology and atom-label validation.
- Added play, pause, resume, stepping, speed, looping, keyboard controls, scripting, and snapshot support for trajectories.

### v0.9: Native Windows Support

- Added a Win32 console backend with virtual-terminal output, native input events, UTF-8 handling, resizing, and RAII restoration.
- Made RCSB process launching and cache paths portable, added MSVC warning settings, and retained redirected snapshot output.

### v0.10: BinaryCIF Support

- Added dependency-free MessagePack and BinaryCIF decoding for numeric, string, packing, delta, run-length, fixed-point, quantized, and mask encodings.
- Added local `.bcif` detection, RCSB Models fetching, command/CLI integration, secondary annotations, offline fixtures, and large-structure smoke tests.

### v0.11: Reliability and Portability

- Preserved insertion codes in PDB, mmCIF, and BinaryCIF secondary ranges, with exact boundary tests for residues such as `100`, `100A`, and `100B`.
- Added bounded BinaryCIF document, container, encoding, and decoded-value handling with adversarial malformed-input tests.
- Added a platform-neutral terminal decoder covering pasted input, split ANSI arrows, repeats, UTF-8 code points, Escape, and resize events.
- Added Linux GCC/Clang, macOS AppleClang, and Windows MSVC CI with Release tests and deterministic snapshots.
- Added tag-triggered Linux, macOS, and Windows release archives, individual SHA-256 files, a consolidated checksum manifest, and a locally runnable packaging target.

## v0.12: Trajectory Workflows

- **Implemented:** Random-access DCD coordinates using the existing topology supplied by PDB, mmCIF, or BinaryCIF.
- **Implemented:** Disk indexing for XYZ and DCD trajectories without retaining every frame in memory.
- Add least-squares trajectory alignment to a reference frame for a selected atom set.
- Add RMSD and RMSF calculations with selection-aware summaries and coloring.
- Add playback range, stride, and frame-number controls for long simulations.

**Complete when:** trajectories larger than available memory can be stepped in either direction, aligned playback is stable, and RMSD/RMSF results match reference fixtures within documented tolerances.

## v0.13: Molecular Semantics and Export

- Read explicit connectivity from PDB `CONECT` and mmCIF/BinaryCIF `struct_conn`, preferring deposited bonds over inferred bonds.
- Add biological-assembly expansion from deposited symmetry and operator tables.
- Support selecting and switching among PDB/mmCIF models and NMR ensembles.
- Export the whole structure or a selection to PDB and mmCIF while preserving identifiers, insertion codes, alternate locations, and coordinates.
- Add residue and chain labels that avoid collisions in dense terminal views.

**Complete when:** explicit connectivity and assembly/model selection survive load-render-export-reload round trips, with no topology changes outside the selected operation.

## v0.14: Map and Session Support

- Load CCP4/MRC density maps and render configurable isosurfaces or slices alongside atomic structures.
- Save and restore viewer sessions, including camera, rendering settings, named selections, trajectory position, and loaded paths.
- Add command macros and optional key bindings in a user configuration file.
- Profile and optimize rendering for very large assemblies and density maps before adding further visual effects.

**Complete when:** a saved session reproduces the same deterministic frame, map transforms align with deposited structures, and representative large-data benchmarks stay interactive.

## Guiding Constraints

- Keep local structure viewing dependency-free unless a format has no reasonable compact implementation.
- Preserve the shared command dispatcher across interactive, scripted, and command-line use.
- Keep topology separate from coordinate frames and avoid renderer-specific data in parsers.
- Require offline fixtures for parsers and analysis algorithms; network tests remain smoke tests only.
- Treat Windows, macOS, and Linux as equal release targets.
