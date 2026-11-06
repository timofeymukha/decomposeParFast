# Faster decomposePar for large meshes (OpenFOAM)

decomposePar gets slower with every processor it writes, and on large meshes it becomes the bottleneck. This patch removes the causes. The decomposition it writes is unchanged.

## Contents

- `0001-ENH-decomposePar-remove-O-nProcs-meshSize-work-and-c.patch` — for the openfoam.com `develop` branch (on top of 9bcf41d9, 2026-11-04). Apply with `git am`.
- `decomposePar-v2606.patch` — the same change for v2606.
- `decomposeParFast/` — the patched v2606 sources, which build as a separate executable `decomposeParFast`.

## Changes

All changes are in `applications/utilities/parallelProcessing/decomposePar`.

1. **Lookups allocated once.** `writeDecomposition()` allocated and filled `pointLookup(nPoints)` and `faceLookup(nFaces)` for every processor. These are now allocated once, together with the zone markers, and only the entries a processor touched are reset afterwards.
2. **Point distribution no longer scans the whole mesh.** `decomposeMesh()` allocated a `bitSet(nPoints)` per processor and called `sortedToc()`, which scans the whole set. The bitSet is now allocated once; points are collected as they are first marked, then sorted.
3. **Processor meshes rebuilt in memory, not read back.**
   - **Before:** the field stage read every processor mesh and its `*ProcAddressing` back from disk. In a collated file, reading block $i$ skips blocks $0..i-1$ by reading them (`istream::ignore`). So processor $i$ costs $O(i)$, and the stage is quadratic.
   - **Now:** the processor `fvMesh` and its addressing are rebuilt in memory with the same code that writes them, moved into `buildProcMesh<MeshType>()`.
   - **Fallbacks:** the read-back is still used with `-fields`, and when the points come from a different time instance. The new option `-readProcMeshes` forces it.

Items 1 and 2 cost $O(N_\text{mesh})$ per processor, so $O(N_\text{mesh}\cdot N_\text{procs})$ in total. Item 3 is $O(N_\text{procs}^2)$.

## Timings

Stock vs patched decomposePar. Hierarchical method, collated format, serial runs.

| Case | Labels | Stock | Patched |
|---|---|---|---|
| plane channel, 2.7M cells → 1024 | Int32 | 31.7 s | 5.0 s |
| plane channel, 670k cells → 4096, field stage | Int32 | 44 s | 1.3 s |
| channel, 46.4M cells → 3072 | Int64 | **1:53:21** | **2:56** |
| channel, 528.5M cells → 34560 | Int64 | stopped after 20 min, still distributing points; estimated about 4 days | **1:30:49** (377 GB), measured without change 2, which accounted for 565 s |

Per-processor cost for the 46.4M-cell case:

| Stage | Stock | Patched |
|---|---|---|
| processor meshes | 0.91 s per processor | 20 ms per processor |
| fields | $17\,\text{ms} + 0.88\,\text{ms}\cdot i$ for processor $i$, rising to 2.3 s | 16.5 ms per processor |

For the 46.4M-cell case, the output (12 collated files, 13 GB) is md5-identical to stock.

## Verification

Every tutorial in `develop` that calls decomposePar was run with stock and patched decomposePar.
- **Method:** each pair of outputs was compared file by file. Each call ran twice: as the tutorial calls it, and with `-fileHandler collated`.
- **Coverage:** 514 comparisons from 239 tutorials, including `-region`, `-allRegions`, `-cellDist`, `-force`, `-constant`, `-latestTime`, `-fields`, `-dry-run` and finite-area cases.
- **Identical:** 486 comparisons.
- **Both failed the same way:** 18 comparisons.
- **Different:** 10 comparisons, from 5 tutorials, and none are caused by the patch. In each of them, stock itself is not reproducible:
  - **`pointLevel` is read out of bounds.** This is a pre-existing bug: after `createBaffles`/`mergeOrSplitBaffles -split`, `pointLevel` is shorter than `points`. All in-range entries are identical between stock and patched.
  - **`dynamicCode/`** contains compiled coded boundary conditions, which embed the case path.
  - **`-dry-run -cellDist`:** only the case name in the VTK header comment differs.

The patch compiles on `develop` without new warnings.
