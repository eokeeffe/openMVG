# GPS / Visual / Hybrid Atomic Clustering

`openMVG_clustering` is the library behind
`openMVG_main_ComputeAtomicClustering` (see
`src/software/SfM/clustering/main_ComputeAtomicClustering.cpp`): it splits a
set of images into overlapping clusters so each can be reconstructed
independently and later merged, using:

- **easyexif** (already vendored at `src/third_party/easyexif`) -- reading
  GPS coordinates from image EXIF data
- **LEMON** (already vendored at `src/third_party/lemon`) -- the
  k-nearest-neighbor graph and Dijkstra shortest-path search used to build
  the "molecule" (the overlap-inducing bridge cluster)
- **OpenCV** -- ORB feature extraction and k-means, for the visual and
  hybrid methods only

All three methods share one pipeline (`atomic_clustering.cpp`), generalized
over a `DistanceFn` -- swap in GPS distance, visual dissimilarity, or a
blend of both, and the atom/molecule logic itself doesn't change.

## Method summary

- **GPS-only**: distance = GPS haversine distance
- **Visual-only**: distance = ORB + visual-vocabulary dissimilarity (GPS is
  only used afterward, to plot each cluster on the map -- not for
  clustering itself)
- **Hybrid GPS/visual**: distance = weighted blend of both (`--hybrid_alpha`)

Adapted from Xie et al., "Hierarchical Clustering-Aligning Framework Based
Fast Large-Scale 3D Reconstruction Using Aerial Imagery" (Remote Sensing,
2019): atoms are a hard partition grown outward from well-spaced seeds;
the molecule is built from shortest paths between adjacent atoms in a
k-NN graph, and its images stay members of their original atom too --
that's where the overlap comes from.

## Building

The GPS-only method has no extra dependencies beyond what OpenMVG already
vendors, so `openMVG_clustering` and `openMVG_main_ComputeAtomicClustering`
build by default. The visual and hybrid methods additionally need OpenCV;
enable them with:

```bash
cmake -DOpenMVG_USE_OPENCV=ON ..
```

If OpenCV isn't found (or `OpenMVG_USE_OPENCV` is left off, the default),
the library and tool still build, but `--method visual` and `--method
hybrid` are rejected at runtime -- only `--method gps` is available.

## Usage

See `openMVG_main_ComputeAtomicClustering --help` for the full option list
(camera sensor database lookup, GPS pose priors, per-cluster `sfm_data`
output, and the Leaflet cluster map). Quick examples:

```bash
openMVG_main_ComputeAtomicClustering -i /path/to/photos -o /path/to/out \
    --method gps --min_atom_size 15 --max_atom_size 50

openMVG_main_ComputeAtomicClustering -i /path/to/photos -o /path/to/out \
    --method visual --vocab_size 300 --orb_features 500

openMVG_main_ComputeAtomicClustering -i /path/to/photos -o /path/to/out \
    --method hybrid --hybrid_alpha 0.5
```

## What's been validated, and what hasn't

I don't have network access in the environment this was built in, so I
could not install OpenCV or LEMON to do a full real build of the visual
pipeline against real ORB/k-means code (LEMON is now covered, since
OpenMVG always vendors and builds it internally). What I *did* do instead:

- Wrote hand-built stub headers that mimic the real public API surface of
  OpenCV (`opencv2/*.hpp`) and ran `g++ -fsyntax-only` against every `.cpp`
  file with those stubs -- this caught and fixed one real bug (an OpenCV
  macro/namespace mismatch) before you'd have hit it.
- Compiled and ran a synthetic-data functional test directly against
  `atomic_clustering.cpp` (now formalized as `atomic_clustering_test.cpp`'s
  `GpsOnly_*` case) -- confirming the atom/molecule logic (farthest-point
  sampling, region growing, k-NN graph construction, shortest-path molecule
  generation) produces correct results: multiple atoms, a bridging
  molecule, real overlap, full coverage.

What this does **not** validate: the real OpenCV headers might differ in
some detail from my stubs (a renamed field, a slightly different template
signature, etc.), and `visual_features.cpp` was only syntax-checked, not
functionally tested against real ORB feature extraction. `atomic_clustering_test.cpp`'s
`VisualOnly_*`/`Hybrid_*` cases (built only when `ATOMIC_CLUSTERING_HAS_VISUAL`
is set) are the first real exercise of that path -- if they fail to compile
or pass against the real OpenCV on your machine, paste the error back and it
can be fixed quickly, since the logic itself has already been proven out.
