# Development guide

The `dev2026` branch starts from `0a3f9dc`. The first reorganization moves
the core sources into `src/` without changing the statistical algorithms,
command names, output schemas, or runtime resource locations.

## Repository layout

| Location | Purpose |
| --- | --- |
| `src/` | Core C++ sources and headers, including existing bundled numerical code |
| `misc/` | Companion utilities and experimental fitting programs |
| `R/` | Existing plotting and analysis scripts; some require dataset-specific paths |
| `test/` | Unit checks, integration scripts, fixtures, and expected checksums |
| `doc/` | Format reference and developer documentation |
| `.github/workflows/` | CI and release packaging |
| `rank2level.txt` | Runtime taxonomy rank mapping; kept at the root for compatibility |
| `Makefile` | Build entry point; executable remains `./metaDMG-cpp` |

Generated objects and dependency files live beside their sources in `src/`.
The generated `version.h` remains at the root. The companion utilities link
against `src/shared.o`; the test scripts use the same source layout.

## Code map

- `metaDMG.cpp`: command dispatch.
- `main_getdamage.cpp`, `profile.cpp`: alignment processing and mismatch counts.
- `ngsLCA.cpp`, `ngsLCA_cli.cpp`: taxonomic classification and CLI options.
- `main_dfit.cpp`, `dfit_optim.cpp`, `pval.cpp`: model fitting, resampling, and statistics.
- `main_pmd.cpp`: per-read PMD scoring.
- `main_aggregate_stat.cpp`: taxonomic aggregation.
- `main_print.cpp`, `main_filter_bdamage.cpp`: output and filtering.
- `regression.cpp`: optional multinomial regression, guarded by `__REGRESSION__`.
- `shared.cpp`, `mrand.cpp`: shared utilities and random number generators.
- `bfgs.cpp`, `matplotlibcpp.h`: existing numerical/plotting support; preserve provenance.

## Build and verification

Run from the repository root:

```sh
make
make test
```

With installed HTSlib and dependencies, use `make HTSSRC=systemwide` and,
if needed, `PREFIX=/path/to/dependency/environment`. Use a compiler and
libraries built for the same architecture. The current test script invokes
`c++` directly, so its compiler environment must also match the build.

`make test` expects the main executable to have been built. It runs the taxonomy
sorting regression followed by `test/testAll2.sh`. Tests generate files in
`test/`; `make testclean` removes the integration outputs. The existing
`make clean` also removes the configured HTSlib source directory when not
using `HTSSRC=systemwide`.

For a structural change, compare against the starting commit using identical
fixtures, seeds, compilers, and libraries. Do not revise expected checksums just
to make a reorganization pass.

## Research workspace

The sibling `metaDMG/` directory is a research archive, not part of this Git
repository. Its guide and inventory describe the manuscript, simulation, and
figure material. Changes there are not protected by this branch. Keep raw
results and Google Docs pointers in place until script dependencies and source
provenance are mapped.

## Follow-up work

1. Reconcile documented CLI defaults with initialization and help output.
2. Specify the likelihood, bootstrap methods, and both Zfit definitions.
3. Make plotting scripts accept paths and parse output columns by name.
4. Separate reusable plotting tools from manuscript-specific analyses.
5. Establish reproducible benchmarks for the manuscript claims.

These are follow-up tasks, not changes included in the source relocation.
