# Initial dev2026 reorganization validation

Validation date: 3 November 2026. Baseline: `0a3f9dc`.

## Scope

Moved 35 root-level C++ source/header files into `src/`. Each relocated file
was compared byte-for-byte against its original Git blob and was identical.
Updated build rules, companion utility links, and test includes for the new paths.
Added layout documentation and ignores for generated build/test artifacts.
No statistical model, default parameter, or expected checksum was changed.

## Local environment

Validation ran on the local Mac using installed Conda `bam-filter` dependencies
(HTSlib 1.21, GSL 2.7) and Apple clang 21 targeting x86_64 to match those libraries.
The inherited Conda clang 17 failed against the installed macOS SDK on both
baseline and reorganized sources. Switching compilers resolved that build issue.

Both versions built with the equivalent of:

```sh
make HTSSRC=systemwide PREFIX=/path/to/bam-filter \
  CXX=/usr/bin/clang++ CC=/usr/bin/clang \
  CXXFLAGS='-O2 -std=c++11 -arch x86_64 -MMD -MP' \
  CFLAGS='-O2 -arch x86_64 -MMD -MP' -j2
```

The test script directly invokes `c++`, so local validation used a temporary
compiler wrapper with the same architecture and an rpath to the dependency
library directory. CPATH and LIBRARY_PATH pointed at that environment's include
and lib directories. This wrapper is environment setup, not a repository change.

## Results

- Baseline and reorganized application and companion utilities built successfully.
- All 35 relocated files were byte-identical to their baseline source blobs.
- Taxonomy sorting regression passed.
- DUST and PMD unit tests passed with the corrected library lookup.
- All 18 files listed in `test/output.md5` matched between baseline and
  reorganized runs, including the input fixture and the global fit output.
- The full suite did not pass its stored checksums: both versions failed the
  same `output/test_dfit_global.dfit.fix` checksum. All other listed checksums
  passed. This is a pre-existing local discrepancy; its numerical cause has
  not been diagnosed. Expected checksums were preserved.
- `git diff --check` passed.

The first full runs also hit missing-rpath errors in the two unit binaries;
those were resolved in the local test compiler configuration. These checks
support a behavior-preserving relocation, not a fresh validation of scientific
accuracy or of optional regression builds.

The final full-suite rerun with the corrected compiler wrapper completed with
only the same global-fit checksum failure; the unit tests and other checks passed.
