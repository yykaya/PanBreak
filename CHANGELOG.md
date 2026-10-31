# Changelog

## 0.1.1 (2026-10-31)

### Added
- `panbreak check`: reports the Python dependencies and external tools that PanBreak can find, with
  versions and the subcommands that need them.
- End-to-end test on a toy graph and an unaligned BAM (`bubbles` → `kmers` → `genotype`), run in CI
  and in the conda package test.
- CI job that installs the conda environment on Linux and macOS (arm64) and builds the conda package
  from the repository.

### Changed
- GFF files with GTF-style attributes (`key "value"`, e.g. liftoff output) are accepted by
  `regions` and `annotate`.
- Insertion clustering across genomes; bubble labels prefer protein-coding genes.
- Graph layout and haplotype flow plots scale to large panels.
- Python 3.10 or newer is required; tested with 3.10, 3.12 and 3.13.
- Package metadata uses an SPDX license expression.

## 0.1.0 (2026-10-22)

First public release.
