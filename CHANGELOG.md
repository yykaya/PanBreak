# Changelog

## Unreleased

### Added
- `panbreak graph --builder minigraph`: locus graphs built natively with minigraph (`-cxggs`), with
  assembly paths from `minigraph -cxasm --vc`; runs on macOS without containers.
- `panbreak graph --builder cactus`: locus graphs built with Minigraph-Cactus (`cactus-pangenome`) in a
  digest-pinned container or a native install.
- `panbreak convert`: imports GFA with `P` lines, `W` lines (Minigraph-Cactus, vg) or minigraph rGFA
  with GAF alignments, and extracts a locus from a large graph with odgi (`--region`).
- `--skip-build` (alias of `--skip-pggb`) reuses the output of any builder.
- `panbreak check` lists minigraph and cactus-pangenome and the available graph builders.

### Changed
- `bubbles` reverses paths that run against the reference orientation.
- `bubbles` ignores path events within `--cover-pad` bp of a path end (unanchored sequence ends);
  results on pggb graphs of the example are unchanged.
- `--trim-up`/`--trim-down` match paths by sample name, so they work with any builder.
- minigraph is part of the conda environment.

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
