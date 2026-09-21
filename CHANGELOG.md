# Changelog

This changelog records changes to the PacBio Kinnex Iso-Seq Pipeline.

## [Unreleased]

## [0.4.0]

### Breaking changes

- Renamed `ReferenceOverrides.barcoded_primers` to `indexed_primers` and
  `ReferenceOverrides.skera_adapters` to `segmentation_adapters`.
- Renamed `kinnex_primers_set` to `segmentation_adapter_set` and changed its
  values from `8fold`, `12fold`, and `16fold` to `8-fold`, `12-fold`, and `16-fold`.
- Replaced the Isocall threshold inputs on `secondary_analysis` and
  standalone `isocall` with `isocall_config_preset` (`default` or `yolo`).
  Standalone `isocall` additionally accepts a custom TOML configuration file,
  which takes precedence over the preset.

### New features

- Added named reference selection through `ref_name`, defaulting to
  `GRCh38_gencode49`. An explicit immutable `reference_container` takes
  precedence, and typed per-file overrides remain supported.
- Updated to Kinnex Iso-Seq resource bundle `0.2.0` support, including separate
  packaged Iso-Seq v2 and Iso-Seq 96 indexed-primer FASTAs. Preprocessing and
  end-to-end workflows can select them with `isoseq_primers_set`; omission
  selects `IsoSeq-v2`.

### Improvements

- Updated Isocall from 1.1.0 to [1.3.0](https://github.com/PacificBiosciences/isocall/releases/tag/1.3.0).
- Updated Slurm examples and development tooling to Sprocket 0.30.1, and added
  Jenkins validation and real-container task-test coverage.

### Fixes

- Fixed Cromwell container-path localization of optional File inputs used by
  Pigeon, Lima, and Isocall.

## [0.3.0]

### Breaking changes

- **New reference inputs:** `ref_map_file` and the separate reference-file
  inputs have been replaced by `reference_container` and
  `reference_overrides`. The provided reference container is specific to
  GRCh38 and was built from the
  [Kinnex Iso-Seq resource bundle on Zenodo](https://zenodo.org/records/10839617).
  You can use all packaged resources, replace only the files you need, or
  provide a fully custom set of resources.

### Improvements

- Added `kinnex_primers_set` with the values `8fold`, `12fold`, and `16fold`.
  When unset, it defaults to `8fold`.
- Parallelized FLNC BAM alignment before merging.

## [0.2.0]

- Updated [isocall from 1.0.0 to 1.1.0](https://github.com/PacificBiosciences/isocall/releases/tag/1.1.0).
- Fixed Cromwell execution of tasks that consume arrays of input files.

## [0.1.0]

- Initial release.
