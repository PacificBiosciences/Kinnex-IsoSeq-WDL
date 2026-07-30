# Changelog

Any changes to PacBio Kinnex Iso-Seq Pipeline are noted below:

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

- You can now choose a Kinnex primer set set for `kinnex_primers_set` to `8fold`,
  `12fold`, or `16fold`. If you leave it unset, the workflow uses `8fold`.
- The workflow now aligns FLNC BAMs in parallel before merging them.
- The Slurm examples and documentation now target Sprocket 0.28.0.

## [0.2.0]

- Updated [isocall from 1.0.0 to 1.1.0](https://github.com/PacificBiosciences/isocall/releases/tag/1.1.0).
- Fixed Cromwell execution of tasks that consume arrays of input files.

## [0.1.0]

- Initial release.
