# preprocessing.wdl inputs and outputs

This entrypoint supports both preprocessing starting states: HiFi demux mode
for one raw HiFi BAM, and cDNA demux-only mode for one or more already
HiFi-demultiplexed BAMs.

Both modes follow the same downstream preprocessing path: `skera split`, cDNA
demultiplexing with `lima`, and `isoseq refine`.

You can run this workflow as a standalone preprocessing entrypoint. The
[end-to-end HiFi BAM workflow](./kinnex_isoseq.md) also calls it before the
FLNC-level `secondary_analysis` stages.

## Inputs

| Type | Name | Description |
| ---- | ---- | ----------- |
| Array[HiFiSource] | hifi_sources | Input HiFi source BAMs. Omit `hifi_barcode` for HiFi demux mode, or provide one `hifi_barcode` per BAM for cDNA demux-only mode. |
| File | [biosample_csv](./biosample_csv.md) | Shared 3-column biosample CSV with exact header `HiFi Barcode,cDNA Barcode,Bio Sample`. |
| String | [ref_name](./reference_container.md#named-reference-selection) | Packaged reference to use when `reference_container` is omitted. Choices: `GRCh38_gencode49`. Default: `"GRCh38_gencode49"`. |
| String? | [reference_container](./reference_container.md#explicit-container-override) | Optional explicit immutable reference-container URI ending in a 64-character lowercase `@sha256:` digest. It takes precedence over `ref_name`. |
| ReferenceOverrides | [reference_overrides](./reference_container.md#typed-overrides) | Optional typed per-file overrides. Default: empty. |
| String? | segmentation_adapter_set | Optional segmentation adapter set used to select packaged Skera adapters when `reference_overrides.segmentation_adapters` is absent: `8-fold`, `12-fold`, or `16-fold`. Omission selects `8-fold`; supplying both inputs is an error. |
| String? | isoseq_primers_set | Optional Iso-Seq primer set used to select packaged indexed primers when `reference_overrides.indexed_primers` is absent: `IsoSeq-v2` or `IsoSeq96`. Omission selects `IsoSeq-v2`; supplying both inputs is an error. |
| Boolean | isoseq_require_polya | Pass `--require-polya` to downstream `isoseq refine`. Default: `true`. |
| String | backend | Execution backend. Default: `"HPC"`; only HPC is currently supported. |
| Int | max_retries | Maximum retries for failed task attempts. Default: `2`. |
| Int | add_memory_mb | Extra memory in MB added to each task request. Default: `0`. |
| String? | container_registry | Optional PacBio registry for registry-relative task images and named reference containers. If omitted, `"quay.io/pacbio"` is used; explicit full image references are not rewritten. |

The [reference-container specification](./reference_container.md) documents
the resources. The `biosample_csv` is run specific and must be prepared for
each run.

`HiFiSource` has the following fields:

| Type | Name | Description |
| ---- | ---- | ----------- |
| File | hifi_bam | Source HiFi BAM. |
| String? | hifi_barcode | HiFi barcode pair for an already HiFi-demuxed BAM. Omit for HiFi demux mode. |

The workflow infers the preprocessing mode from `hifi_sources`:

- HiFi demux mode: every `hifi_sources` entry omits `hifi_barcode`; exactly one
  source BAM is required.
- cDNA demux-only mode: every `hifi_sources` entry provides `hifi_barcode`; one
  or more already HiFi-demuxed source BAMs are allowed.

A single workflow run cannot mix these two starting states.

HiFi demux mode example:

```json
{
  "preprocessing.hifi_sources": [
    {
      "hifi_bam": "raw_movie.bam"
    }
  ]
}
```

cDNA demux-only mode example:

```json
{
  "preprocessing.hifi_sources": [
    {
      "hifi_bam": "movie.bcM0001.bam",
      "hifi_barcode": "bcM0001--bcM0001"
    },
    {
      "hifi_bam": "movie.bcM0002.bam",
      "hifi_barcode": "bcM0002--bcM0002"
    }
  ]
}
```

### `biosample_csv` sample sheet

`biosample_csv` is a shared workflow-level sample sheet and is always required.
It maps each outer SMRTbell barcode and cDNA Barcode combination to the
biological sample name that should be written to final FLNC BAM `SM` tags. See
[the biosample CSV specification](./biosample_csv.md) for the CSV schema,
examples, and mode-specific barcode rules.

The workflow derives lima-compatible 2-column CSVs from the 3-column sheet. For
both source states, it filters the sheet to one outer barcode and maps
`cDNA Barcode -> Bio Sample`. Downstream cDNA `lima` always receives this
derived per-outer CSV and overwrites `SM` tags with the CSV sample names.

### Validation

The workflow validates its inputs before upstream HiFi demux, `skera`, or cDNA
`lima` starts:

- source BAM basenames must be unique so source-dataset reporting identities
  remain unambiguous
- HiFi demux mode requires exactly one `hifi_sources` entry
- cDNA demux-only mode requires every `hifi_sources` entry to include
  `hifi_barcode`
- every `hifi_barcode` value must exist in `biosample_csv`
- in both modes, every source BAM must contain exactly one valid `@RG PU` movie
  name, and all source BAMs must share that movie name
- all cDNA Barcode pairs in `biosample_csv` are checked against the
  `indexed_primers`
- in HiFi demux mode, all outer barcode pairs in `biosample_csv` are
  checked against the `hifi_demux_barcodes` and must be symmetric
- in cDNA demux-only mode, all `hifi_barcode` values are checked against
  the `hifi_demux_barcodes`
- in cDNA demux-only mode, a recognized HiFi barcode in a source BAM basename
  that conflicts with the supplied `hifi_barcode` produces a warning
- `Bio Sample` names must be 40 characters or fewer and may contain only
  alphanumeric characters, underscores, and hyphens

The validated `PU` movie name, not the source BAM basename, supplies the output
prefix. Lima's split-named barcode suffixes and Refine's naming step produce
`<movie>.<hifi-barcode>.<cdna-barcode>.flnc.bam` names.

When preprocessing outputs are passed to `secondary_analysis`, FLNC grouping
derives sanitized filename prefixes from the final BAM `SM` values and fails if
two different sample names collide after sanitization.

HiFi Barcode names are exact lima barcode-pair names. The workflow does not
infer or normalize names such as `bcM0001` into `bcM0001--bcM0001`.

### HiFi demux and cDNA demux-only modes

- HiFi demux mode uses the `hifi_demux_barcodes`.
  The workflow runs `lima` with the symmetric-adapter preset, then scatters
  over each demuxed HiFi BAM and extracts the outer barcode from the filename
  suffix produced by `lima --split-named`. The `HiFi Barcode` values in
  `biosample_csv` therefore need to match the outer-barcode names coming from
  the barcode FASTA and must be symmetric pairs such as `bcM0001--bcM0001`.
- cDNA demux-only mode skips the upstream lima step entirely.
  Each `hifi_sources` entry must contain one `hifi_barcode` value. Each value
  must equal one of the outer `HiFi Barcode` values present in `biosample_csv`
  and selects the slice of the biosample CSV that applies to the corresponding
  already-demuxed input BAM. Multiple source BAMs may use the same
  `hifi_barcode` value.

### Dataset naming

Datasets normalized by the internal `hifi_demux` step are named as:

- `<movie-from-PU>.<barcode-pair>`

Final FLNC BAMs are named as:

- `<movie-from-PU>.<hifi-barcode-pair>.<cdna-barcode-pair>.flnc.bam`

These names are stable if an input BAM is renamed because the movie component
comes from the BAM header rather than its basename.

For example, an input BAM with `PU:m21003_240927_230709` and outer barcode
`bcM0001--bcM0001` becomes downstream dataset
`m21003_240927_230709.bcM0001--bcM0001`, regardless of the source BAM basename.
The same rule applies in HiFi-demux and cDNA-demux-only modes.

The workflow gathers the filter-summary report from every successful
`isoseq refine` call into one table-oriented PacBio report. Each row represents
one technical partition and retains the original `Bio Sample`, source BAM
basename, `HiFi Barcode`, and `cDNA Barcode`. Repeated `Bio Sample` names remain
unchanged and are disambiguated by those technical identity columns; report
metrics are not aggregated across repeated samples.

## Outputs

| Type | Name | Description |
| ---- | ---- | ----------- |
| String | workflow_name | Constant workflow identifier. |
| String | workflow_version | Workflow release version. |
| String? | reference_container_uri | Exact immutable reference-container URI used to resolve defaults, if any. |
| String | reference_mode | Effective source mode: `container`, `hybrid`, or `custom`. |
| String? | base_resource_bundle_version | Resource-bundle version from the reference container, if used. |
| Array[String] | source_dataset_names | Derived source dataset identifiers in scatter order. Each value is `basename(hifi_bam, ".bam")`. |
| Array[String] | dataset_names | Flattened downstream dataset names passed into the shared preprocessing stage. |
| Array[File] | hifi_demux_datasets | Flattened internal `hifi_demux` dataset XML outputs. |
| Array[File] | hifi_demux_bams | Flattened optional HiFi BAMs produced by internal `hifi_demux`. |
| Array[String] | flnc_names | Flattened FLNC output basenames across all normalized downstream datasets. |
| Array[File] | flnc_bams | Flattened FLNC BAMs across all normalized downstream datasets. |
| Array[File] | flnc_bam_pbis | Flattened PacBio BAM indexes for FLNC BAMs across all normalized downstream datasets. |
| File | refine_summary_report | Combined `isoseq refine` filter-summary report with one table row per successful technical partition. |

Tasks write tool logs in their execution directories for debugging. The
workflow does not expose those logs as outputs.
