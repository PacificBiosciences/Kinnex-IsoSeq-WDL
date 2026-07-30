# preprocessing.wdl inputs and outputs

This entrypoint implements the preprocessing workflow for both supported
starting states: either HiFi demux mode for one raw HiFi BAM, or cDNA
demux-only mode for one or more already HiFi-demuxed BAMs.

After this normalization step, all datasets follow the same downstream
preprocessing path: `skera split`, cDNA demultiplexing with `lima`, and
`isoseq refine`.

This workflow remains available as a standalone preprocessing entrypoint. It is
also called by [the end-to-end HiFi BAM workflow](./kinnex_isoseq.md) before
the FLNC-level `secondary_analysis` stages.

## Inputs

| Type | Name | Description |
| ---- | ---- | ----------- |
| Array[HiFiSource] | hifi_sources | Input HiFi source BAMs. Omit `hifi_barcode` for HiFi demux mode, or provide one `hifi_barcode` per BAM for cDNA demux-only mode. |
| File | [biosample_csv](./biosample_csv.md) | Shared 3-column biosample CSV with exact header `HiFi Barcode,cDNA Barcode,Bio Sample`. |
| String? | [reference_container](./reference_container.md) | Optional immutable reference-container URI ending in a 64-character lowercase `@sha256:` digest. It supplies defaults for resources not overridden below. |
| ReferenceOverrides | [reference_overrides](./reference_container.md#typed-overrides) | Optional typed per-file overrides. Default: empty. |
| String? | kinnex_primers_set | Optional Kinnex primer set used to select packaged Skera adapters when `reference_overrides.skera_adapters` is absent: `8fold`, `12fold`, or `16fold`. Omission selects `8fold`; supplying both inputs is an error. |
| Array[File]? | consensusreadset_xmls | Optional SMRT Link/internal input ConsensusReadSet XMLs for FLNC dataset XML generation. |
| Boolean | isoseq_require_polya | Pass `--require-polya` to downstream `isoseq refine`. Default: `true`. |
| String | backend | Backend where the workflow will be executed. Default: `"HPC"`; only HPC is currently supported. |
| Int | max_retries | Maximum retries for failed task attempts. Default: `2`. |
| Int | add_memory_mb | Add Task Memory (MB). Increasing this number allocates extra memory per task when submitting jobs to the compute backend. Default: `0`. |
| String? | container_registry | Optional PacBio registry for registry-relative task images. If omitted, `"quay.io/pacbio"` is used; full image references are not rewritten. |

The resources are documented in the
[reference-container contract](./reference_container.md). The `biosample_csv`
is run-specific and must be prepared separately for each run.

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

The workflow internally derives lima-compatible 2-column CSVs from the
3-column sheet. For both source states, the cDNA step uses a 2-column CSV
filtered to the rows belonging to a single outer barcode, mapping
`cDNA Barcode -> Bio Sample`. The downstream cDNA `lima` step always receives this
derived per-outer CSV and overwrites `SM` tags with the CSV sample names.

### Validation

The workflow validates the simplified input contract before upstream HiFi demux,
`skera`, or cDNA `lima` starts:

- source BAM basenames must be unique because they become output prefixes
- HiFi demux mode requires exactly one `hifi_sources` entry
- cDNA demux-only mode requires every `hifi_sources` entry to include
  `hifi_barcode`
- every `hifi_barcode` value must exist in `biosample_csv`
- already HiFi-demuxed source BAMs must each contain exactly one distinct
  `@RG PU` value, and all BAMs must share the same `PU`
- all cDNA Barcode pairs in `biosample_csv` are checked against the
  `barcoded_primers`
- in HiFi demux mode, all outer barcode pairs in `biosample_csv` are
  checked against the `hifi_demux_barcodes` and must be symmetric
- in cDNA demux-only mode, all `hifi_barcode` values are checked against
  the `hifi_demux_barcodes`
- `Bio Sample` names must be 40 characters or fewer and may contain only
  alphanumeric characters, underscores, and hyphens

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

- `<basename(hifi_bam, ".bam")>.<barcode-pair>`

For example, a source BAM named `m21003_240927_230709.hifi_reads.bam` plus a
HiFi-demux BAM named with `bcM0001--bcM0001` becomes downstream dataset
`m21003_240927_230709.hifi_reads.bcM0001--bcM0001`.

Datasets in cDNA demux-only mode use `basename(hifi_bam, ".bam")` as the
downstream dataset name.

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
| File? | flnc_dataset_xml | Optional generated top-level FLNC ConsensusReadSet XML when `consensusreadset_xmls` is provided. |
| Array[File]? | flnc_child_dataset_xmls | Optional generated per-biosample child FLNC ConsensusReadSet XMLs. |
| Array[File]? | flnc_dataset_bams | Optional FLNC BAMs referenced by generated dataset XMLs. |
| Array[File]? | flnc_dataset_bam_pbis | Optional FLNC BAM indexes referenced by generated dataset XMLs. |

Tool log files are generated inside task execution directories for debugging,
but are not exposed as workflow outputs.
