# isoform_classification.wdl inputs and outputs

This workflow implements the isoform classification stage. It prepares the
annotation and any enabled, available CAGE and junction coverage resources for
Pigeon, normalizes and sorts the joint isoform GTF, runs `pigeon classify`, and
then runs `pigeon filter`.

This stage is separate from isocall. It consumes the joint isoform GTF and
matching isocall count matrix, and passes the count matrix to
`pigeon classify --flnc`. Pigeon classify support resources come from the
resolved reference set. The workflow validates the annotation GTF order and
generates Pigeon `.pgi` indices for CAGE and junction files that are both
available and enabled. Boolean workflow flags control whether available
resources are passed with `--poly-a`, `--cage-peak`, and `--coverage`.

## Inputs

| Type | Name | Description |
| ---- | ---- | ----------- |
| File | isoforms_gtf | Joint isoform GTF from the isocall stage. Gzipped input is accepted and normalized internally. |
| File | isocall_count_matrix | Matching `isocall call` supporting-read count matrix passed to `pigeon classify --flnc`. |
| String? | [reference_container](./reference_container.md) | Optional immutable reference-container URI ending in a 64-character lowercase `@sha256:` digest. It supplies defaults for resources not overridden below. |
| ReferenceOverrides | [reference_overrides](./reference_container.md#typed-overrides) | Optional typed per-file overrides. Default: empty. |
| Boolean | pigeon_use_polya | Whether to pass the resolved `pigeon_poly_a` resource to `pigeon classify --poly-a` when that resource is available. Default: `true`. |
| Boolean | pigeon_use_cage_peak | Whether to prepare and pass the resolved `pigeon_cage_peak_bed` resource to `pigeon classify --cage-peak` when that resource is available. Default: `true`. |
| Boolean | pigeon_use_junction | Whether to prepare and pass the resolved `pigeon_junction_coverage` resource to `pigeon classify --coverage` when that resource is available. Default: `true`. |
| Int | pigeon_min_ref_length | `pigeon classify --min-ref-length` parameter. Default: `100`. |
| String | output_prefix | Prefix for classification outputs. Default: `"joint"`. |
| String | backend | Backend where the workflow will be executed. Default: `"HPC"`; only HPC is currently supported. |
| Int | max_retries | Maximum retries for failed task attempts. Default: `2`. |
| Int | add_memory_mb | Add Task Memory (MB). Increasing this number allocates extra memory per task when submitting jobs to the compute backend. Default: `0`. |
| String? | container_registry | Optional PacBio registry for registry-relative task images. If omitted, `"quay.io/pacbio"` is used; full image references are not rewritten. |

Task runtime blocks own CPU, base memory, retry, and Docker image selection.
`add_memory_mb` adds to each task memory request. The `pigeon` tasks use a
pinned PacBio image from `"quay.io/pacbio"` by default.

Overriding `genome_fasta` also requires an explicit matching
`annotation_gtf_gz` override and disables fallback to the packaged Pigeon
supports. This prevents accidental mixing of the GRCh38 bundle with mouse or
another custom genome. `pigeon classify` and `pigeon filter` still run when no
supplemental supports are provided; see the
[reference-resource contract](./reference_container.md#typed-overrides) for
the scientific implications.

## Outputs

| Type | Name | Description |
| ---- | ---- | ----------- |
| String | workflow_name | Constant workflow identifier. |
| String | workflow_version | Workflow release version. |
| String? | reference_container_uri | Exact immutable reference-container URI used to resolve defaults, if any. |
| String | reference_mode | Effective source mode: `container`, `hybrid`, or `custom`. |
| String? | base_resource_bundle_version | Resource-bundle version from the reference container, if used. |
| File | pigeon_classification | Output from `pigeon classify`. |
| File | filtered_isoforms_gtf | Filtered lite isoform GTF from `pigeon filter`. |
| File | pigeon_filtered_classification | Filtered lite classification table from `pigeon filter`. |

Tool log files are generated inside task execution directories for debugging,
but are not exposed as workflow outputs.
