# isoform_classification.wdl inputs and outputs

The isoform classification stage prepares the annotation for Pigeon, normalizes
and sorts the joint isoform GTF, then runs `pigeon classify` and `pigeon filter`.

This stage is separate from isocall. It takes the joint isoform GTF and
matching isocall count matrix, and passes the count matrix to
`pigeon classify --flnc`. Support resources come from the resolved reference
set. The workflow validates the annotation GTF order. When a CAGE or junction
coverage resource is available and enabled, it prepares the file and generates
its Pigeon `.pgi` index. Boolean workflow flags control whether Pigeon receives
available resources through `--poly-a`, `--cage-peak`, and `--coverage`.

## Inputs

| Type | Name | Description |
| ---- | ---- | ----------- |
| File | isoforms_gtf | Joint isoform GTF from the isocall stage. Gzipped input is accepted and normalized internally. |
| File | isocall_count_matrix | Matching `isocall call` supporting-read count matrix passed to `pigeon classify --flnc`. |
| String | [ref_name](./reference_container.md#named-reference-selection) | Packaged reference to use when `reference_container` is omitted. Choices: `GRCh38_gencode49`. Default: `"GRCh38_gencode49"`. |
| String? | [reference_container](./reference_container.md#explicit-container-override) | Optional explicit immutable reference-container URI ending in a 64-character lowercase `@sha256:` digest. It takes precedence over `ref_name`. |
| ReferenceOverrides | [reference_overrides](./reference_container.md#typed-overrides) | Optional typed per-file overrides. Default: empty. |
| Boolean | pigeon_use_polya | Whether to pass the resolved `pigeon_poly_a` resource to `pigeon classify --poly-a` when that resource is available. Default: `true`. |
| Boolean | pigeon_use_cage_peak | Whether to prepare and pass the resolved `pigeon_cage_peak_bed` resource to `pigeon classify --cage-peak` when that resource is available. Default: `true`. |
| Boolean | pigeon_use_junction | Whether to prepare and pass the resolved `pigeon_junction_coverage` resource to `pigeon classify --coverage` when that resource is available. Default: `true`. |
| Int | pigeon_min_ref_length | `pigeon classify --min-ref-length` parameter. Default: `100`. |
| String | output_prefix | Prefix for classification outputs. Default: `"joint"`. |
| String | backend | Execution backend. Default: `"HPC"`; only HPC is currently supported. |
| Int | max_retries | Maximum retries for failed task attempts. Default: `2`. |
| Int | add_memory_mb | Extra memory in MB added to each task request. Default: `0`. |
| String? | container_registry | Optional PacBio registry for registry-relative task images and named reference containers. If omitted, `"quay.io/pacbio"` is used; explicit full image references are not rewritten. |

Task runtime blocks specify CPU, base memory, retry behavior, and Docker images.
`add_memory_mb` adds to each task's memory request. The `pigeon` tasks use a
pinned PacBio image from `"quay.io/pacbio"` by default.

A `genome_fasta` override also requires an explicit matching
`annotation_gtf_gz` override and disables fallback to the packaged Pigeon
support files. This prevents accidentally mixing the GRCh38 bundle with a
custom genome. You can still run `pigeon classify` and `pigeon filter` without
supplemental support files. See the
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
| File | pigeon_summary_report | Transcript Mapping and Classification report, `isoseq_mapping.report.json`, with per-sample summary metrics, joint classification tables and plots, and filtered gene saturation. |
| Array[File] | pigeon_summary_report_resources | Companion PNGs referenced by the report: six full-size plot images and six `_thumb.png` thumbnails. |

Each task writes tool logs to its execution directory for debugging. The
workflow does not expose those logs as outputs.

## Transcript mapping and classification report

The Transcript Mapping and Classification report (`isoseq_mapping.report.json`)
is available from standalone classification, secondary analysis, and the
end-to-end workflow. It contains:

- **Per-sample summary metrics:** gene and isoform counts before and after
  filtering, filtered known-gene and known-isoform counts, and filtered novel
  isoforms with expression greater than 1 TPM.
- **Joint classification tables and plots:** isoform counts by structural
  category, CAGE and polyA support, category percentages, and isoform-length
  distributions, before and after filtering.
- **Gene-saturation plots:** all-gene and known-gene detection as supporting
  reads are subsampled from the filtered classification, pooled across samples.

The summary has one row per biological sample, including samples with no
retained support. Joint classification statistics count each isoform once,
regardless of how many samples support it. Known-isoform summary metrics
deduplicate reference transcript identities, so their totals can differ from
the joint category counts.

CAGE and polyA support are assessed only when the corresponding resources are
available and enabled. Missing support resources do not mean biological
support is absent.

Gene saturation excludes singletons: a transcript identity contributes to gene
detection only after it has more than one sampled supporting read. This fixed
policy applies only to saturation; it does not change classification,
filtering, or summary metrics.

Keep `pigeon_summary_report` with all twelve `pigeon_summary_report_resources`
PNGs (six plots and six thumbnails). Empty classifications still produce the
report and labeled empty-state plots.
