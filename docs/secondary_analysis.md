# secondary_analysis.wdl inputs and outputs

This workflow validates FLNC BAMs and groups them by sample. It aligns each
input BAM independently, merges same-sample alignments, then runs isocall
profiling and joint calling. Isoform classification uses the resulting joint
isoform GTF and count matrix.

For an entrypoint that also includes preprocessing, use
[the end-to-end HiFi BAM workflow](./kinnex_isoseq.md).

The [FLNC alignment tasks](../workflows/secondary_analysis/alignment/flnc_alignment_tasks.wdl)
handle grouping, `pbmm2 align`, and merging aligned BAMs from the same sample.
This workflow coordinates those tasks but does not expose a standalone public
grouping or alignment entrypoint.

The workflow builds one `pbmm2` reference index with the `ISOSEQ` preset outside
the alignment scatter. Every per-input alignment uses that index, and BAMs in
the same sample group can align in parallel. After alignment, `pbsamoa` merges
multi-BAM sample groups. One-BAM groups pass their `pbmm2` BAM and BAI through
directly.

## Inputs

| Type | Name | Description |
| ---- | ---- | ----------- |
| Array[File] | flnc_bams | Input FLNC BAMs. Each BAM must contain `@RG` records, every `@RG` record must contain a non-empty `SM` value no longer than 100 characters, and all `SM` values within one BAM must be identical. BAMs are associated by the same exact original `SM` before each BAM is aligned independently. |
| String | [ref_name](./reference_container.md#named-reference-selection) | Packaged reference to use when `reference_container` is omitted. Choices: `GRCh38_gencode49`. Default: `"GRCh38_gencode49"`. |
| String? | [reference_container](./reference_container.md#explicit-container-override) | Optional explicit immutable reference-container URI ending in a 64-character lowercase `@sha256:` digest. It takes precedence over `ref_name`. |
| ReferenceOverrides | [reference_overrides](./reference_container.md#typed-overrides) | Optional typed per-file overrides. Default: empty. |
| Boolean | pigeon_use_polya | Whether to pass the resolved `pigeon_poly_a` resource to `pigeon classify --poly-a` when that resource is available. Default: `true`. |
| Boolean | pigeon_use_cage_peak | Whether to prepare and pass the resolved `pigeon_cage_peak_bed` resource to `pigeon classify --cage-peak` when that resource is available. Default: `true`. |
| Boolean | pigeon_use_junction | Whether to prepare and pass the resolved `pigeon_junction_coverage` resource to `pigeon classify --coverage` when that resource is available. Default: `true`. |
| File? | isocall_extra_merged_profile | Optional merged profile added to the generated isocall profiles before joint calling. |
| String | isocall_config_preset | `isocall call` configuration preset: `default` enables the standard filters, while `yolo` disables optional filters but retains read-support thresholds. Default: `"default"`. |
| Int | pigeon_min_ref_length | `pigeon classify --min-ref-length` parameter. Default: `100`. |
| String | output_prefix | Shared prefix for the joint-calling and classification outputs. Default: `"joint"`. |
| String | backend | Execution backend. Default: `"HPC"`; only HPC is currently supported. |
| Int | max_retries | Maximum retries for failed task attempts. Default: `2`. |
| Int | add_memory_mb | Extra memory in MB added to each task request. Default: `0`. |
| String? | container_registry | Optional PacBio registry for registry-relative task images and named reference containers. If omitted, `"quay.io/pacbio"` is used; explicit full image references are not rewritten. |

For external FLNC BAM inputs, `secondary_analysis` has less restrictive
sample-name rules than preprocessing. Original `SM` values may contain punctuation and spaces,
but cannot be empty, exceed 100 characters, or contain tab, carriage-return, or
newline characters. The workflow preserves the exact original `SM` as the
sample name and in read-group metadata. For filenames, it derives a sanitized
sample prefix by replacing every run of non-alphanumeric characters with `_`,
removing leading and trailing `_`, and truncating the result to 40 characters.
The workflow fails if the prefix is empty or if two different original `SM`
values produce the same sanitized prefix.

Only the isoform-classification stage uses the Pigeon resources from the
resolved reference set. The workflow checks that the annotation GTF is already
in Pigeon input order. It prepares CAGE and junction resources and generates
their Pigeon `.pgi` indices only when the resource is available and its
`pigeon_use_*` boolean is `true`. It passes polyA only when that resource is
available and enabled, and forwards the joint `isocall call` count matrix to
`pigeon classify --flnc`.

A `genome_fasta` override also requires an explicit matching
`annotation_gtf_gz` override and disables fallback to all packaged Pigeon
support resources. You can then enable each support resource independently by
supplying its override. See the [reference-resource specification](./reference_container.md#typed-overrides).

Task runtime blocks specify CPU, base memory, retry behavior, and Docker images.
Public inputs expose only the shared backend settings above; `add_memory_mb`
adds to each task's memory request.

## Outputs

| Type | Name | Description |
| ---- | ---- | ----------- |
| String | workflow_name | Constant workflow identifier. |
| String | workflow_version | Workflow release version. |
| String? | reference_container_uri | Exact immutable reference-container URI used to resolve defaults, if any. |
| String | reference_mode | Effective source mode: `container`, `hybrid`, or `custom`. |
| String? | base_resource_bundle_version | Resource-bundle version from the reference container, if used. |
| Array[String] | sample_names | Exact `SM` values from input FLNC BAM headers, in first-seen input order. |
| Array[String] | sample_prefixes | Filename-safe sanitized prefixes derived from `sample_names`, in sample order. These prefixes are used for final aligned BAM filenames. |
| Array[Int] | group_sizes | Number of input FLNC BAMs in each sample group. |
| Array[File] | aligned_bams | One coordinate-sorted aligned BAM per sample. |
| Array[File] | aligned_bam_bais | BAM BAI files for `aligned_bams`, in sample order. |
| File | isocall_isoforms_gtf | Joint `isocall call` GTF output. |
| File | isocall_count_matrix | Joint `isocall call` per-sample supporting-read count matrix. |
| File | isocall_closest_known | Joint `isocall call` novel-to-known nearest-isoform table. |
| File | pigeon_classification | Output from `pigeon classify`. |
| File | filtered_isoforms_gtf | Filtered lite isoform GTF from `pigeon filter`. |
| File | pigeon_filtered_classification | Filtered lite classification table from `pigeon filter`. |
| File | pigeon_summary_report | Transcript Mapping and Classification report, `isoseq_mapping.report.json`, with per-sample native metrics, joint classification tables and plots, and filtered pooled gene saturation. |
| Array[File] | pigeon_summary_report_resources | Twelve PNG companions referenced by the report: six full-size images and six `_thumb.png` thumbnails. |

Each task writes tool logs to its execution directory for debugging. The
workflow does not expose those logs as outputs.

## Transcript mapping and classification reporting

Secondary analysis publishes the report and companion resources generated by
the isoform-classification stage. The report has three tables and three
two-plot groups, labeled `(All Samples)`. Only the summary table is per sample.
It matches native classify/filter metrics by exact sample name and includes
samples with zero retained support and additional merged-profile samples.
Classification tables and plots count joint isoforms once, while saturation
pools supporting reads from the filtered classification.

Saturation always runs native `pigeon report --exclude-singletons`. A
transcript identity contributes to gene detection only after it has more than
one sampled supporting read. Known reference transcript identities pool
support; novel isoform identities remain separate. All sampled reads remain on
the x-axis. This fixed saturation policy does not change classification,
filtering, native summary metrics, or category/length plots. There is no
configurable input for the policy. The native saturation
TSV remains internal. Keep the JSON with all twelve PNG resources, including
empty-state plots and thumbnails.

See the [classification report overview](./isoform_classification.md#transcript-mapping-and-classification-report)
for report contents and interpretation.
