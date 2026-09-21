# isocall.wdl inputs and outputs

The isocall stage runs `isocall profile` on each aligned FLNC BAM and prepares
known isoforms from the resolved compressed genome annotation with
`isocall prep-isoforms`. It then merges the generated profiles with
`isocall merge` and runs joint `isocall call`. You can include an
`isocall_extra_merged_profile` in the profile merge.

## Inputs

| Type | Name | Description |
| ---- | ---- | ----------- |
| Array[File] | aligned_bams | Aligned FLNC BAMs. Isocall uses each BAM's unique `@RG SM` value as the profile sample name when present; otherwise it falls back to its own filename-based sample naming. |
| Array[File] | aligned_bam_bais | BAM BAI files for `aligned_bams`, in the same order. These are localized next to each BAM before `isocall profile` so isocall can inspect the BAM header. |
| File? | isocall_extra_merged_profile | Optional merged profile added to the generated profiles before merging. It does not replace the generated profiles. |
| String | [ref_name](./reference_container.md#named-reference-selection) | Packaged reference to use when `reference_container` is omitted. Choices: `GRCh38_gencode49`. Default: `"GRCh38_gencode49"`. |
| String? | [reference_container](./reference_container.md#explicit-container-override) | Optional explicit immutable reference-container URI ending in a 64-character lowercase `@sha256:` digest. It takes precedence over `ref_name`. |
| ReferenceOverrides | [reference_overrides](./reference_container.md#typed-overrides) | Optional typed per-file overrides. Default: empty. |
| String | isocall_config_preset | `isocall call` configuration preset: `default` enables the standard filters, while `yolo` disables optional filters but retains read-support thresholds. Default: `"default"`. |
| File? | isocall_config_file | Optional custom TOML configuration file for `isocall call`. When supplied, this file takes precedence over `isocall_config_preset`. |
| String | output_prefix | Prefix for joint isocall outputs. Default: `"joint"`. |
| String | backend | Execution backend. Default: `"HPC"`; only HPC is currently supported. |
| Int | max_retries | Maximum retries for failed task attempts. Default: `2`. |
| Int | add_memory_mb | Extra memory in MB added to each task request. Default: `0`. |
| String? | container_registry | Optional PacBio registry for registry-relative task images and named reference containers. If omitted, `"quay.io/pacbio"` is used; explicit full image references are not rewritten. |

This standalone stage does not validate or sanitize BAM `SM` values. Developers
who need deterministic sample names should provide aligned BAMs with one unique
non-empty `@RG SM` value per BAM, or use `secondary_analysis` so FLNC grouping
handles sample-name validation and filename-prefix sanitization before
alignment.

Isocall configures all calling-method parameters through `--config`.
Custom configuration files can override sampling thresholds and filter
settings; omitted TOML sections retain their defaults.

Overriding `genome_fasta` also requires an explicit matching
`annotation_gtf_gz` override. The packaged annotation never falls back across a
custom-genome boundary.

## Outputs

| Type | Name | Description |
| ---- | ---- | ----------- |
| String | workflow_name | Constant workflow identifier. |
| String | workflow_version | Workflow release version. |
| String? | reference_container_uri | Exact immutable reference-container URI used to resolve defaults, if any. |
| String | reference_mode | Effective source mode: `container`, `hybrid`, or `custom`. |
| String? | base_resource_bundle_version | Resource-bundle version from the reference container, if used. |
| File | isocall_isoforms_gtf | Joint `isocall call` GTF output. |
| File | isocall_count_matrix | Joint `isocall call` per-sample supporting-read count matrix. |
| File | isocall_closest_known | Joint `isocall call` novel-to-known nearest-isoform table. |

Tasks write tool logs in their execution directories for debugging. The
workflow does not expose those logs as outputs.
