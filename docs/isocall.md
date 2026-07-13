# isocall.wdl inputs and outputs

This entrypoint implements the isocall stage. It runs `isocall profile` on each
aligned FLNC BAM, prepares known isoforms from the compressed genome annotation
selected through the reference map with `isocall prep-isoforms`, merges the
generated profiles with `isocall merge`, and runs joint `isocall call`. An
optional `isocall_extra_merged_profile` can be included during profile merging.

## Inputs

| Type | Name | Description |
| ---- | ---- | ----------- |
| Array[File] | aligned_bams | Aligned FLNC BAMs. Isocall uses each BAM's unique `@RG SM` value as the profile sample name when present; otherwise it falls back to its own filename-based sample naming. |
| Array[File] | aligned_bam_bais | BAM BAI files for `aligned_bams`, in the same order. These are localized next to each BAM before `isocall profile` so isocall can inspect the BAM header. |
| File? | isocall_extra_merged_profile | Optional merged profile. If provided, it is added to the generated profiles and merged together with them. It does not replace the generated profiles. |
| File | [ref_map_file](./ref_map.md) | Reference-map TSV. Isocall requires `name`, `annotation_gtf_gz`, `genome_fasta`, and `genome_fasta_index`; shared maps may also include other stage keys. |
| Float | isocall_min_read_fraction | `isocall call` parameter. Default: `0.99`. |
| Int | isocall_max_bundles_per_gene | `isocall call` parameter. Default: `10000`. |
| Int | isocall_min_reads_per_isoform | `isocall call` parameter. Default: `3`. |
| String | output_prefix | Prefix for joint isocall outputs. Default: `"joint"`. |
| String | backend | Backend where the workflow will be executed. Default: `"HPC"`; only HPC is currently supported. |
| Int | max_retries | Maximum retries for failed task attempts. Default: `2`. |
| Int | add_memory_mb | Add Task Memory (MB). Increasing this number allocates extra memory per task when submitting jobs to the compute backend. Default: `0`. |
| String? | container_registry | Optional PacBio registry for registry-relative task images. If omitted, `"quay.io/pacbio"` is used; full image references are not rewritten. |

This standalone stage does not validate or sanitize BAM `SM` values. Callers that
need deterministic sample names should provide aligned BAMs with one unique
non-empty `@RG SM` value per BAM, or use `secondary_analysis` so FLNC grouping handles
sample-name validation and filename-prefix sanitization before alignment.

## Outputs

| Type | Name | Description |
| ---- | ---- | ----------- |
| String | workflow_name | Constant workflow identifier. |
| String | workflow_version | Workflow release version. |
| String | reference_name | Reference name from `ref_map_file`. |
| File | isocall_isoforms_gtf | Joint `isocall call` GTF output. |
| File | isocall_count_matrix | Joint `isocall call` per-sample supporting-read count matrix. |
| File | isocall_closest_known | Joint `isocall call` novel-to-known nearest-isoform table. |

Tool log files are generated inside task execution directories for debugging,
but are not exposed as workflow outputs.
