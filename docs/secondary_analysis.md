# secondary_analysis.wdl inputs and outputs

This entrypoint chains the currently implemented RNA stages into a single
workflow: it groups and merges FLNC BAMs by sample, aligns one grouped FLNC BAM
per sample, runs isocall profiling and joint isocall calling, and then runs
isoform classification on the joint isoform GTF and count matrix.

The workflow starts from FLNC BAMs. Use
[the end-to-end HiFi BAM workflow](./kinnex_isoseq.md) for an entrypoint
that includes preprocessing.

FLNC grouping, same-sample merging, and `pbmm2 align` are implemented in
[the FLNC alignment tasks](../workflows/secondary_analysis/alignment/flnc_alignment_tasks.wdl).
These tasks are orchestrated by this workflow and are not exposed as a
standalone public grouping or alignment entrypoint.

## Inputs

| Type | Name | Description |
| ---- | ---- | ----------- |
| Array[File] | flnc_bams | Input FLNC BAMs. Each BAM must contain `@RG` records, every `@RG` record must contain a non-empty `SM` value, and all `SM` values within one BAM must be identical. BAMs with the same exact original `SM` are grouped before alignment. |
| File | [ref_map_file](./ref_map.md) | Reference-map TSV containing the union of keys required by FLNC alignment, isocall, and isoform classification. |
| Boolean | pigeon_use_polya | Whether to pass the `pigeon_poly_a` resource from `ref_map_file` to `pigeon classify --poly-a`. Default: `true`. |
| Boolean | pigeon_use_cage_peak | Whether to pass the prepared `pigeon_cage_peak_bed` resource to `pigeon classify --cage-peak`. Default: `true`. |
| Boolean | pigeon_use_junction | Whether to pass the prepared `pigeon_junction_coverage` resource to `pigeon classify --coverage`. Default: `true`. |
| File? | isocall_extra_merged_profile | Optional merged profile. If provided, it is merged together with generated isocall profiles before joint calling. |
| Float | isocall_min_read_fraction | `isocall call` parameter. Default: `0.99`. |
| Int | isocall_max_bundles_per_gene | `isocall call` parameter. Default: `10000`. |
| Int | isocall_min_reads_per_isoform | `isocall call` parameter. Default: `3`. |
| Int | pigeon_min_ref_length | `pigeon classify --min-ref-length` parameter. Default: `100`. |
| String | output_prefix | Shared prefix for the joint-calling and classification outputs. Default: `"joint"`. |
| String | backend | Backend where the workflow will be executed. Default: `"HPC"`; only HPC is currently supported. |
| Int | max_retries | Maximum retries for failed task attempts. Default: `2`. |
| Int | add_memory_mb | Add Task Memory (MB). Increasing this number allocates extra memory per task when submitting jobs to the compute backend. Default: `0`. |
| String? | container_registry | Optional PacBio registry for registry-relative task images. If omitted, `"quay.io/pacbio"` is used; full image references are not rewritten. |

For external FLNC BAM inputs, `secondary_analysis` is more permissive than
preprocessing: original `SM` values may contain punctuation and spaces, but must
not be empty or contain tab, carriage-return, or newline characters. The exact
original `SM` is preserved as the sample name and in read-group metadata. For
filenames, the workflow derives a sanitized sample prefix by replacing every run
of non-alphanumeric characters with `_`, removing leading/trailing `_`, and
truncating the result to 40 characters. If this produces an empty prefix, or if
two different original `SM` values produce the same sanitized prefix, the
workflow fails.

The Pigeon-specific resources live in [`ref_map_file`](./ref_map.md) and are
used only by the isoform-classification stage. The workflow validates that the
annotation GTF is already in Pigeon input order, materializes CAGE and junction
coverage resources, and generates Pigeon `.pgi` indices internally before
classification. The `pigeon_use_*` booleans decide whether `pigeon classify`
receives the corresponding arguments for a given run. The joint `isocall call`
count matrix is forwarded to `pigeon classify --flnc`.

Task runtime blocks own CPU, base memory, retry, and Docker image selection.
Public inputs expose only the shared backend knobs above; `add_memory_mb` adds
to each task memory request.

## Outputs

| Type | Name | Description |
| ---- | ---- | ----------- |
| String | workflow_name | Constant workflow identifier. |
| String | workflow_version | Workflow release version. |
| String | reference_name | Reference name from `ref_map_file`. |
| Array[String] | sample_names | Exact `SM` values from grouped FLNC BAM headers, in first-seen input order. |
| Array[String] | sample_prefixes | Filename-safe sanitized prefixes derived from `sample_names`, in sample order. These prefixes are used for grouped FLNC and aligned BAM filenames. |
| Array[Int] | group_sizes | Number of input FLNC BAMs in each sample group. |
| Array[File] | grouped_flnc_bams | FLNC BAMs after per-sample grouping, named `<sample_prefix>.flnc.bam`. Every group, including one-BAM groups, is materialized with `samtools merge`. |
| Array[File] | aligned_bams | Sorted `pbmm2` alignment outputs. |
| Array[File] | aligned_bam_bais | BAM BAI files for the aligned BAMs. |
| File | isocall_isoforms_gtf | Joint `isocall call` GTF output. |
| File | isocall_count_matrix | Joint `isocall call` per-sample supporting-read count matrix. |
| File | isocall_closest_known | Joint `isocall call` novel-to-known nearest-isoform table. |
| File | pigeon_classification | Output from `pigeon classify`. |
| File | filtered_isoforms_gtf | Filtered lite isoform GTF from `pigeon filter`. |
| File | pigeon_filtered_classification | Filtered lite classification table from `pigeon filter`. |

Tool log files are generated inside task execution directories for debugging,
but are not exposed as workflow outputs.
