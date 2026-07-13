# isoform_classification.wdl inputs and outputs

This entrypoint implements the isoform classification stage. It prepares the
annotation, CAGE, and junction coverage resources for Pigeon, normalizes and
sorts the joint isoform GTF, runs `pigeon classify`, and then runs
`pigeon filter`.

This stage is separate from isocall. It consumes the joint isoform GTF and
matching isocall count matrix, and passes the count matrix to
`pigeon classify --flnc`. Pigeon classify support resources are selected through
`ref_map_file`; the workflow validates the annotation GTF order, CAGE file, and
junction coverage file, and generates the required Pigeon `.pgi` indices.
Boolean workflow flags control whether `--poly-a`, `--cage-peak`, and
`--coverage` are enabled for a given run.

## Inputs

| Type | Name | Description |
| ---- | ---- | ----------- |
| File | isoforms_gtf | Joint isoform GTF from the isocall stage. Gzipped input is accepted and normalized internally. |
| File | isocall_count_matrix | Matching `isocall call` supporting-read count matrix passed to `pigeon classify --flnc`. |
| File | [ref_map_file](./ref_map.md) | Reference-map TSV. Isoform classification requires `name`, `annotation_gtf_gz`, `genome_fasta`, `genome_fasta_index`, `pigeon_poly_a`, `pigeon_cage_peak_bed`, and `pigeon_junction_coverage`. Shared maps may also carry isocall keys. |
| Boolean | pigeon_use_polya | Whether to pass the `pigeon_poly_a` resource from `ref_map_file` to `pigeon classify --poly-a`. Default: `true`. |
| Boolean | pigeon_use_cage_peak | Whether to pass the prepared `pigeon_cage_peak_bed` resource to `pigeon classify --cage-peak`. Default: `true`. |
| Boolean | pigeon_use_junction | Whether to pass the prepared `pigeon_junction_coverage` resource to `pigeon classify --coverage`. Default: `true`. |
| Int | pigeon_min_ref_length | `pigeon classify --min-ref-length` parameter. Default: `100`. |
| String | output_prefix | Prefix for classification outputs. Default: `"joint"`. |
| String | backend | Backend where the workflow will be executed. Default: `"HPC"`; only HPC is currently supported. |
| Int | max_retries | Maximum retries for failed task attempts. Default: `2`. |
| Int | add_memory_mb | Add Task Memory (MB). Increasing this number allocates extra memory per task when submitting jobs to the compute backend. Default: `0`. |
| String? | container_registry | Optional PacBio registry for registry-relative task images. If omitted, `"quay.io/pacbio"` is used; full image references are not rewritten. |

Task runtime blocks own CPU, base memory, retry, and Docker image selection.
`add_memory_mb` adds to each task memory request. The `pigeon` tasks use a
pinned PacBio image from `"quay.io/pacbio"` by default.

## Outputs

| Type | Name | Description |
| ---- | ---- | ----------- |
| String | workflow_name | Constant workflow identifier. |
| String | workflow_version | Workflow release version. |
| String | reference_name | Reference name from `ref_map_file`. |
| File | pigeon_classification | Output from `pigeon classify`. |
| File | filtered_isoforms_gtf | Filtered lite isoform GTF from `pigeon filter`. |
| File | pigeon_filtered_classification | Filtered lite classification table from `pigeon filter`. |

Tool log files are generated inside task execution directories for debugging,
but are not exposed as workflow outputs.
