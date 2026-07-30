# secondary_analysis.wdl inputs and outputs

This entrypoint chains the currently implemented RNA stages into a single
workflow: it validates and groups FLNC BAMs by sample, aligns each input FLNC
BAM independently, merges same-sample aligned BAMs, runs isocall profiling and
joint isocall calling, and then runs isoform classification on the joint
isoform GTF and count matrix.

The workflow starts from FLNC BAMs. Use
[the end-to-end HiFi BAM workflow](./kinnex_isoseq.md) for an entrypoint
that includes preprocessing.

FLNC grouping, `pbmm2 align`, and same-sample aligned-BAM merging are
implemented in
[the FLNC alignment tasks](../workflows/secondary_analysis/alignment/flnc_alignment_tasks.wdl).
These tasks are orchestrated by this workflow and are not exposed as a
standalone public grouping or alignment entrypoint.

The workflow builds one `pbmm2` reference index with the `ISOSEQ` preset and
passes that index to every per-input alignment. Reference indexing is outside
the alignment scatter, so it is performed once per workflow run. Input BAMs
within the same sample group can therefore align in parallel. Multi-BAM sample
groups are merged with `pbsamoa` after alignment, while one-BAM groups pass
their `pbmm2` BAM and BAI through directly.

## Inputs

| Type | Name | Description |
| ---- | ---- | ----------- |
| Array[File] | flnc_bams | Input FLNC BAMs. Each BAM must contain `@RG` records, every `@RG` record must contain a non-empty `SM` value, and all `SM` values within one BAM must be identical. BAMs are associated by the same exact original `SM` before each BAM is aligned independently. |
| String? | [reference_container](./reference_container.md) | Optional immutable reference-container URI ending in a 64-character lowercase `@sha256:` digest. It supplies defaults for resources not overridden below. |
| ReferenceOverrides | [reference_overrides](./reference_container.md#typed-overrides) | Optional typed per-file overrides. Default: empty. |
| Boolean | pigeon_use_polya | Whether to pass the resolved `pigeon_poly_a` resource to `pigeon classify --poly-a` when that resource is available. Default: `true`. |
| Boolean | pigeon_use_cage_peak | Whether to prepare and pass the resolved `pigeon_cage_peak_bed` resource to `pigeon classify --cage-peak` when that resource is available. Default: `true`. |
| Boolean | pigeon_use_junction | Whether to prepare and pass the resolved `pigeon_junction_coverage` resource to `pigeon classify --coverage` when that resource is available. Default: `true`. |
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

The Pigeon-specific resources come from the resolved reference set and are used
only by the isoform-classification stage. The workflow validates that the
annotation GTF is already in Pigeon input order. It materializes CAGE and
junction resources and generates their Pigeon `.pgi` indices only when the
resource is available and its `pigeon_use_*` boolean is `true`. PolyA is
similarly passed only when available and enabled. The joint `isocall call`
count matrix is forwarded to `pigeon classify --flnc`.

Overriding `genome_fasta` also requires an explicit matching
`annotation_gtf_gz` override. It disables fallback to all packaged Pigeon
support resources; each support is then independently enabled by supplying its
override. See the [reference-resource specification](./reference_container.md#typed-overrides).

Task runtime blocks own CPU, base memory, retry, and Docker image selection.
Public inputs expose only the shared backend knobs above; `add_memory_mb` adds
to each task memory request.

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

Tool log files are generated inside task execution directories for debugging,
but are not exposed as workflow outputs.
