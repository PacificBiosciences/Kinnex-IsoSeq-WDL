# kinnex_isoseq.wdl inputs and outputs

This entrypoint chains preprocessing with the FLNC-level secondary analysis
workflow.
It starts from HiFi source BAMs, runs preprocessing to generate FLNC BAMs,
passes those FLNC BAMs into secondary analysis, groups, merges, and aligns FLNC
BAMs by sample, runs isocall profiling and joint isocall calling, and then runs
isoform classification on the joint isoform GTF and count matrix.

The [FLNC-level secondary analysis](./secondary_analysis.md) entrypoint remains
the FLNC-to-classification workflow for callers that already have FLNC BAMs.

## Inputs

| Type | Name | Description |
| ---- | ---- | ----------- |
| Array[HiFiSource] | hifi_sources | Input HiFi source BAMs. Omit `hifi_barcode` for HiFi demux mode, or provide one `hifi_barcode` per BAM for cDNA demux-only mode. |
| File | [biosample_csv](./biosample_csv.md) | Shared 3-column biosample CSV with exact header `HiFi Barcode,cDNA Barcode,Bio Sample`. |
| File | hifi_demux_barcodes | SMRTbell barcode FASTA. Required for every run. Used for upstream HiFi demux in HiFi demux mode and barcode validation in both modes. |
| File | skera_adapters | Adapter FASTA for `skera split`. |
| File | barcoded_primers | Primer FASTA shared by downstream cDNA `lima` and `isoseq refine`. |
| Array[File]? | consensusreadset_xmls | Optional SMRT Link/internal input ConsensusReadSet XMLs for FLNC dataset XML generation. |
| File | [ref_map_file](./ref_map.md) | Reference-map TSV containing the union of keys required by FLNC alignment, isocall, and isoform classification. |
| Boolean | isoseq_require_polya | Pass `--require-polya` to downstream `isoseq refine`. Default: `true`. |
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

The preprocessing barcode, adapter, and primer files are explicit workflow
inputs separate from `ref_map_file`, which describes genome, annotation, and
classification resources. See the [resource bundle layout](./resource_bundle.md)
for the default bundled paths. The `biosample_csv` is run-specific and is not
part of the shared resource bundle.

The `HiFiSource` fields, mode examples, source BAM validation rules, and
sample-sheet rules are described in
[the preprocessing workflow interface](./preprocessing.md) and
[the biosample CSV specification](./biosample_csv.md).

## Outputs

| Type | Name | Description |
| ---- | ---- | ----------- |
| String | workflow_name | Constant workflow identifier. |
| String | workflow_version | Workflow release version. |
| String | preprocessing_workflow_name | Child preprocessing workflow identifier. |
| String | secondary_analysis_workflow_name | Child secondary analysis workflow identifier. |
| String | reference_name | Reference name from `ref_map_file`. |
| Array[String] | source_dataset_names | Derived source dataset identifiers in preprocessing scatter order. |
| Array[String] | preprocessing_dataset_names | Flattened downstream dataset names passed into the preprocessing stage. |
| Array[File] | hifi_demux_datasets | Flattened internal `hifi_demux` dataset XML outputs. |
| Array[File] | hifi_demux_bams | Flattened optional HiFi BAMs produced by internal `hifi_demux`. |
| Array[String] | flnc_names | Flattened FLNC output basenames across all normalized downstream datasets. |
| Array[File] | flnc_bams | Flattened FLNC BAMs across all normalized downstream datasets. |
| Array[File] | flnc_bam_pbis | Flattened PacBio BAM indexes for FLNC BAMs across all normalized downstream datasets. |
| Array[File] | refine_summary_reports | Flattened `isoseq refine` filter-summary JSON reports across all normalized downstream datasets. |
| File? | flnc_dataset_xml | Optional generated top-level FLNC ConsensusReadSet XML when `consensusreadset_xmls` is provided. |
| Array[File]? | flnc_child_dataset_xmls | Optional generated per-biosample child FLNC ConsensusReadSet XMLs. |
| Array[File]? | flnc_dataset_bams | Optional packaged FLNC BAMs referenced by generated dataset XMLs. |
| Array[File]? | flnc_dataset_bam_pbis | Optional packaged FLNC BAM indexes referenced by generated dataset XMLs. |
| Array[String] | sample_names | Exact `SM` values from grouped FLNC BAM headers, in first-seen input order. |
| Array[String] | sample_prefixes | Filename-safe sanitized prefixes derived from `sample_names`, in sample order. For preprocessing-produced FLNC BAMs, strict `Bio Sample` names usually make these prefixes match the sample names except that non-alphanumeric separators such as `_` and `-` are normalized to `_`. |
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
