# kinnex_isoseq.wdl inputs and outputs

This entrypoint chains preprocessing with the FLNC-level secondary analysis
workflow.

It starts from HiFi source BAMs, runs preprocessing to generate FLNC BAMs,
passes those FLNC BAMs into secondary analysis, validates and groups them by
sample, aligns each FLNC BAM independently, merges same-sample aligned BAMs,
runs isocall profiling and joint isocall calling, and then runs isoform
classification on the joint isoform GTF and count matrix.

The [FLNC-level secondary analysis](./secondary_analysis.md) entrypoint remains
the FLNC-to-classification workflow for callers that already have FLNC BAMs.

## Inputs

| Type | Name | Description |
| ---- | ---- | ----------- |
| Array[HiFiSource] | hifi_sources | Input HiFi source BAMs. Omit `hifi_barcode` for HiFi demux mode, or provide one `hifi_barcode` per BAM for cDNA demux-only mode. |
| File | [biosample_csv](./biosample_csv.md) | Shared 3-column biosample CSV with exact header `HiFi Barcode,cDNA Barcode,Bio Sample`. |
| String? | [reference_container](./reference_container.md) | Optional immutable reference-container URI ending in a 64-character lowercase `@sha256:` digest. It supplies defaults for resources not overridden below. |
| ReferenceOverrides | [reference_overrides](./reference_container.md#typed-overrides) | Optional typed per-file overrides. Default: empty. |
| String? | kinnex_primers_set | Optional Kinnex primer set used to select packaged Skera adapters when `reference_overrides.skera_adapters` is absent: `8fold`, `12fold`, or `16fold`. Omission selects `8fold`; supplying both inputs is an error. |
| Array[File]? | consensusreadset_xmls | Optional SMRT Link/internal input ConsensusReadSet XMLs for FLNC dataset XML generation. |
| String | output_prefix | Shared prefix for the joint-calling and classification outputs. Default: `"joint"`. |
| String | backend | Backend where the workflow will be executed. Default: `"HPC"`; only HPC is currently supported. |
| Int | max_retries | Maximum retries for failed task attempts. Default: `2`. |
| Int | add_memory_mb | Add Task Memory (MB). Increasing this number allocates extra memory per task when submitting jobs to the compute backend. Default: `0`. |
| String? | container_registry | Optional PacBio registry for registry-relative task images. If omitted, `"quay.io/pacbio"` is used; full image references are not rewritten. |

The effective resources can come entirely from the reference container,
entirely from typed overrides, or from overrides layered on the container
defaults. The `biosample_csv`, input BAMs, and optional dataset XMLs remain
run-specific filesystem inputs. See the
[reference-resource contract](./reference_container.md) for required fields,
override precedence, and validation rules.

Overriding `genome_fasta` also requires an explicit `annotation_gtf_gz`
override and disables fallback to the packaged Pigeon support files. This
supports non-GRCh38 genomes, including mouse, without silently mixing
genome-specific resources.

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
| String? | reference_container_uri | Exact immutable reference-container URI used to resolve defaults, if any. |
| String | reference_mode | Effective source mode: `container`, `hybrid`, or `custom`. |
| String? | base_resource_bundle_version | Resource-bundle version from the reference container, if used. |
| Array[String] | source_dataset_names | Derived source dataset identifiers in preprocessing scatter order. |
| Array[String] | preprocessing_dataset_names | Flattened downstream dataset names passed into the preprocessing stage. |
| Array[File] | hifi_demux_datasets | Flattened internal `hifi_demux` dataset XML outputs. |
| Array[File] | hifi_demux_bams | Flattened optional HiFi BAMs produced by internal `hifi_demux`. |
| Array[String] | flnc_names | Flattened FLNC output basenames across all normalized downstream datasets. |
| Array[File] | flnc_bams | Flattened FLNC BAMs across all normalized downstream datasets. |
| Array[File] | flnc_bam_pbis | Flattened PacBio BAM indexes for FLNC BAMs across all normalized downstream datasets. |
| File | refine_summary_report | Combined `isoseq refine` filter-summary report with one table row per successful technical partition, identified by Bio Sample, source dataset, HiFi barcode, and cDNA barcode. |
| File? | flnc_dataset_xml | Optional generated top-level FLNC ConsensusReadSet XML when `consensusreadset_xmls` is provided. |
| Array[File]? | flnc_child_dataset_xmls | Optional generated per-biosample child FLNC ConsensusReadSet XMLs. |
| Array[File]? | flnc_dataset_bams | Optional FLNC BAMs referenced by generated dataset XMLs. |
| Array[File]? | flnc_dataset_bam_pbis | Optional FLNC BAM indexes referenced by generated dataset XMLs. |
| Array[String] | sample_names | Exact `SM` values from input FLNC BAM headers, in first-seen input order. |
| Array[String] | sample_prefixes | Filename-safe sanitized prefixes derived from `sample_names`, in sample order. For preprocessing-produced FLNC BAMs, strict `Bio Sample` names usually make these prefixes match the sample names except that non-alphanumeric separators such as `_` and `-` are normalized to `_`. |
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
