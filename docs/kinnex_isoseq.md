# kinnex_isoseq.wdl inputs and outputs

This entrypoint chains preprocessing with FLNC-level secondary analysis.

Starting from HiFi source BAMs, it generates FLNC BAMs in preprocessing and
passes them to secondary analysis. Secondary analysis validates and groups the
BAMs by sample, aligns each BAM independently, merges same-sample alignments,
runs isocall profiling and joint calling, and then runs isoform classification
on the joint isoform GTF and count matrix.

Call [FLNC-level secondary analysis](./secondary_analysis.md) directly if you
already have FLNC BAMs.

## Inputs

| Type | Name | Description |
| ---- | ---- | ----------- |
| Array[HiFiSource] | hifi_sources | Input HiFi source BAMs. Omit `hifi_barcode` for HiFi demux mode, or provide one `hifi_barcode` per BAM for cDNA demux-only mode. |
| File | [biosample_csv](./biosample_csv.md) | Shared 3-column biosample CSV with exact header `HiFi Barcode,cDNA Barcode,Bio Sample`. |
| String | [ref_name](./reference_container.md#named-reference-selection) | Packaged reference to use when `reference_container` is omitted. Choices: `GRCh38_gencode49`. Default: `"GRCh38_gencode49"`. |
| String? | [reference_container](./reference_container.md#explicit-container-override) | Optional explicit immutable reference-container URI ending in a 64-character lowercase `@sha256:` digest. It takes precedence over `ref_name`. |
| ReferenceOverrides | [reference_overrides](./reference_container.md#typed-overrides) | Optional typed per-file overrides. Default: empty. |
| String? | segmentation_adapter_set | Optional segmentation adapter set used to select packaged Skera adapters when `reference_overrides.segmentation_adapters` is absent: `8-fold`, `12-fold`, or `16-fold`. Omission selects `8-fold`; supplying both inputs is an error. |
| String? | isoseq_primers_set | Optional Iso-Seq primer set used to select packaged indexed primers when `reference_overrides.indexed_primers` is absent: `IsoSeq-v2` or `IsoSeq96`. Omission selects `IsoSeq-v2`; supplying both inputs is an error. |
| String | isocall_config_preset | `isocall call` configuration preset: `default` enables the standard filters, while `yolo` disables optional filters but retains read-support thresholds. Default: `"default"`. |
| String | output_prefix | Shared prefix for the joint-calling and classification outputs. Default: `"joint"`. |
| String | backend | Execution backend. Default: `"HPC"`; only HPC is currently supported. |
| Int | max_retries | Maximum retries for failed task attempts. Default: `2`. |
| Int | add_memory_mb | Extra memory in MB added to each task request. Default: `0`. |
| String? | container_registry | Optional PacBio registry for registry-relative task images and named reference containers. If omitted, `"quay.io/pacbio"` is used; explicit full image references are not rewritten. |

Effective resources can come entirely from the reference container or typed
overrides, or from overrides layered on the container defaults. The
`biosample_csv` and input BAMs remain run-specific filesystem inputs. See the
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

Tasks write tool logs in their execution directories for debugging. The
workflow does not expose those logs as outputs.
