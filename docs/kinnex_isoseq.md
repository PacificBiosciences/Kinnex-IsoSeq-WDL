# kinnex_isoseq.wdl inputs and outputs

This workflow starts with HiFi source BAMs and runs preprocessing to generate
FLNC BAMs. It passes those BAMs to secondary analysis, which validates and
groups them by sample, aligns each BAM independently, and merges same-sample
alignments. Secondary analysis then runs isocall profiling and joint calling,
followed by isoform classification on the joint isoform GTF and count matrix.

Call [FLNC-level secondary analysis](./secondary_analysis.md) directly if you
already have FLNC BAMs.

Preprocessing produces read segmentation and Iso-Seq Refine reports covering
Skera metrics, primer assignment, refine results, and FLNC read lengths.
See the [preprocessing report overview](./preprocessing.md#reports) for contents
and interpretation. Keep each report with its companion PNGs and thumbnails.

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

Resources can come entirely from the reference container, entirely from typed
overrides, or from overrides combined with container defaults. The
`biosample_csv` and input BAMs remain run-specific filesystem inputs. See the
[reference-resource contract](./reference_container.md) for required fields,
override precedence, and validation rules.

A `genome_fasta` override also requires an explicit `annotation_gtf_gz`
override and disables fallback to the packaged Pigeon support files. This lets
you use non-GRCh38 genomes, including mouse, without silently mixing
genome-specific resources.

For `HiFiSource` fields, mode examples, source BAM validation rules, and
sample-sheet rules, see
[the preprocessing workflow interface](./preprocessing.md) and
[the biosample CSV specification](./biosample_csv.md).

Read segmentation metrics combine inputs by HiFi barcode. Refine report tables
retain source and barcode identities; their plots combine sources only when
the HiFi barcode, cDNA barcode, and biological sample all match.

The end-to-end wrapper also publishes the Transcript Mapping and Classification
report from secondary analysis. Filtered gene saturation always runs native
`pigeon report --exclude-singletons`. A transcript identity contributes to gene
detection only after it has more than one sampled supporting read. Known
reference transcript identities pool support, while novel isoform identities
remain separate. All sampled reads remain on the x-axis. This fixed saturation
policy does not change classification/filter outputs, native summary metrics,
or category/length plots.

The report has per-sample native summary metrics, unfiltered and filtered joint
classification tables, their category/length plots, and pooled filtered
all-gene/known-gene saturation plots. These make up three tables and three
two-plot groups. Keep `pigeon_summary_report` (`isoseq_mapping.report.json`) with
all twelve `pigeon_summary_report_resources` PNGs (six full images and six
thumbnails), including empty-state plots. The raw saturation TSV remains
internal. See the [classification report overview](./isoform_classification.md#transcript-mapping-and-classification-report)
for report contents and interpretation.

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
| File | segmentation_summary_csv | `segmentation.summary.csv`, with one weighted Skera metric row per distinct HiFi barcode in first-seen normalized-input order. |
| File | segmentation_summary_report | PacBio Read Segmentation report `segmentation.summary.report.json`, containing the same barcode table, all-input metrics, and four plots per distinct HiFi barcode with S-reads. |
| Array[File] | segmentation_summary_report_resources | Per-barcode Read Segmentation PNG images and `_thumb.png` thumbnails referenced by `segmentation_summary_report`; empty when the total S-read count is zero. |
| File | refine_summary_report | Iso-Seq Refine summary metrics, primer data, and read-length plots in `isoseq_refine.report.json`. |
| Array[File] | refine_summary_report_resources | PNG images and thumbnails referenced by the Iso-Seq Refine report; empty when no plots are produced. |
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
| File | pigeon_summary_report | Transcript Mapping and Classification report, `isoseq_mapping.report.json`, from the secondary-analysis classification stage. |
| Array[File] | pigeon_summary_report_resources | Twelve companion PNGs referenced by the report: six full-size plot images and six `_thumb.png` thumbnails. |

Each task writes tool logs to its execution directory for debugging. The
workflow does not expose those logs as outputs.
