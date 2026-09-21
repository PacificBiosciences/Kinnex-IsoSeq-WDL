version 1.0

import "backend_configuration.wdl" as BackendConfiguration
import "preprocessing/preprocessing_core.wdl" as PreprocessingCore
import "reference_resources/reference_resources.wdl" as ReferenceResources
import "secondary_analysis/secondary_analysis_core.wdl" as SecondaryAnalysisCore

workflow kinnex_isoseq {
  meta {
    description: "PacBio Kinnex Iso-Seq end-to-end pipeline: run preprocessing from HiFi BAMs to FLNC BAMs, then FLNC grouping/alignment, joint isocall, and isoform classification."
    outputs: {
      workflow_name: {
        description: "Workflow name"
      },
      workflow_version: {
        description: "Workflow version"
      },
      preprocessing_workflow_name: {
        description: "Preprocessing workflow name"
      },
      secondary_analysis_workflow_name: {
        description: "Secondary analysis workflow name"
      },
      reference_container_uri: {
        description: "Immutable reference-container URI when container defaults were used"
      },
      reference_mode: {
        description: "Reference selection mode: container, hybrid, or custom"
      },
      base_resource_bundle_version: {
        description: "Base container resource-bundle version when container defaults were used"
      },
      source_dataset_names: {
        description: "Source dataset names"
      },
      preprocessing_dataset_names: {
        description: "Preprocessing dataset names"
      },
      hifi_demux_datasets: {
        description: "HiFi-demuxed ConsensusReadSet XMLs"
      },
      hifi_demux_bams: {
        description: "HiFi-demuxed BAMs"
      },
      flnc_names: {
        description: "FLNC names"
      },
      flnc_bams: {
        description: "FLNC BAMs"
      },
      flnc_bam_pbis: {
        description: "FLNC BAM PBIs"
      },
      refine_summary_report: {
        description: "Combined isoseq refine summary report"
      },
      sample_names: {
        description: "Sample names"
      },
      sample_prefixes: {
        description: "Sample prefixes"
      },
      group_sizes: {
        description: "FLNC BAM counts per sample group"
      },
      aligned_bams: {
        description: "Aligned FLNC BAMs"
      },
      aligned_bam_bais: {
        description: "Aligned FLNC BAM BAIs"
      },
      isocall_isoforms_gtf: {
        description: "Isocall isoforms GTF"
      },
      isocall_count_matrix: {
        description: "Isocall count matrix"
      },
      isocall_closest_known: {
        description: "Closest known isoforms"
      },
      pigeon_classification: {
        description: "Pigeon classification table"
      },
      filtered_isoforms_gtf: {
        description: "Filtered isoforms GTF"
      },
      pigeon_filtered_classification: {
        description: "Filtered pigeon classification table"
      }
    }
  }

  parameter_meta {
    hifi_sources: {
      description: "Input HiFi source BAMs. Omit hifi_barcode for HiFi demux mode, or provide one hifi_barcode per BAM for cDNA demux-only mode"
    }
    biosample_csv: {
      description: "Shared 3-column biosample CSV with header `HiFi Barcode,cDNA Barcode,Bio Sample`"
    }
    ref_name: {
      description: "Packaged reference to use when reference_container is omitted",
      choices: [
        "GRCh38_gencode49"
      ]
    }
    reference_container: {
      description: "Optional explicit immutable reference-container URI; takes precedence over ref_name and is not rewritten by container_registry"
    }
    reference_overrides: {
      description: "Typed optional overrides for reference files; a genome override requires an annotation override and disables packaged Pigeon support fallbacks"
    }
    segmentation_adapter_set: {
      description: "Optional segmentation adapter set used to select the packaged segmentation adapter FASTA when no custom override is supplied; omission selects 8-fold",
      choices: [
        "8-fold",
        "12-fold",
        "16-fold"
      ]
    }
    isoseq_primers_set: {
      description: "Optional Iso-Seq primer set used to select the packaged indexed-primer FASTA when no custom override is supplied; omission selects IsoSeq-v2",
      choices: [
        "IsoSeq-v2",
        "IsoSeq96"
      ]
    }
    isocall_config_preset: {
      description: "Isocall calling configuration preset",
      choices: [
        "default",
        "yolo"
      ]
    }
    output_prefix: {
      description: "Shared prefix for joint calling and classification outputs"
    }
    backend: {
      description: "Backend where the workflow will be executed. Only HPC is supported."
    }
    max_retries: {
      description: "Maximum retries for failed task attempts"
    }
    add_memory_mb: {
      name: "Add Task Memory (MB)",
      description: "Increasing this number allocates extra memory per task when submitting jobs to the compute backend."
    }
    nproc: {
      description: "Maximum CPU threads per task for development and backend throttling",
      hidden: true
    }
    container_registry: {
      description: "Optional PacBio registry for registry-relative task images and named reference containers. If omitted, quay.io/pacbio is used."
    }
  }

  input {
    Array[HiFiSource] hifi_sources
    File biosample_csv
    String ref_name = "GRCh38_gencode49"
    String? reference_container
    ReferenceOverrides reference_overrides = object {
    }
    String? segmentation_adapter_set
    String? isoseq_primers_set
    String isocall_config_preset = "default"
    String output_prefix = "joint"

    # Backend configuration
    String backend = "HPC"
    Int max_retries = 2
    Int add_memory_mb = 0
    Int nproc = 32
    String? container_registry
  }

  call BackendConfiguration.backend_configuration { input:
    backend = backend,
    max_retries = max_retries,
    add_memory_mb = add_memory_mb,
    nproc = nproc,
    container_registry = container_registry
  }

  RuntimeAttributes default_runtime_attributes = backend_configuration.runtime_attributes

  call ReferenceResources.resolve_reference_resources { input:
    resolution_profile = "kinnex_isoseq",
    ref_name = ref_name,
    reference_container = reference_container,
    reference_overrides = reference_overrides,
    segmentation_adapter_set = segmentation_adapter_set,
    isoseq_primers_set = isoseq_primers_set,
    runtime_attributes = default_runtime_attributes
  }

  ResolvedReferenceResources reference_resources = resolve_reference_resources.resources

  scatter (hifi_source in hifi_sources) {
    File hifi_bam = hifi_source.hifi_bam
    String? optional_hifi_barcode = hifi_source.hifi_barcode
  }

  Array[String] hifi_bam_barcodes = select_all(optional_hifi_barcode)
  Boolean needs_hifi_demux = length(hifi_bam_barcodes) == 0

  call PreprocessingCore.preprocessing_core as preprocessing_core { input:
    hifi_bams = hifi_bam,
    hifi_bam_barcodes = hifi_bam_barcodes,
    needs_hifi_demux = needs_hifi_demux,
    biosample_csv = biosample_csv,
    hifi_demux_barcodes = select_first([
      reference_resources.hifi_demux_barcodes
    ]),
    segmentation_adapters = select_first([
      reference_resources.segmentation_adapters
    ]),
    indexed_primers = select_first([
      reference_resources.indexed_primers
    ]),
    runtime_attributes = default_runtime_attributes
  }

  call SecondaryAnalysisCore.secondary_analysis_core as secondary_analysis_core { input:
    flnc_bams = preprocessing_core.flnc_bams,
    genome_fasta = select_first([
      reference_resources.genome_fasta
    ]),
    genome_fasta_index = select_first([
      reference_resources.genome_fasta_index
    ]),
    annotation_gtf_gz = select_first([
      reference_resources.annotation_gtf_gz
    ]),
    pigeon_poly_a = reference_resources.pigeon_poly_a,
    pigeon_cage_peak_bed = reference_resources.pigeon_cage_peak_bed,
    pigeon_junction_coverage = reference_resources.pigeon_junction_coverage,
    runtime_attributes = default_runtime_attributes,
    isocall_config_preset = isocall_config_preset,
    output_prefix = output_prefix
  }

  output {
    String workflow_name = "kinnex_isoseq"
    String workflow_version = "0.4.0"
    String preprocessing_workflow_name = "preprocessing"
    String secondary_analysis_workflow_name = "secondary_analysis"
    String? reference_container_uri = resolve_reference_resources.reference_container_uri
    String reference_mode = resolve_reference_resources.reference_mode
    String? base_resource_bundle_version = resolve_reference_resources.base_resource_bundle_version
    Array[String] source_dataset_names = preprocessing_core.source_dataset_names
    Array[String] preprocessing_dataset_names = preprocessing_core.dataset_names
    Array[File] hifi_demux_datasets = preprocessing_core.hifi_demux_datasets
    Array[File] hifi_demux_bams = preprocessing_core.hifi_demux_bams
    Array[String] flnc_names = preprocessing_core.flnc_names
    Array[File] flnc_bams = preprocessing_core.flnc_bams
    Array[File] flnc_bam_pbis = preprocessing_core.flnc_bam_pbis
    File refine_summary_report = preprocessing_core.refine_summary_report
    Array[String] sample_names = secondary_analysis_core.sample_names
    Array[String] sample_prefixes = secondary_analysis_core.sample_prefixes
    Array[Int] group_sizes = secondary_analysis_core.group_sizes
    Array[File] aligned_bams = secondary_analysis_core.aligned_bams
    Array[File] aligned_bam_bais = secondary_analysis_core.aligned_bam_bais
    File isocall_isoforms_gtf = secondary_analysis_core.isocall_isoforms_gtf
    File isocall_count_matrix = secondary_analysis_core.isocall_count_matrix
    File isocall_closest_known = secondary_analysis_core.isocall_closest_known
    File pigeon_classification = secondary_analysis_core.pigeon_classification
    File filtered_isoforms_gtf = secondary_analysis_core.filtered_isoforms_gtf
    File pigeon_filtered_classification = secondary_analysis_core.pigeon_filtered_classification
  }
}
