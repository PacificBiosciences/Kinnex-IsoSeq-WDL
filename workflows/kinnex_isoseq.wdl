version 1.0

import "backend_configuration.wdl" as BackendConfiguration
import "preprocessing/preprocessing_core.wdl" as PreprocessingCore
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
      reference_name: {
        description: "Reference name"
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
      refine_summary_reports: {
        description: "isoseq refine filter summary reports"
      },
      flnc_dataset_xml: {
        description: "FLNC ConsensusReadSet XML"
      },
      flnc_child_dataset_xmls: {
        description: "Child FLNC ConsensusReadSet XMLs"
      },
      flnc_dataset_bams: {
        description: "Packaged FLNC BAMs"
      },
      flnc_dataset_bam_pbis: {
        description: "Packaged FLNC BAM PBIs"
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
      grouped_flnc_bams: {
        description: "Grouped FLNC BAMs"
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
    hifi_demux_barcodes: {
      description: "SMRTbell barcode FASTA used by optional upstream HiFi demux and barcode validation"
    }
    skera_adapters: {
      description: "Adapter FASTA for skera split"
    }
    barcoded_primers: {
      description: "Primer FASTA shared by lima and isoseq refine"
    }
    consensusreadset_xmls: {
      description: "SMRT Link/internal metadata: optional input ConsensusReadSet XMLs for FLNC dataset XML generation."
    }
    ref_map_file: {
      description: "TSV containing reference genome information for FLNC alignment, isocall, and isoform classification"
    }
    isoseq_require_polya: {
      description: "Pass --require-polya to isoseq refine"
    }
    pigeon_use_polya: {
      description: "Whether to pass the ref_map pigeon_poly_a resource to pigeon classify"
    }
    pigeon_use_cage_peak: {
      description: "Whether to pass the ref_map pigeon_cage_peak_bed resources to pigeon classify"
    }
    pigeon_use_junction: {
      description: "Whether to pass the ref_map pigeon_junction_coverage resources to pigeon classify"
    }
    isocall_extra_merged_profile: {
      description: "Optional merged profile that will be merged in addition to the generated isocall profiles"
    }
    isocall_min_read_fraction: {
      description: "Minimum read fraction for joint isocall calling"
    }
    isocall_max_bundles_per_gene: {
      description: "Maximum bundles per gene for joint isocall calling"
    }
    isocall_min_reads_per_isoform: {
      description: "Minimum reads per isoform for joint isocall calling"
    }
    pigeon_min_ref_length: {
      description: "Minimum reference length for pigeon classify"
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
      description: "Optional PacBio registry for registry-relative task images. If omitted, quay.io/pacbio is used."
    }
  }

  input {
    Array[HiFiSource] hifi_sources
    File biosample_csv
    File hifi_demux_barcodes
    File skera_adapters
    File barcoded_primers
    Array[File]? consensusreadset_xmls
    File ref_map_file
    Boolean isoseq_require_polya = true
    Boolean pigeon_use_polya = true
    Boolean pigeon_use_cage_peak = true
    Boolean pigeon_use_junction = true
    File? isocall_extra_merged_profile
    Float isocall_min_read_fraction = 0.99
    Int isocall_max_bundles_per_gene = 10000
    Int isocall_min_reads_per_isoform = 3
    Int pigeon_min_ref_length = 100
    String output_prefix = "joint"

    # Backend configuration
    String backend = "HPC"
    Int max_retries = 2
    Int add_memory_mb = 0
    Int nproc = 16
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
    hifi_demux_barcodes = hifi_demux_barcodes,
    skera_adapters = skera_adapters,
    barcoded_primers = barcoded_primers,
    consensusreadset_xmls = consensusreadset_xmls,
    runtime_attributes = default_runtime_attributes,
    isoseq_require_polya = isoseq_require_polya
  }

  call SecondaryAnalysisCore.secondary_analysis_core as secondary_analysis_core { input:
    flnc_bams = preprocessing_core.flnc_bams,
    ref_map_file = ref_map_file,
    runtime_attributes = default_runtime_attributes,
    pigeon_use_polya = pigeon_use_polya,
    pigeon_use_cage_peak = pigeon_use_cage_peak,
    pigeon_use_junction = pigeon_use_junction,
    isocall_extra_merged_profile = isocall_extra_merged_profile,
    isocall_min_read_fraction = isocall_min_read_fraction,
    isocall_max_bundles_per_gene = isocall_max_bundles_per_gene,
    isocall_min_reads_per_isoform = isocall_min_reads_per_isoform,
    pigeon_min_ref_length = pigeon_min_ref_length,
    output_prefix = output_prefix
  }

  output {
    String workflow_name = "kinnex_isoseq"
    String workflow_version = "0.2.0"
    String preprocessing_workflow_name = "preprocessing"
    String secondary_analysis_workflow_name = "secondary_analysis"
    String reference_name = secondary_analysis_core.reference_name
    Array[String] source_dataset_names = preprocessing_core.source_dataset_names
    Array[String] preprocessing_dataset_names = preprocessing_core.dataset_names
    Array[File] hifi_demux_datasets = preprocessing_core.hifi_demux_datasets
    Array[File] hifi_demux_bams = preprocessing_core.hifi_demux_bams
    Array[String] flnc_names = preprocessing_core.flnc_names
    Array[File] flnc_bams = preprocessing_core.flnc_bams
    Array[File] flnc_bam_pbis = preprocessing_core.flnc_bam_pbis
    Array[File] refine_summary_reports = preprocessing_core.refine_summary_reports
    File? flnc_dataset_xml = preprocessing_core.flnc_dataset_xml
    Array[File]? flnc_child_dataset_xmls = preprocessing_core.flnc_child_dataset_xmls
    Array[File]? flnc_dataset_bams = preprocessing_core.flnc_dataset_bams
    Array[File]? flnc_dataset_bam_pbis = preprocessing_core.flnc_dataset_bam_pbis
    Array[String] sample_names = secondary_analysis_core.sample_names
    Array[String] sample_prefixes = secondary_analysis_core.sample_prefixes
    Array[Int] group_sizes = secondary_analysis_core.group_sizes
    Array[File] grouped_flnc_bams = secondary_analysis_core.grouped_flnc_bams
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
