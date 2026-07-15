version 1.0

import "backend_configuration.wdl" as BackendConfiguration
import "secondary_analysis/secondary_analysis_core.wdl" as SecondaryAnalysisCore

workflow secondary_analysis {
  meta {
    description: "PacBio Kinnex Iso-Seq pipeline: align FLNC BAMs, run joint isocall, then run isoform classification."
    outputs: {
      workflow_name: {
        description: "Workflow name"
      },
      workflow_version: {
        description: "Workflow version"
      },
      reference_name: {
        description: "Reference name"
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
    flnc_bams: {
      description: "Input FLNC BAMs"
    }
    ref_map_file: {
      description: "TSV containing reference genome information for all implemented RNA stages"
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
      description: "Backend where the workflow will be executed; only HPC is supported"
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
      description: "Optional PacBio registry for registry-relative task images; if omitted, quay.io/pacbio is used"
    }
  }

  input {
    Array[File] flnc_bams
    File ref_map_file
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

  call SecondaryAnalysisCore.secondary_analysis_core as secondary_analysis_core { input:
    flnc_bams = flnc_bams,
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
    String workflow_name = "secondary_analysis"
    String workflow_version = "0.2.0"
    String reference_name = secondary_analysis_core.reference_name
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
