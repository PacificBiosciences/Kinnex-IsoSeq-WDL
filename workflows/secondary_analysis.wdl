version 1.0

import "backend_configuration.wdl" as BackendConfiguration
import "reference_resources/reference_resources.wdl" as ReferenceResources
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
      reference_container_uri: {
        description: "Immutable reference-container URI when container defaults were used"
      },
      reference_mode: {
        description: "Reference selection mode: container, hybrid, or custom"
      },
      base_resource_bundle_version: {
        description: "Base container resource-bundle version when container defaults were used"
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
    flnc_bams: {
      description: "Input FLNC BAMs"
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
    pigeon_use_polya: {
      description: "Whether to pass the polyA resource to pigeon classify when the resource is available"
    }
    pigeon_use_cage_peak: {
      description: "Whether to prepare and pass the CAGE peak resource when it is available"
    }
    pigeon_use_junction: {
      description: "Whether to prepare and pass the junction-coverage resource when it is available"
    }
    isocall_extra_merged_profile: {
      description: "Optional merged profile that will be merged in addition to the generated isocall profiles"
    }
    isocall_config_preset: {
      description: "Isocall calling configuration preset",
      choices: [
        "default",
        "yolo"
      ]
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
      description: "Optional PacBio registry for registry-relative task images and named reference containers; if omitted, quay.io/pacbio is used"
    }
  }

  input {
    Array[File] flnc_bams
    String ref_name = "GRCh38_gencode49"
    String? reference_container
    ReferenceOverrides reference_overrides = object {
    }
    Boolean pigeon_use_polya = true
    Boolean pigeon_use_cage_peak = true
    Boolean pigeon_use_junction = true
    File? isocall_extra_merged_profile
    String isocall_config_preset = "default"
    Int pigeon_min_ref_length = 100
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
    resolution_profile = "secondary_analysis",
    ref_name = ref_name,
    reference_container = reference_container,
    reference_overrides = reference_overrides,
    runtime_attributes = default_runtime_attributes
  }

  ResolvedReferenceResources reference_resources = resolve_reference_resources.resources

  call SecondaryAnalysisCore.secondary_analysis_core as secondary_analysis_core { input:
    flnc_bams = flnc_bams,
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
    pigeon_use_polya = pigeon_use_polya,
    pigeon_use_cage_peak = pigeon_use_cage_peak,
    pigeon_use_junction = pigeon_use_junction,
    isocall_extra_merged_profile = isocall_extra_merged_profile,
    isocall_config_preset = isocall_config_preset,
    pigeon_min_ref_length = pigeon_min_ref_length,
    output_prefix = output_prefix
  }

  output {
    String workflow_name = "secondary_analysis"
    String workflow_version = "0.4.0"
    String? reference_container_uri = resolve_reference_resources.reference_container_uri
    String reference_mode = resolve_reference_resources.reference_mode
    String? base_resource_bundle_version = resolve_reference_resources.base_resource_bundle_version
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
