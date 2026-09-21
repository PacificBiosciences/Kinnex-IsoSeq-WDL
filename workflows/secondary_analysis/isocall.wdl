version 1.0

import "../backend_configuration.wdl" as BackendConfiguration
import "../reference_resources/reference_resources.wdl" as ReferenceResources
import "isocall/isocall_core.wdl" as IsocallCore

workflow isocall {
  meta {
    description: "PacBio Kinnex Iso-Seq isocall workflow: profile aligned FLNC BAMs, merge profiles, and run joint isocall call."
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
      isocall_isoforms_gtf: {
        description: "Isocall isoforms GTF"
      },
      isocall_count_matrix: {
        description: "Isocall count matrix"
      },
      isocall_closest_known: {
        description: "Closest known isoforms"
      }
    }
  }

  parameter_meta {
    aligned_bams: {
      description: "Aligned FLNC BAMs"
    }
    aligned_bam_bais: {
      description: "Aligned FLNC BAM BAIs"
    }
    isocall_extra_merged_profile: {
      description: "Optional merged profile that will be merged in addition to the generated profiles"
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
      description: "Typed optional overrides for reference files; a genome override requires an annotation override"
    }
    isocall_config_preset: {
      description: "Isocall calling configuration preset",
      choices: [
        "default",
        "yolo"
      ]
    }
    isocall_config_file: {
      description: "Optional custom Isocall TOML configuration file, which takes precedence over isocall_config_preset"
    }
    output_prefix: {
      description: "Basename prefix for joint isocall outputs"
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
    Array[File] aligned_bams
    Array[File] aligned_bam_bais
    File? isocall_extra_merged_profile
    String ref_name = "GRCh38_gencode49"
    String? reference_container
    ReferenceOverrides reference_overrides = object {
    }
    String isocall_config_preset = "default"
    File? isocall_config_file
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
    resolution_profile = "isocall",
    ref_name = ref_name,
    reference_container = reference_container,
    reference_overrides = reference_overrides,
    runtime_attributes = default_runtime_attributes
  }

  ResolvedReferenceResources reference_resources = resolve_reference_resources.resources

  call IsocallCore.isocall_core as isocall_core { input:
    aligned_bams = aligned_bams,
    aligned_bam_bais = aligned_bam_bais,
    isocall_extra_merged_profile = isocall_extra_merged_profile,
    annotation_gtf_gz = select_first([
      reference_resources.annotation_gtf_gz
    ]),
    genome_fasta = select_first([
      reference_resources.genome_fasta
    ]),
    genome_fasta_index = select_first([
      reference_resources.genome_fasta_index
    ]),
    runtime_attributes = default_runtime_attributes,
    isocall_config_preset = isocall_config_preset,
    isocall_config_file = isocall_config_file,
    output_prefix = output_prefix
  }

  output {
    String workflow_name = "isocall"
    String workflow_version = "0.4.0"
    String? reference_container_uri = resolve_reference_resources.reference_container_uri
    String reference_mode = resolve_reference_resources.reference_mode
    String? base_resource_bundle_version = resolve_reference_resources.base_resource_bundle_version
    File isocall_isoforms_gtf = isocall_core.isocall_isoforms_gtf
    File isocall_count_matrix = isocall_core.isocall_count_matrix
    File isocall_closest_known = isocall_core.isocall_closest_known
  }
}
