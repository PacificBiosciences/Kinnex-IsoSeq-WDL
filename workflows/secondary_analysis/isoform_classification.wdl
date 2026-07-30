version 1.0

import "../backend_configuration.wdl" as BackendConfiguration
import "../reference_resources/reference_resources.wdl" as ReferenceResources
import "isoform_classification/isoform_classification_core.wdl" as IsoformClassificationCore
import "isoform_classification/tasks.wdl" as IsoformClassificationTasks

workflow isoform_classification {
  meta {
    description: "PacBio Kinnex Iso-Seq isoform classification workflow: prepare isoforms, run pigeon classify, then run pigeon filter."
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
    isoforms_gtf: {
      description: "Isoform GTF from the isocall stage"
    }
    isocall_count_matrix: {
      description: "Supporting-read count matrix from the isocall stage"
    }
    reference_container: {
      description: "Optional immutable reference-container URI used for defaults"
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
    pigeon_min_ref_length: {
      description: "Minimum reference length passed to pigeon classify"
    }
    output_prefix: {
      description: "Basename prefix for isoform classification outputs"
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
    File isoforms_gtf
    File isocall_count_matrix
    String? reference_container
    ReferenceOverrides reference_overrides = object {
    }
    Boolean pigeon_use_polya = true
    Boolean pigeon_use_cage_peak = true
    Boolean pigeon_use_junction = true
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
    resolution_profile = "isoform_classification",
    reference_container = reference_container,
    reference_overrides = reference_overrides,
    runtime_attributes = default_runtime_attributes
  }

  ResolvedReferenceResources reference_resources = resolve_reference_resources.resources

  call IsoformClassificationTasks.prepare_pigeon_resources { input:
    annotation_gtf_gz = select_first([
      reference_resources.annotation_gtf_gz
    ]),
    genome_fasta = select_first([
      reference_resources.genome_fasta
    ]),
    genome_fasta_index = select_first([
      reference_resources.genome_fasta_index
    ]),
    pigeon_poly_a = reference_resources.pigeon_poly_a,
    pigeon_cage_peak_bed = reference_resources.pigeon_cage_peak_bed,
    pigeon_junction_coverage = reference_resources.pigeon_junction_coverage,
    pigeon_use_cage_peak = pigeon_use_cage_peak,
    pigeon_use_junction = pigeon_use_junction,
    runtime_attributes = default_runtime_attributes
  }

  call IsoformClassificationCore.isoform_classification_core as isoform_classification_core { input:
    isoforms_gtf = isoforms_gtf,
    isocall_count_matrix = isocall_count_matrix,
    pigeon_resources = prepare_pigeon_resources.pigeon_resources,
    pigeon_use_polya = pigeon_use_polya,
    pigeon_use_cage_peak = pigeon_use_cage_peak,
    pigeon_use_junction = pigeon_use_junction,
    runtime_attributes = default_runtime_attributes,
    pigeon_min_ref_length = pigeon_min_ref_length,
    output_prefix = output_prefix
  }

  output {
    String workflow_name = "isoform_classification"
    String workflow_version = "0.3.0"
    String? reference_container_uri = resolve_reference_resources.reference_container_uri
    String reference_mode = resolve_reference_resources.reference_mode
    String? base_resource_bundle_version = resolve_reference_resources.base_resource_bundle_version
    File pigeon_classification = isoform_classification_core.pigeon_classification
    File filtered_isoforms_gtf = isoform_classification_core.filtered_isoforms_gtf
    File pigeon_filtered_classification = isoform_classification_core.pigeon_filtered_classification
  }
}
