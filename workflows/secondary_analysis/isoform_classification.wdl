version 1.0

import "../backend_configuration.wdl" as BackendConfiguration
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
      reference_name: {
        description: "Reference name"
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
    ref_map_file: {
      description: "TSV containing the required annotation and reference genome files"
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
    File ref_map_file
    Boolean pigeon_use_polya = true
    Boolean pigeon_use_cage_peak = true
    Boolean pigeon_use_junction = true
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

  Map[String, String] reference_mapping = read_map(ref_map_file)

  call IsoformClassificationTasks.prepare_pigeon_resources { input:
    annotation_gtf_gz = reference_mapping["annotation_gtf_gz"],  # !FileCoercion
    genome_fasta = reference_mapping["genome_fasta"],  # !FileCoercion
    genome_fasta_index = reference_mapping["genome_fasta_index"],  # !FileCoercion
    pigeon_poly_a = reference_mapping["pigeon_poly_a"],  # !FileCoercion
    pigeon_cage_peak_bed = reference_mapping["pigeon_cage_peak_bed"],  # !FileCoercion
    pigeon_junction_coverage = reference_mapping["pigeon_junction_coverage"],  # !FileCoercion
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
    String workflow_version = "0.1.0"
    String reference_name = reference_mapping["name"]
    File pigeon_classification = isoform_classification_core.pigeon_classification
    File filtered_isoforms_gtf = isoform_classification_core.filtered_isoforms_gtf
    File pigeon_filtered_classification = isoform_classification_core.pigeon_filtered_classification
  }
}
