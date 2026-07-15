version 1.0

import "../backend_configuration.wdl" as BackendConfiguration
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
      reference_name: {
        description: "Reference name"
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
    ref_map_file: {
      description: "TSV containing reference genome information for isocall, including annotation_gtf_gz"
    }
    isocall_min_read_fraction: {
      description: "Minimum fraction of reads required for a reported isoform"
    }
    isocall_max_bundles_per_gene: {
      description: "Maximum number of splice bundles to evaluate per gene"
    }
    isocall_min_reads_per_isoform: {
      description: "Minimum reads required for a reported isoform"
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
      description: "Optional PacBio registry for registry-relative task images; if omitted, quay.io/pacbio is used"
    }
  }

  input {
    Array[File] aligned_bams
    Array[File] aligned_bam_bais
    File? isocall_extra_merged_profile
    File ref_map_file
    Float isocall_min_read_fraction = 0.99
    Int isocall_max_bundles_per_gene = 10000
    Int isocall_min_reads_per_isoform = 3
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

  call IsocallCore.isocall_core as isocall_core { input:
    aligned_bams = aligned_bams,
    aligned_bam_bais = aligned_bam_bais,
    isocall_extra_merged_profile = isocall_extra_merged_profile,
    annotation_gtf_gz = reference_mapping["annotation_gtf_gz"],  # !FileCoercion
    genome_fasta = reference_mapping["genome_fasta"],  # !FileCoercion
    genome_fasta_index = reference_mapping["genome_fasta_index"],  # !FileCoercion
    runtime_attributes = default_runtime_attributes,
    isocall_min_read_fraction = isocall_min_read_fraction,
    isocall_max_bundles_per_gene = isocall_max_bundles_per_gene,
    isocall_min_reads_per_isoform = isocall_min_reads_per_isoform,
    output_prefix = output_prefix
  }

  output {
    String workflow_name = "isocall"
    String workflow_version = "0.2.0"
    String reference_name = reference_mapping["name"]
    File isocall_isoforms_gtf = isocall_core.isocall_isoforms_gtf
    File isocall_count_matrix = isocall_core.isocall_count_matrix
    File isocall_closest_known = isocall_core.isocall_closest_known
  }
}
