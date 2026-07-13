version 1.0

import "../../rna_structs.wdl"
import "tasks.wdl" as Tasks

workflow isocall_core {
  meta {
    description: "Isocall stage for PacBio Kinnex Iso-Seq: profile aligned FLNC BAMs, prepare known isoforms, merge profiles, then run joint isocall call."
    outputs: {
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
      description: "Optional extra merged isocall profile"
    }
    annotation_gtf_gz: {
      description: "Compressed annotation GTF"
    }
    genome_fasta: {
      description: "Reference FASTA"
    }
    genome_fasta_index: {
      description: "Reference FASTA index"
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
    output_prefix: {
      description: "Output prefix"
    }
    isocall_profile_threads: {
      description: "CPU threads for isocall profile"
    }
    isocall_profile_mem_gb: {
      description: "Memory allocation in GB for isocall profile"
    }
    isocall_prep_isoforms_threads: {
      description: "CPU threads for isocall prep-isoforms"
    }
    isocall_prep_isoforms_mem_gb: {
      description: "Memory allocation in GB for isocall prep-isoforms"
    }
    isocall_merge_profiles_threads: {
      description: "CPU threads for isocall merge"
    }
    isocall_merge_profiles_mem_gb: {
      description: "Memory allocation in GB for isocall merge"
    }
    isocall_call_threads: {
      description: "CPU threads for isocall call"
    }
    isocall_call_mem_gb: {
      description: "Memory allocation in GB for isocall call"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    Array[File] aligned_bams
    Array[File] aligned_bam_bais
    File? isocall_extra_merged_profile
    File annotation_gtf_gz
    File genome_fasta
    File genome_fasta_index
    Float isocall_min_read_fraction = 0.99
    Int isocall_max_bundles_per_gene = 10000
    Int isocall_min_reads_per_isoform = 3
    String output_prefix = "joint"
    Int isocall_profile_threads = 4
    Int isocall_profile_mem_gb = 32
    Int isocall_prep_isoforms_threads = 4
    Int isocall_prep_isoforms_mem_gb = 16
    Int isocall_merge_profiles_threads = 4
    Int isocall_merge_profiles_mem_gb = 16
    Int isocall_call_threads = 16
    Int isocall_call_mem_gb = 32
    RuntimeAttributes runtime_attributes
  }

  scatter (sample_index in range(length(aligned_bams))) {
    File aligned_bam = aligned_bams[sample_index]
    File aligned_bam_index = aligned_bam_bais[sample_index]

    call Tasks.isocall_profile { input:
      aligned_bam = aligned_bam,
      aligned_bam_index = aligned_bam_index,
      threads = isocall_profile_threads,
      mem_gb = isocall_profile_mem_gb,
      runtime_attributes = runtime_attributes
    }
  }

  call Tasks.isocall_prep_isoforms { input:
    annotation_gtf_gz = annotation_gtf_gz,
    threads = isocall_prep_isoforms_threads,
    mem_gb = isocall_prep_isoforms_mem_gb,
    runtime_attributes = runtime_attributes,
    output_prefix = output_prefix
  }

  call Tasks.isocall_merge_profiles { input:
    profiles = isocall_profile.profile,
    isocall_extra_merged_profile = isocall_extra_merged_profile,
    threads = isocall_merge_profiles_threads,
    mem_gb = isocall_merge_profiles_mem_gb,
    runtime_attributes = runtime_attributes,
    output_prefix = output_prefix
  }

  call Tasks.isocall_call as isocall_call_joint { input:
    merged_profile = isocall_merge_profiles.isocall_merged_profile,
    known_isoforms_model = isocall_prep_isoforms.known_isoforms_model,
    genome_fasta = genome_fasta,
    genome_fasta_index = genome_fasta_index,
    isocall_min_read_fraction = isocall_min_read_fraction,
    isocall_max_bundles_per_gene = isocall_max_bundles_per_gene,
    isocall_min_reads_per_isoform = isocall_min_reads_per_isoform,
    threads = isocall_call_threads,
    mem_gb = isocall_call_mem_gb,
    runtime_attributes = runtime_attributes,
    output_prefix = output_prefix
  }

  output {
    File isocall_isoforms_gtf = isocall_call_joint.isocall_isoforms_gtf
    File isocall_count_matrix = isocall_call_joint.isocall_count_matrix
    File isocall_closest_known = isocall_call_joint.isocall_closest_known
  }
}
