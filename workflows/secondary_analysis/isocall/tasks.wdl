version 1.0

import "../../rna_structs.wdl"

task isocall_profile {
  meta {
    description: "Generate one isocall profile from one aligned FLNC BAM."
    outputs: {
      profile: {
        description: "Isocall profile"
      }
    }
  }

  parameter_meta {
    aligned_bam: {
      description: "Aligned FLNC BAM"
    }
    aligned_bam_index: {
      description: "Aligned FLNC BAM index"
    }
    threads: {
      description: "CPU threads"
    }
    mem_gb: {
      description: "Memory allocation in GB"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    File aligned_bam
    File aligned_bam_index
    Int threads = 4
    Int mem_gb = 32
    RuntimeAttributes runtime_attributes
  }

  Int total_mem_mb = (mem_gb * 1024) + runtime_attributes.add_memory_mb
  Int effective_threads = if (threads > runtime_attributes.nproc)
    then runtime_attributes.nproc
    else threads

  String output_prefix = basename(aligned_bam, ".bam")

  command <<<
    set -euo pipefail

    ln -sf "~{aligned_bam}" profile_input.bam
    ln -sf "~{aligned_bam_index}" profile_input.bam.bai

    isocall profile \
      --reads profile_input.bam \
      --output "~{output_prefix}.isocall_profile.gz"
  >>>

  output {
    File profile = "~{output_prefix}.isocall_profile.gz"
  }

  runtime {
    cpu: effective_threads
    memory: total_mem_mb + " MB"
    docker: runtime_attributes.container_registry + "/isocall@sha256:34cbf179e342515ddfdd65e7f170ab76e4f6fd3993fcd06f236c13712eb8e882"  # 1.3.0_build1
    maxRetries: runtime_attributes.max_retries
  }
}

task isocall_prep_isoforms {
  meta {
    description: "Prepare known isoforms from a genome annotation before joint isocall calling."
    outputs: {
      known_isoforms_model: {
        description: "Known isoforms model"
      }
    }
  }

  parameter_meta {
    annotation_gtf_gz: {
      description: "Compressed annotation GTF"
    }
    output_prefix: {
      description: "Output prefix"
    }
    threads: {
      description: "CPU threads"
    }
    mem_gb: {
      description: "Memory allocation in GB"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    File annotation_gtf_gz
    String output_prefix = "joint"
    Int threads = 4
    Int mem_gb = 16
    RuntimeAttributes runtime_attributes
  }

  Int total_mem_mb = (mem_gb * 1024) + runtime_attributes.add_memory_mb
  Int effective_threads = if (threads > runtime_attributes.nproc)
    then runtime_attributes.nproc
    else threads

  command <<<
    set -euo pipefail

    isocall prep-isoforms \
      --gtf "~{annotation_gtf_gz}" \
      --output "~{output_prefix}.isocall_model.gz" \
      2> "~{output_prefix}.prep_isoforms.log"
  >>>

  output {
    File known_isoforms_model = "~{output_prefix}.isocall_model.gz"
  }

  runtime {
    cpu: effective_threads
    memory: total_mem_mb + " MB"
    docker: runtime_attributes.container_registry + "/isocall@sha256:34cbf179e342515ddfdd65e7f170ab76e4f6fd3993fcd06f236c13712eb8e882"  # 1.3.0_build1
    maxRetries: runtime_attributes.max_retries
  }
}

task isocall_merge_profiles {
  meta {
    description: "Merge isocall profile outputs, optionally including a provided merged profile."
    outputs: {
      isocall_merged_profile: {
        description: "Merged isocall profile"
      }
    }
  }

  parameter_meta {
    profiles: {
      description: "Isocall profiles"
    }
    isocall_extra_merged_profile: {
      description: "Optional extra merged isocall profile"
    }
    output_prefix: {
      description: "Output prefix"
    }
    threads: {
      description: "CPU threads"
    }
    mem_gb: {
      description: "Memory allocation in GB"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    Array[File] profiles
    File? isocall_extra_merged_profile
    String output_prefix = "joint"
    Int threads = 4
    Int mem_gb = 16
    RuntimeAttributes runtime_attributes
  }

  Int total_mem_mb = (mem_gb * 1024) + runtime_attributes.add_memory_mb
  Int effective_threads = if (threads > runtime_attributes.nproc)
    then runtime_attributes.nproc
    else threads

  Array[File] all_profiles = flatten([
    profiles,
    select_all([
      isocall_extra_merged_profile
    ])
  ])

  command <<<
    set -euo pipefail

    printf '%s\n' "~{sep="\" \"" all_profiles}" > all_profiles.txt
    mapfile -t all_profiles < all_profiles.txt

    isocall merge \
      --profiles "${all_profiles[@]}" \
      --output "~{output_prefix}.merged_profile.gz" \
      2> "~{output_prefix}.isocall_merge.log"
  >>>

  output {
    File isocall_merged_profile = "~{output_prefix}.merged_profile.gz"
  }

  runtime {
    cpu: effective_threads
    memory: total_mem_mb + " MB"
    docker: runtime_attributes.container_registry + "/isocall@sha256:34cbf179e342515ddfdd65e7f170ab76e4f6fd3993fcd06f236c13712eb8e882"  # 1.3.0_build1
    maxRetries: runtime_attributes.max_retries
  }
}

task isocall_call {
  meta {
    description: "Run isocall calling from a merged profile and prepared isoform model."
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
    merged_profile: {
      description: "Merged isocall profile"
    }
    known_isoforms_model: {
      description: "Known isoforms model"
    }
    genome_fasta: {
      description: "Reference FASTA"
    }
    genome_fasta_index: {
      description: "Reference FASTA index"
    }
    isocall_config_preset: {
      description: "Isocall calling configuration preset"
    }
    isocall_config_file: {
      description: "Optional custom Isocall TOML configuration file, which takes precedence over the preset"
    }
    output_prefix: {
      description: "Output prefix"
    }
    threads: {
      description: "CPU threads"
    }
    mem_gb: {
      description: "Memory allocation in GB"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    File merged_profile
    File known_isoforms_model
    File genome_fasta
    File genome_fasta_index
    String isocall_config_preset = "default"
    File? isocall_config_file
    String output_prefix = "joint"
    Int threads = 16
    Int mem_gb = 32
    RuntimeAttributes runtime_attributes
  }

  Int total_mem_mb = (mem_gb * 1024) + runtime_attributes.add_memory_mb
  Int effective_threads = if (threads > runtime_attributes.nproc)
    then runtime_attributes.nproc
    else threads

  command <<<
    set -euo pipefail

    # isocall expects the FASTA index next to the reference as <reference>.fai.
    ln --symbolic "~{genome_fasta}" "reference.fa"
    ln --symbolic "~{genome_fasta_index}" "reference.fa.fai"

    isocall_config="~{isocall_config_preset}"
    if [[ "~{defined(isocall_config_file)}" == "true" ]]; then
      isocall_config="~{isocall_config_file}"
    fi

    isocall call \
      --merged-profile "~{merged_profile}" \
      --threads "~{effective_threads}" \
      --known-isoforms "~{known_isoforms_model}" \
      --output-prefix "~{output_prefix}" \
      --reference "reference.fa" \
      --config "${isocall_config}" \
      > "~{output_prefix}.isoforms.gtf.gz" \
      2> "~{output_prefix}.isocall_call.log"
  >>>

  output {
    File isocall_isoforms_gtf = "~{output_prefix}.isoforms.gtf.gz"
    File isocall_count_matrix = "~{output_prefix}.count_matrix.txt"
    File isocall_closest_known = "~{output_prefix}.closest_known.txt"
  }

  runtime {
    cpu: effective_threads
    memory: total_mem_mb + " MB"
    docker: runtime_attributes.container_registry + "/isocall@sha256:34cbf179e342515ddfdd65e7f170ab76e4f6fd3993fcd06f236c13712eb8e882"  # 1.3.0_build1
    maxRetries: runtime_attributes.max_retries
  }
}
