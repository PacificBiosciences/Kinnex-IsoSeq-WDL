version 1.0

import "../../rna_structs.wdl"

task skera_split_hifi {
  meta {
    description: "Split one HiFi BAM into S-reads with skera."
    outputs: {
      dataset_name_out: {
        description: "Dataset name"
      },
      segmented_bam: {
        description: "Segmented BAM"
      }
    }
  }

  parameter_meta {
    dataset_name: {
      description: "Dataset name"
    }
    hifi_bam: {
      description: "HiFi BAM"
    }
    segmentation_adapters: {
      description: "Segmentation adapter FASTA"
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
    String dataset_name
    File hifi_bam
    File segmentation_adapters
    Int threads = 16
    Int mem_gb = 32
    RuntimeAttributes runtime_attributes
  }

  Int total_mem_mb = (mem_gb * 1024) + runtime_attributes.add_memory_mb
  Int effective_threads = if (threads > runtime_attributes.nproc)
    then runtime_attributes.nproc
    else threads

  String output_prefix = dataset_name

  command <<<
    set -euo pipefail

    skera split \
      --num-threads ~{effective_threads} \
      --log-level INFO \
      --log-file "~{output_prefix}.skera.log" \
      "~{hifi_bam}" \
      "~{segmentation_adapters}" \
      "~{output_prefix}.sreads.bam"
  >>>

  output {
    String dataset_name_out = dataset_name
    File segmented_bam = "~{output_prefix}.sreads.bam"
  }

  runtime {
    cpu: effective_threads
    memory: total_mem_mb + " MB"
    docker: runtime_attributes.container_registry + "/skera@sha256:f823bd64f2beec351cae82e67c9c6f257b4896ddd59f8b2050568d59f3a165f6"  # 1.4.0_build3
    maxRetries: runtime_attributes.max_retries
  }
}
