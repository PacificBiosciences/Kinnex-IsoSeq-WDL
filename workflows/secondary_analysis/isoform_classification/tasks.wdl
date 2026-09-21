version 1.0

import "../../rna_structs.wdl"

task prepare_pigeon_resources {
  meta {
    description: "Prepare the annotation and any enabled, available supplemental resources for pigeon classification."
    outputs: {
      prepared_pigeon_cage_peak_bed: {
        description: "Prepared Pigeon CAGE peak BED when enabled and available"
      },
      prepared_pigeon_cage_peak_bed_index: {
        description: "Prepared Pigeon CAGE peak index when enabled and available"
      },
      prepared_pigeon_junction_coverage: {
        description: "Prepared Pigeon junction coverage when enabled and available"
      },
      prepared_pigeon_junction_coverage_index: {
        description: "Prepared Pigeon junction-coverage index when enabled and available"
      },
      pigeon_resources: {
        description: "Prepared pigeon classification resources"
      }
    }
  }

  parameter_meta {
    annotation_gtf_gz: {
      description: "Compressed annotation GTF"
    }
    genome_fasta: {
      description: "Reference FASTA"
    }
    genome_fasta_index: {
      description: "Reference FASTA index"
    }
    pigeon_poly_a: {
      description: "Optional Pigeon polyA resource"
    }
    pigeon_cage_peak_bed: {
      description: "Optional Pigeon CAGE peak BED"
    }
    pigeon_junction_coverage: {
      description: "Optional Pigeon junction coverage"
    }
    pigeon_use_cage_peak: {
      description: "Whether to prepare the available CAGE peak resource"
    }
    pigeon_use_junction: {
      description: "Whether to prepare the available junction-coverage resource"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    File annotation_gtf_gz
    File genome_fasta
    File genome_fasta_index
    File? pigeon_poly_a
    File? pigeon_cage_peak_bed
    File? pigeon_junction_coverage
    Boolean pigeon_use_cage_peak = true
    Boolean pigeon_use_junction = true
    RuntimeAttributes runtime_attributes
  }

  Int threads = 1
  Int mem_gb = 8
  Int total_mem_mb = (mem_gb * 1024) + runtime_attributes.add_memory_mb
  Boolean prepare_pigeon_cage_peak = pigeon_use_cage_peak && defined(pigeon_cage_peak_bed)
  Boolean prepare_pigeon_junction = pigeon_use_junction && defined(pigeon_junction_coverage)

  command <<<
    set -euo pipefail

    mkdir -p pigeon_reference_inputs
    gzip -dc "~{annotation_gtf_gz}" > pigeon_reference_inputs/annotation.gtf

    pigeon_prepare_args=()
    if [[ "~{prepare_pigeon_cage_peak}" == "true" ]]; then
      cp -- "~{pigeon_cage_peak_bed}" pigeon_reference_inputs/cage_peak.bed
      pigeon_prepare_args+=(pigeon_reference_inputs/cage_peak.bed)
    fi
    if [[ "~{prepare_pigeon_junction}" == "true" ]]; then
      cp -- "~{pigeon_junction_coverage}" pigeon_reference_inputs/junction_coverage.tsv
      pigeon_prepare_args+=(pigeon_reference_inputs/junction_coverage.tsv)
    fi

    pigeon prepare \
      --log-level INFO \
      --log-file pigeon.prepare.log \
      pigeon_reference_inputs/annotation.gtf \
      "${pigeon_prepare_args[@]}"
  >>>

  output {
    File? prepared_pigeon_cage_peak_bed = "pigeon_reference_inputs/cage_peak.sorted.bed"
    File? prepared_pigeon_cage_peak_bed_index = "pigeon_reference_inputs/cage_peak.sorted.bed.pgi"
    File? prepared_pigeon_junction_coverage = "pigeon_reference_inputs/junction_coverage.sorted.tsv"
    File? prepared_pigeon_junction_coverage_index = "pigeon_reference_inputs/junction_coverage.sorted.tsv.pgi"
    PigeonResources pigeon_resources = object {
      annotation_gtf: "pigeon_reference_inputs/annotation.sorted.gtf",
      annotation_gtf_pgi: "pigeon_reference_inputs/annotation.sorted.gtf.pgi",
      genome_fasta: genome_fasta,
      genome_fasta_index: genome_fasta_index,
      pigeon_poly_a: pigeon_poly_a,
      pigeon_cage_peak_bed: prepared_pigeon_cage_peak_bed,
      pigeon_cage_peak_bed_index: prepared_pigeon_cage_peak_bed_index,
      pigeon_junction_coverage: prepared_pigeon_junction_coverage,
      pigeon_junction_coverage_index: prepared_pigeon_junction_coverage_index
    }
  }

  runtime {
    cpu: threads
    memory: total_mem_mb + " MB"
    docker: runtime_attributes.container_registry + "/pigeon@sha256:94cf3cd64f600ae974a7956177be2231c1c9ecb468a1cc4b70b05dcacf59a81c"  # 26.2.0_build2
    maxRetries: runtime_attributes.max_retries
  }
}

task prepare_isoforms_gtf {
  meta {
    description: "Normalize and sort isoform GTF input before pigeon classification."
    outputs: {
      prepared_isoforms_gtf: {
        description: "Prepared isoforms GTF"
      }
    }
  }

  parameter_meta {
    isoforms_gtf: {
      description: "Isoform GTF"
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
    File isoforms_gtf
    String output_prefix = "joint"
    Int threads = 8
    Int mem_gb = 32
    RuntimeAttributes runtime_attributes
  }

  Int total_mem_mb = (mem_gb * 1024) + runtime_attributes.add_memory_mb
  Int effective_threads = if (threads > runtime_attributes.nproc)
    then runtime_attributes.nproc
    else threads

  command <<<
    set -euo pipefail

    isoforms_input="~{isoforms_gtf}"
    unsorted_gtf="~{output_prefix}.isoforms.unsorted.gtf"

    if [[ "${isoforms_input}" == *.gz ]]; then
      gzip -dc "${isoforms_input}" > "${unsorted_gtf}"
    else
      cp "${isoforms_input}" "${unsorted_gtf}"
    fi

    pigeon sort -o "~{output_prefix}.isoforms.gtf" "${unsorted_gtf}"
    rm -f "${unsorted_gtf}"
  >>>

  output {
    File prepared_isoforms_gtf = "~{output_prefix}.isoforms.gtf"
  }

  runtime {
    cpu: effective_threads
    memory: total_mem_mb + " MB"
    docker: runtime_attributes.container_registry + "/pigeon@sha256:94cf3cd64f600ae974a7956177be2231c1c9ecb468a1cc4b70b05dcacf59a81c"  # 26.2.0_build2
    maxRetries: runtime_attributes.max_retries
  }
}

task pigeon_classify_isoforms {
  meta {
    description: "Run pigeon classify on one joint isoform GTF."
    outputs: {
      pigeon_classification: {
        description: "Pigeon classification table"
      },
      junctions: {
        description: "Pigeon junctions table"
      }
    }
  }

  parameter_meta {
    isoforms_gtf: {
      description: "Prepared isoforms GTF"
    }
    isocall_count_matrix: {
      description: "Isocall count matrix"
    }
    pigeon_resources: {
      description: "Prepared pigeon classification resources"
    }
    pigeon_use_polya: {
      description: "Whether to pass the polyA resource to pigeon classify when the resource is available"
    }
    pigeon_use_cage_peak: {
      description: "Whether to pass the CAGE peak resource to pigeon classify when the resource is available"
    }
    pigeon_use_junction: {
      description: "Whether to pass the junction-coverage resource to pigeon classify when the resource is available"
    }
    pigeon_min_ref_length: {
      description: "Minimum reference length for pigeon classify"
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
    File isoforms_gtf
    File isocall_count_matrix
    PigeonResources pigeon_resources
    Boolean pigeon_use_polya = true
    Boolean pigeon_use_cage_peak = true
    Boolean pigeon_use_junction = true
    Int pigeon_min_ref_length
    String output_prefix = "joint"
    Int threads = 16
    Int mem_gb = 64
    RuntimeAttributes runtime_attributes
  }

  Int total_mem_mb = (mem_gb * 1024) + runtime_attributes.add_memory_mb
  Int effective_threads = if (threads > runtime_attributes.nproc)
    then runtime_attributes.nproc
    else threads

  Boolean effective_pigeon_use_polya = pigeon_use_polya && defined(pigeon_resources.pigeon_poly_a)
  Boolean effective_pigeon_use_cage_peak = (pigeon_use_cage_peak && defined(pigeon_resources.pigeon_cage_peak_bed) && defined(pigeon_resources.pigeon_cage_peak_bed_index))
  Boolean effective_pigeon_use_junction = (pigeon_use_junction && defined(pigeon_resources.pigeon_junction_coverage) && defined(pigeon_resources.pigeon_junction_coverage_index))
  String annotation_gtf_basename = basename(pigeon_resources.annotation_gtf)
  String genome_fasta_basename = basename(pigeon_resources.genome_fasta)
  String pigeon_cage_peak_bed_basename = if effective_pigeon_use_cage_peak
    then basename(select_first([
      pigeon_resources.pigeon_cage_peak_bed
    ]))
    else ""
  String pigeon_junction_coverage_basename = if effective_pigeon_use_junction
    then basename(select_first([
      pigeon_resources.pigeon_junction_coverage
    ]))
    else ""

  command <<<
    set -euo pipefail

    # pigeon classify expects annotation/reference indices next to the primary files.
    ln --symbolic "~{pigeon_resources.annotation_gtf}" .
    ln --symbolic "~{pigeon_resources.annotation_gtf_pgi}" "~{annotation_gtf_basename}.pgi"
    ln --symbolic "~{pigeon_resources.genome_fasta}" .
    ln --symbolic "~{pigeon_resources.genome_fasta_index}" "~{genome_fasta_basename}.fai"

    pigeon_classify_args=()
    if [[ "~{effective_pigeon_use_polya}" == "true" ]]; then
      pigeon_classify_args+=(--poly-a "~{pigeon_resources.pigeon_poly_a}")
    fi
    if [[ "~{effective_pigeon_use_cage_peak}" == "true" ]]; then
      ln --symbolic "~{pigeon_resources.pigeon_cage_peak_bed}" .
      ln --symbolic "~{pigeon_resources.pigeon_cage_peak_bed_index}" "~{pigeon_cage_peak_bed_basename}.pgi"
      pigeon_classify_args+=(--cage-peak "~{pigeon_cage_peak_bed_basename}")
    fi
    if [[ "~{effective_pigeon_use_junction}" == "true" ]]; then
      ln --symbolic "~{pigeon_resources.pigeon_junction_coverage}" .
      ln --symbolic "~{pigeon_resources.pigeon_junction_coverage_index}" "~{pigeon_junction_coverage_basename}.pgi"
      pigeon_classify_args+=(--coverage "~{pigeon_junction_coverage_basename}")
    fi

    pigeon classify \
      --num-threads "~{effective_threads}" \
      --log-level INFO \
      --log-file "~{output_prefix}.pigeon.log" \
      --flnc "~{isocall_count_matrix}" \
      "~{isoforms_gtf}" \
      "~{annotation_gtf_basename}" \
      "~{genome_fasta_basename}" \
      --out-prefix "~{output_prefix}.pigeon" \
      "${pigeon_classify_args[@]}" \
      --min-ref-length "~{pigeon_min_ref_length}"
  >>>

  output {
    File pigeon_classification = "~{output_prefix}.pigeon_classification.txt"
    File junctions = "~{output_prefix}.pigeon_junctions.txt"
  }

  runtime {
    cpu: effective_threads
    memory: total_mem_mb + " MB"
    docker: runtime_attributes.container_registry + "/pigeon@sha256:94cf3cd64f600ae974a7956177be2231c1c9ecb468a1cc4b70b05dcacf59a81c"  # 26.2.0_build2
    maxRetries: runtime_attributes.max_retries
  }
}

task pigeon_filter_isoforms {
  meta {
    description: "Run pigeon filter to produce the filtered lite isoform set."
    outputs: {
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
      description: "Prepared isoforms GTF"
    }
    classification: {
      description: "Pigeon classification table"
    }
    junctions: {
      description: "Pigeon junctions table"
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
    File isoforms_gtf
    File classification
    File junctions
    String output_prefix = "joint"
    Int threads = 8
    Int mem_gb = 32
    RuntimeAttributes runtime_attributes
  }

  Int total_mem_mb = (mem_gb * 1024) + runtime_attributes.add_memory_mb
  Int effective_threads = if (threads > runtime_attributes.nproc)
    then runtime_attributes.nproc
    else threads

  command <<<
    set -euo pipefail

    isoforms_input="~{isoforms_gtf}"
    classification_input="~{classification}"
    junctions_input="~{junctions}"
    expected_isoforms="~{output_prefix}.isoforms.gtf"
    expected_classification="~{output_prefix}.pigeon_classification.txt"
    expected_junctions="~{output_prefix}.pigeon_junctions.txt"

    # pigeon filter derives the sibling junctions input and filtered output
    # basenames from these filenames, so expose stable local names via links.
    ln -sf "${isoforms_input}" "${expected_isoforms}"
    ln -sf "${classification_input}" "${expected_classification}"
    ln -sf "${junctions_input}" "${expected_junctions}"

    pigeon filter \
      --num-threads "~{effective_threads}" \
      --log-level INFO \
      --log-file "~{output_prefix}.pigeon.filter.log" \
      --isoforms "${expected_isoforms}" \
      "${expected_classification}"
  >>>

  output {
    File filtered_isoforms_gtf = "~{output_prefix}.isoforms.filtered_lite.gtf"
    File pigeon_filtered_classification = "~{output_prefix}.pigeon_classification.filtered_lite_classification.txt"
  }

  runtime {
    cpu: effective_threads
    memory: total_mem_mb + " MB"
    docker: runtime_attributes.container_registry + "/pigeon@sha256:94cf3cd64f600ae974a7956177be2231c1c9ecb468a1cc4b70b05dcacf59a81c"  # 26.2.0_build2
    maxRetries: runtime_attributes.max_retries
  }
}
