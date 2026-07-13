version 1.0

import "../../rna_structs.wdl"
import "tasks.wdl" as Tasks

workflow isoform_classification_core {
  meta {
    description: "Isoform classification stage for PacBio Kinnex Iso-Seq: prepare isoforms for pigeon, then classify and filter."
    outputs: {
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
      description: "Isoform GTF"
    }
    isocall_count_matrix: {
      description: "Isocall count matrix"
    }
    pigeon_resources: {
      description: "Prepared pigeon classification resources"
    }
    pigeon_use_polya: {
      description: "Pass the pigeon polyA resource to pigeon classify"
    }
    pigeon_use_cage_peak: {
      description: "Pass the pigeon CAGE peak resource to pigeon classify"
    }
    pigeon_use_junction: {
      description: "Pass the pigeon junction coverage resource to pigeon classify"
    }
    pigeon_min_ref_length: {
      description: "Minimum reference length for pigeon classify"
    }
    output_prefix: {
      description: "Output prefix"
    }
    prepare_isoforms_threads: {
      description: "CPU threads for preparing isoform GTF inputs"
    }
    prepare_isoforms_mem_gb: {
      description: "Memory allocation in GB for preparing isoform GTF inputs"
    }
    pigeon_classify_threads: {
      description: "CPU threads for pigeon classify"
    }
    pigeon_classify_mem_gb: {
      description: "Memory allocation in GB for pigeon classify"
    }
    pigeon_filter_threads: {
      description: "CPU threads for pigeon filter"
    }
    pigeon_filter_mem_gb: {
      description: "Memory allocation in GB for pigeon filter"
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
    Int pigeon_min_ref_length = 100
    String output_prefix = "joint"
    Int prepare_isoforms_threads = 8
    Int prepare_isoforms_mem_gb = 32
    Int pigeon_classify_threads = 16
    Int pigeon_classify_mem_gb = 64
    Int pigeon_filter_threads = 8
    Int pigeon_filter_mem_gb = 32
    RuntimeAttributes runtime_attributes
  }

  call Tasks.prepare_isoforms_gtf { input:
    isoforms_gtf = isoforms_gtf,
    threads = prepare_isoforms_threads,
    mem_gb = prepare_isoforms_mem_gb,
    runtime_attributes = runtime_attributes,
    output_prefix = output_prefix
  }

  call Tasks.pigeon_classify_isoforms { input:
    isoforms_gtf = prepare_isoforms_gtf.prepared_isoforms_gtf,
    isocall_count_matrix = isocall_count_matrix,
    pigeon_resources = pigeon_resources,
    pigeon_use_polya = pigeon_use_polya,
    pigeon_use_cage_peak = pigeon_use_cage_peak,
    pigeon_use_junction = pigeon_use_junction,
    pigeon_min_ref_length = pigeon_min_ref_length,
    threads = pigeon_classify_threads,
    mem_gb = pigeon_classify_mem_gb,
    runtime_attributes = runtime_attributes,
    output_prefix = output_prefix
  }

  call Tasks.pigeon_filter_isoforms { input:
    isoforms_gtf = prepare_isoforms_gtf.prepared_isoforms_gtf,
    classification = pigeon_classify_isoforms.pigeon_classification,
    junctions = pigeon_classify_isoforms.junctions,
    threads = pigeon_filter_threads,
    mem_gb = pigeon_filter_mem_gb,
    runtime_attributes = runtime_attributes,
    output_prefix = output_prefix
  }

  output {
    File pigeon_classification = pigeon_classify_isoforms.pigeon_classification
    File filtered_isoforms_gtf = pigeon_filter_isoforms.filtered_isoforms_gtf
    File pigeon_filtered_classification = pigeon_filter_isoforms.pigeon_filtered_classification
  }
}
