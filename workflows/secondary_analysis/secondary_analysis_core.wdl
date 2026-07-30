version 1.0

import "alignment/flnc_alignment_tasks.wdl" as FlncAlignmentTasks
import "isocall/isocall_core.wdl" as IsocallCore
import "isoform_classification/isoform_classification_core.wdl" as IsoformClassificationCore
import "isoform_classification/tasks.wdl" as IsoformClassificationTasks

workflow secondary_analysis_core {
  meta {
    description: "PacBio Kinnex Iso-Seq pipeline implementation: align FLNC BAMs, run joint isocall, then run isoform classification."
    outputs: {
      workflow_name: {
        description: "Workflow name"
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
      description: "FLNC BAMs"
    }
    genome_fasta: {
      description: "Reference genome FASTA"
    }
    genome_fasta_index: {
      description: "Reference FASTA index"
    }
    annotation_gtf_gz: {
      description: "Compressed annotation GTF"
    }
    pigeon_poly_a: {
      description: "Optional Pigeon polyA motif list"
    }
    pigeon_cage_peak_bed: {
      description: "Optional Pigeon CAGE peak BED"
    }
    pigeon_junction_coverage: {
      description: "Optional Pigeon junction coverage table"
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
      description: "Optional extra merged isocall profile"
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
    pigeon_min_ref_length: {
      description: "Minimum reference length for pigeon classify"
    }
    output_prefix: {
      description: "Output prefix"
    }
    pbsamoa_merge_compression: {
      description: "Compression level for pbsamoa merge"
    }
    pbsamoa_merge_threads: {
      description: "CPU threads for pbsamoa merge"
    }
    pbsamoa_merge_mem_gb: {
      description: "Memory allocation in GB for pbsamoa merge"
    }
    pbmm2_index_threads: {
      description: "CPU threads for pbmm2 reference indexing"
    }
    pbmm2_index_mem_gb: {
      description: "Memory allocation in GB for pbmm2 reference indexing"
    }
    pbmm2_align_threads: {
      description: "CPU threads for pbmm2 alignment"
    }
    pbmm2_align_mem_gb: {
      description: "Memory allocation in GB for pbmm2 alignment"
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
    Array[File] flnc_bams
    File genome_fasta
    File genome_fasta_index
    File annotation_gtf_gz
    File? pigeon_poly_a
    File? pigeon_cage_peak_bed
    File? pigeon_junction_coverage
    Boolean pigeon_use_polya = true
    Boolean pigeon_use_cage_peak = true
    Boolean pigeon_use_junction = true
    File? isocall_extra_merged_profile
    Float isocall_min_read_fraction = 0.99
    Int isocall_max_bundles_per_gene = 10000
    Int isocall_min_reads_per_isoform = 3
    Int pigeon_min_ref_length = 100
    String output_prefix = "joint"
    Int pbsamoa_merge_compression = 6
    Int pbsamoa_merge_threads = 32
    Int pbsamoa_merge_mem_gb = 16
    Int pbmm2_index_threads = 8
    Int pbmm2_index_mem_gb = 32
    Int pbmm2_align_threads = 32
    Int pbmm2_align_mem_gb = 64
    Int isocall_profile_threads = 4
    Int isocall_profile_mem_gb = 32
    Int isocall_prep_isoforms_threads = 4
    Int isocall_prep_isoforms_mem_gb = 16
    Int isocall_merge_profiles_threads = 4
    Int isocall_merge_profiles_mem_gb = 16
    Int isocall_call_threads = 16
    Int isocall_call_mem_gb = 32
    Int prepare_isoforms_threads = 8
    Int prepare_isoforms_mem_gb = 32
    Int pigeon_classify_threads = 16
    Int pigeon_classify_mem_gb = 64
    Int pigeon_filter_threads = 8
    Int pigeon_filter_mem_gb = 32
    RuntimeAttributes runtime_attributes
  }

  call IsoformClassificationTasks.prepare_pigeon_resources { input:
    annotation_gtf_gz = annotation_gtf_gz,
    genome_fasta = genome_fasta,
    genome_fasta_index = genome_fasta_index,
    pigeon_poly_a = pigeon_poly_a,
    pigeon_cage_peak_bed = pigeon_cage_peak_bed,
    pigeon_junction_coverage = pigeon_junction_coverage,
    pigeon_use_cage_peak = pigeon_use_cage_peak,
    pigeon_use_junction = pigeon_use_junction,
    runtime_attributes = runtime_attributes
  }

  call FlncAlignmentTasks.group_flnc_bams_by_sm { input:
    flnc_bams = flnc_bams,
    runtime_attributes = runtime_attributes
  }

  call FlncAlignmentTasks.create_pbmm2_index { input:
    genome_fasta = genome_fasta,
    threads = pbmm2_index_threads,
    mem_gb = pbmm2_index_mem_gb,
    runtime_attributes = runtime_attributes
  }

  scatter (sample_index in range(length(group_flnc_bams_by_sm.sample_names))) {
    String sample_prefix = group_flnc_bams_by_sm.sample_prefixes[sample_index]
    Array[Int] group_bam_indices = group_flnc_bams_by_sm.group_bam_indices[sample_index]
    Int group_size = length(group_bam_indices)

    scatter (group_position in range(group_size)) {
      Int group_bam_index = group_bam_indices[group_position]
      File group_flnc_bam = flnc_bams[group_bam_index]
      String alignment_prefix = if group_size == 1
        then sample_prefix
        else "~{sample_prefix}.part-~{group_position}"

      call FlncAlignmentTasks.pbmm2_align_flnc { input:
        output_prefix = alignment_prefix,
        flnc_bam = group_flnc_bam,
        pbmm2_index = create_pbmm2_index.pbmm2_index,
        threads = pbmm2_align_threads,
        mem_gb = pbmm2_align_mem_gb,
        runtime_attributes = runtime_attributes
      }
    }

    if (group_size > 1) {
      call FlncAlignmentTasks.pbsamoa_merge_aligned_bams { input:
        bams = pbmm2_align_flnc.aligned_bam,
        out_prefix = sample_prefix,
        compression = pbsamoa_merge_compression,
        threads = pbsamoa_merge_threads,
        mem_gb = pbsamoa_merge_mem_gb,
        runtime_attributes = runtime_attributes
      }
    }

    File sample_aligned_bam = if group_size == 1
      then pbmm2_align_flnc.aligned_bam[0]
      else select_first([
        pbsamoa_merge_aligned_bams.merged_bam
      ])
    File sample_aligned_bam_index = if group_size == 1
      then pbmm2_align_flnc.aligned_bam_index[0]
      else select_first([
        pbsamoa_merge_aligned_bams.merged_bam_index
      ])
  }

  call IsocallCore.isocall_core as isocall_core { input:
    aligned_bams = sample_aligned_bam,
    aligned_bam_bais = sample_aligned_bam_index,
    isocall_extra_merged_profile = isocall_extra_merged_profile,
    annotation_gtf_gz = annotation_gtf_gz,
    genome_fasta = genome_fasta,
    genome_fasta_index = genome_fasta_index,
    runtime_attributes = runtime_attributes,
    isocall_min_read_fraction = isocall_min_read_fraction,
    isocall_max_bundles_per_gene = isocall_max_bundles_per_gene,
    isocall_min_reads_per_isoform = isocall_min_reads_per_isoform,
    isocall_profile_threads = isocall_profile_threads,
    isocall_profile_mem_gb = isocall_profile_mem_gb,
    isocall_prep_isoforms_threads = isocall_prep_isoforms_threads,
    isocall_prep_isoforms_mem_gb = isocall_prep_isoforms_mem_gb,
    isocall_merge_profiles_threads = isocall_merge_profiles_threads,
    isocall_merge_profiles_mem_gb = isocall_merge_profiles_mem_gb,
    isocall_call_threads = isocall_call_threads,
    isocall_call_mem_gb = isocall_call_mem_gb,
    output_prefix = output_prefix
  }

  call IsoformClassificationCore.isoform_classification_core as isoform_classification_core { input:
    isoforms_gtf = isocall_core.isocall_isoforms_gtf,
    isocall_count_matrix = isocall_core.isocall_count_matrix,
    pigeon_resources = prepare_pigeon_resources.pigeon_resources,
    pigeon_use_polya = pigeon_use_polya,
    pigeon_use_cage_peak = pigeon_use_cage_peak,
    pigeon_use_junction = pigeon_use_junction,
    runtime_attributes = runtime_attributes,
    pigeon_min_ref_length = pigeon_min_ref_length,
    prepare_isoforms_threads = prepare_isoforms_threads,
    prepare_isoforms_mem_gb = prepare_isoforms_mem_gb,
    pigeon_classify_threads = pigeon_classify_threads,
    pigeon_classify_mem_gb = pigeon_classify_mem_gb,
    pigeon_filter_threads = pigeon_filter_threads,
    pigeon_filter_mem_gb = pigeon_filter_mem_gb,
    output_prefix = output_prefix
  }

  output {
    String workflow_name = "secondary_analysis_core"
    Array[String] sample_names = group_flnc_bams_by_sm.sample_names
    Array[String] sample_prefixes = group_flnc_bams_by_sm.sample_prefixes
    Array[Int] group_sizes = group_flnc_bams_by_sm.group_sizes
    Array[File] aligned_bams = sample_aligned_bam
    Array[File] aligned_bam_bais = sample_aligned_bam_index
    File isocall_isoforms_gtf = isocall_core.isocall_isoforms_gtf
    File isocall_count_matrix = isocall_core.isocall_count_matrix
    File isocall_closest_known = isocall_core.isocall_closest_known
    File pigeon_classification = isoform_classification_core.pigeon_classification
    File filtered_isoforms_gtf = isoform_classification_core.filtered_isoforms_gtf
    File pigeon_filtered_classification = isoform_classification_core.pigeon_filtered_classification
  }
}
