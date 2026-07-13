version 1.0

import "lima_tasks.wdl" as Lima
import "tasks.wdl" as Tasks

workflow preprocessing_stage {
  meta {
    description: "Preprocessing stage for PacBio Kinnex Iso-Seq: skera split, lima demultiplex, then isoseq refine."
    outputs: {
      dataset_names: {
        description: "Preprocessing dataset names"
      },
      flnc_names: {
        description: "FLNC names"
      },
      flnc_bams: {
        description: "FLNC BAMs"
      },
      flnc_bam_pbis: {
        description: "FLNC BAM PBIs"
      },
      refine_summary_reports: {
        description: "isoseq refine filter summary reports"
      }
    }
  }

  parameter_meta {
    datasets: {
      description: "Outer-barcode HiFi datasets"
    }
    skera_adapters: {
      description: "skera adapter FASTA"
    }
    barcoded_primers: {
      description: "Barcoded primer FASTA"
    }
    isoseq_require_polya: {
      description: "Pass --require-polya to isoseq refine"
    }
    skera_split_threads: {
      description: "CPU threads for skera split"
    }
    skera_split_mem_gb: {
      description: "Memory allocation in GB for skera split"
    }
    cdna_lima_threads: {
      description: "CPU threads for cDNA lima demultiplexing"
    }
    cdna_lima_mem_gb: {
      description: "Memory allocation in GB for cDNA lima demultiplexing"
    }
    isoseq_refine_threads: {
      description: "CPU threads for isoseq refine"
    }
    isoseq_refine_mem_gb: {
      description: "Memory allocation in GB for isoseq refine"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    Array[OuterBarcodeHiFiDataset] datasets
    File skera_adapters
    File barcoded_primers
    Boolean isoseq_require_polya = true
    Int skera_split_threads = 8
    Int skera_split_mem_gb = 32
    Int cdna_lima_threads = 16
    Int cdna_lima_mem_gb = 64
    Int isoseq_refine_threads = 16
    Int isoseq_refine_mem_gb = 32
    RuntimeAttributes runtime_attributes
  }

  scatter (dataset in datasets) {
    call Tasks.skera_split_hifi { input:
      dataset_name = dataset.dataset_name,
      hifi_bam = dataset.hifi_bam,
      skera_adapters = skera_adapters,
      threads = skera_split_threads,
      mem_gb = skera_split_mem_gb,
      runtime_attributes = runtime_attributes
    }

    call Lima.run_lima as cdna_lima { input:
      bam = skera_split_hifi.segmented_bam,
      barcode_file = barcoded_primers,
      biosample_csv = dataset.cdna_biosample_csv,
      isoseq = true,
      split_named = true,
      overwrite_biosample_names = true,
      output_prefix = dataset.dataset_name,
      expected_barcode_pairs = dataset.cdna_barcode_pairs,
      threads = cdna_lima_threads,
      mem_gb = cdna_lima_mem_gb,
      runtime_attributes = runtime_attributes
    }

    scatter (demuxed_bam in cdna_lima.demuxed_bams) {
      call Tasks.isoseq_refine { input:
        dataset_name = dataset.dataset_name,
        demuxed_bam = demuxed_bam,
        barcoded_primers = barcoded_primers,
        isoseq_require_polya = isoseq_require_polya,
        threads = isoseq_refine_threads,
        mem_gb = isoseq_refine_mem_gb,
        runtime_attributes = runtime_attributes
      }
    }
  }

  output {
    Array[String] dataset_names = skera_split_hifi.dataset_name_out
    Array[String] flnc_names = flatten(isoseq_refine.flnc_name_out)
    Array[File] flnc_bams = flatten(isoseq_refine.flnc_bam)
    Array[File] flnc_bam_pbis = flatten(isoseq_refine.flnc_bam_pbi)
    Array[File] refine_summary_reports = flatten(isoseq_refine.filter_summary_report)
  }
}
