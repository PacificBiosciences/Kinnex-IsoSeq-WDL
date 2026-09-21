version 1.0

import "../rna_structs.wdl"
import "hifi_demux.wdl" as HiFiDemux
import "stage.wdl" as PreProcessing
import "tasks/demux_setup.wdl" as DemuxSetupTasks
import "tasks/validation.wdl" as ValidationTasks

workflow preprocessing_core {
  meta {
    description: "PacBio Kinnex Iso-Seq preprocessing implementation: optional HiFi demux, then skera desegmentation, cDNA lima demultiplexing, and isoseq refine to FLNC."
    outputs: {
      workflow_name: {
        description: "Workflow name"
      },
      source_dataset_names: {
        description: "Source dataset names"
      },
      dataset_names: {
        description: "Preprocessing dataset names"
      },
      hifi_demux_datasets: {
        description: "HiFi-demuxed ConsensusReadSet XMLs"
      },
      hifi_demux_bams: {
        description: "HiFi-demuxed BAMs"
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
      refine_summary_report: {
        description: "Combined isoseq refine summary report"
      }
    }
  }

  parameter_meta {
    hifi_bams: {
      description: "Input HiFi BAMs"
    }
    hifi_bam_barcodes: {
      description: "HiFi BAM barcode pairs"
    }
    needs_hifi_demux: {
      description: "Run upstream HiFi demux"
    }
    biosample_csv: {
      description: "Biosample CSV"
    }
    hifi_demux_barcodes: {
      description: "HiFi demux barcode FASTA"
    }
    segmentation_adapters: {
      description: "Segmentation adapter FASTA"
    }
    indexed_primers: {
      description: "Iso-Seq indexed-primer FASTA"
    }
    isoseq_require_polya: {
      description: "Pass --require-polya to isoseq refine"
    }
    hifi_demux_lima_threads: {
      description: "CPU threads for upstream HiFi lima demultiplexing"
    }
    hifi_demux_lima_mem_gb: {
      description: "Memory allocation in GB for upstream HiFi lima demultiplexing"
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
    Array[File] hifi_bams
    Array[String] hifi_bam_barcodes = []
    Boolean needs_hifi_demux
    File biosample_csv
    File hifi_demux_barcodes
    File segmentation_adapters
    File indexed_primers
    Boolean isoseq_require_polya = true
    Int hifi_demux_lima_threads = 16
    Int hifi_demux_lima_mem_gb = 64
    Int skera_split_threads = 16
    Int skera_split_mem_gb = 32
    Int cdna_lima_threads = 16
    Int cdna_lima_mem_gb = 64
    Int isoseq_refine_threads = 16
    Int isoseq_refine_mem_gb = 32
    RuntimeAttributes runtime_attributes
  }

  scatter (hifi_bam in hifi_bams) {
    String source_dataset_name = basename(hifi_bam, ".bam")
    String empty_hifi_bam_barcode = ""
  }

  Array[String] hifi_bam_barcodes_for_normalization = if needs_hifi_demux || length(hifi_bam_barcodes) != length(hifi_bams)
    then empty_hifi_bam_barcode
    else hifi_bam_barcodes

  call ValidationTasks.validate_preprocessing_inputs { input:
    hifi_bams = hifi_bams,
    source_dataset_names = source_dataset_name,
    hifi_bam_barcodes = hifi_bam_barcodes,
    needs_hifi_demux = needs_hifi_demux,
    biosample_csv = biosample_csv,
    hifi_demux_barcodes = hifi_demux_barcodes,
    indexed_primers = indexed_primers,
    runtime_attributes = runtime_attributes
  }

  File validated_biosample_csv_for_preprocessing = validate_preprocessing_inputs.validated_biosample_csv

  if (needs_hifi_demux) {
    call HiFiDemux.hifi_demux { input:
      hifi_bam = hifi_bams[0],
      movie_name = validate_preprocessing_inputs.movie_name,
      biosample_csv = validated_biosample_csv_for_preprocessing,
      hifi_demux_barcodes = hifi_demux_barcodes,
      expected_hifi_barcode_pairs = validate_preprocessing_inputs.outer_barcode_pairs,
      hifi_demux_lima_threads = hifi_demux_lima_threads,
      hifi_demux_lima_mem_gb = hifi_demux_lima_mem_gb,
      runtime_attributes = runtime_attributes
    }
  }

  scatter (i in range(length(hifi_bams))) {
    if (!needs_hifi_demux) {
      call DemuxSetupTasks.derive_cdna_biosample_csv as derive_cdna_for_passthrough { input:
        three_col_csv = validated_biosample_csv_for_preprocessing,
        outer_barcode = hifi_bam_barcodes_for_normalization[i],
        runtime_attributes = runtime_attributes
      }

      HiFiDemuxedDataset passthrough_dataset = object {
        source_dataset_name: source_dataset_name[i],
        hifi_barcode: hifi_bam_barcodes_for_normalization[i],
        dataset_name: validate_preprocessing_inputs.movie_name + "." + hifi_bam_barcodes_for_normalization[i],
        hifi_bam: hifi_bams[i],
        cdna_biosample_csv: derive_cdna_for_passthrough.cdna_biosample_csv,
        cdna_barcode_pairs: derive_cdna_for_passthrough.cdna_barcode_pairs
      }
    }
  }

  Array[HiFiDemuxedDataset] normalized_datasets = flatten([
    flatten(select_all([
      hifi_demux.normalized_datasets
    ])),
    select_all(passthrough_dataset)
  ])

  call PreProcessing.preprocessing_stage { input:
    datasets = normalized_datasets,
    segmentation_adapters = segmentation_adapters,
    indexed_primers = indexed_primers,
    skera_split_threads = skera_split_threads,
    skera_split_mem_gb = skera_split_mem_gb,
    cdna_lima_threads = cdna_lima_threads,
    cdna_lima_mem_gb = cdna_lima_mem_gb,
    isoseq_refine_threads = isoseq_refine_threads,
    isoseq_refine_mem_gb = isoseq_refine_mem_gb,
    runtime_attributes = runtime_attributes,
    isoseq_require_polya = isoseq_require_polya
  }

  output {
    String workflow_name = "preprocessing_core"
    Array[String] source_dataset_names = source_dataset_name
    Array[String] dataset_names = preprocessing_stage.dataset_names
    Array[File] hifi_demux_datasets = flatten(select_all([
      hifi_demux.demuxed_hifi_datasets
    ]))
    Array[File] hifi_demux_bams = flatten(select_all([
      hifi_demux.demuxed_hifi_bams
    ]))
    Array[String] flnc_names = preprocessing_stage.flnc_names
    Array[File] flnc_bams = preprocessing_stage.flnc_bams
    Array[File] flnc_bam_pbis = preprocessing_stage.flnc_bam_pbis
    File refine_summary_report = preprocessing_stage.refine_summary_report
  }
}
