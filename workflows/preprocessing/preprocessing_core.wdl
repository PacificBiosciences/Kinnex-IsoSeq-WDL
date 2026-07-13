version 1.0

import "../rna_structs.wdl"
import "hifi_demux.wdl" as HiFiDemux
import "stage.wdl" as PreProcessing
import "tasks.wdl" as Tasks

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
      refine_summary_reports: {
        description: "isoseq refine filter summary reports"
      },
      flnc_dataset_xml: {
        description: "FLNC ConsensusReadSet XML"
      },
      flnc_child_dataset_xmls: {
        description: "Child FLNC ConsensusReadSet XMLs"
      },
      flnc_dataset_bams: {
        description: "Packaged FLNC BAMs"
      },
      flnc_dataset_bam_pbis: {
        description: "Packaged FLNC BAM PBIs"
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
    skera_adapters: {
      description: "skera adapter FASTA"
    }
    barcoded_primers: {
      description: "Barcoded primer FASTA"
    }
    consensusreadset_xmls: {
      description: "Optional ConsensusReadSet XMLs"
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
    File skera_adapters
    File barcoded_primers
    Array[File]? consensusreadset_xmls
    Boolean isoseq_require_polya = true
    Int hifi_demux_lima_threads = 16
    Int hifi_demux_lima_mem_gb = 64
    Int skera_split_threads = 8
    Int skera_split_mem_gb = 32
    Int cdna_lima_threads = 16
    Int cdna_lima_mem_gb = 64
    Int isoseq_refine_threads = 16
    Int isoseq_refine_mem_gb = 32
    RuntimeAttributes runtime_attributes
  }

  scatter (hifi_bam in hifi_bams) {
    File source_hifi_bam_for_validation = hifi_bam
    String source_dataset_name_for_validation = basename(hifi_bam, ".bam")
    String empty_hifi_bam_barcode = ""
  }

  Array[String] hifi_bam_barcodes_for_normalization = if needs_hifi_demux || length(hifi_bam_barcodes) != length(hifi_bams)
    then empty_hifi_bam_barcode
    else hifi_bam_barcodes

  call Tasks.validate_preprocessing_inputs { input:
    hifi_bams = hifi_bams,
    source_dataset_names = source_dataset_name_for_validation,
    hifi_bam_barcodes = hifi_bam_barcodes,
    needs_hifi_demux = needs_hifi_demux,
    biosample_csv = biosample_csv,
    hifi_demux_barcodes = hifi_demux_barcodes,
    barcoded_primers = barcoded_primers,
    runtime_attributes = runtime_attributes
  }

  if (defined(consensusreadset_xmls)) {
    Array[File] defined_consensusreadset_xmls_for_validation = select_first([
      consensusreadset_xmls
    ])

    call Tasks.validate_consensusreadset_xmls { input:
      consensusreadset_xmls = defined_consensusreadset_xmls_for_validation,
      source_hifi_bams = source_hifi_bam_for_validation,
      validated_biosample_csv = validate_preprocessing_inputs.validated_biosample_csv,
      runtime_attributes = runtime_attributes
    }
  }

  File validated_biosample_csv_for_preprocessing = select_first([
    validate_consensusreadset_xmls.validated_biosample_csv_out,
    validate_preprocessing_inputs.validated_biosample_csv
  ])

  scatter (i in range(length(hifi_bams))) {
    call HiFiDemux.hifi_demux { input:
      hifi_bam = hifi_bams[i],
      hifi_bam_barcode = hifi_bam_barcodes_for_normalization[i],
      needs_hifi_demux = needs_hifi_demux,
      biosample_csv = validated_biosample_csv_for_preprocessing,
      hifi_demux_barcodes = hifi_demux_barcodes,
      expected_hifi_barcode_pairs = validate_preprocessing_inputs.outer_barcode_pairs,
      hifi_demux_lima_threads = hifi_demux_lima_threads,
      hifi_demux_lima_mem_gb = hifi_demux_lima_mem_gb,
      runtime_attributes = runtime_attributes
    }
  }

  Array[OuterBarcodeHiFiDataset] normalized_datasets = flatten(hifi_demux.normalized_datasets)

  call PreProcessing.preprocessing_stage { input:
    datasets = normalized_datasets,
    skera_adapters = skera_adapters,
    barcoded_primers = barcoded_primers,
    skera_split_threads = skera_split_threads,
    skera_split_mem_gb = skera_split_mem_gb,
    cdna_lima_threads = cdna_lima_threads,
    cdna_lima_mem_gb = cdna_lima_mem_gb,
    isoseq_refine_threads = isoseq_refine_threads,
    isoseq_refine_mem_gb = isoseq_refine_mem_gb,
    runtime_attributes = runtime_attributes,
    isoseq_require_polya = isoseq_require_polya
  }

  if (defined(consensusreadset_xmls)) {
    File first_consensusreadset_xml = select_first([
      validate_consensusreadset_xmls.first_consensusreadset_xml
    ])
    File collection_contexts = select_first([
      validate_consensusreadset_xmls.collection_contexts
    ])

    call Tasks.populate_flnc_dataset_xml { input:
      consensusreadset_xml = first_consensusreadset_xml,
      collection_contexts = collection_contexts,
      validated_biosample_csv = validated_biosample_csv_for_preprocessing,
      flnc_bams = preprocessing_stage.flnc_bams,
      flnc_bam_pbis = preprocessing_stage.flnc_bam_pbis,
      runtime_attributes = runtime_attributes
    }
  }

  output {
    String workflow_name = "preprocessing_core"
    Array[String] source_dataset_names = hifi_demux.source_dataset_name
    Array[String] dataset_names = preprocessing_stage.dataset_names
    Array[File] hifi_demux_datasets = flatten(hifi_demux.demuxed_hifi_datasets)
    Array[File] hifi_demux_bams = flatten(hifi_demux.demuxed_hifi_bams)
    Array[String] flnc_names = preprocessing_stage.flnc_names
    Array[File] flnc_bams = preprocessing_stage.flnc_bams
    Array[File] flnc_bam_pbis = preprocessing_stage.flnc_bam_pbis
    Array[File] refine_summary_reports = preprocessing_stage.refine_summary_reports
    File? flnc_dataset_xml = populate_flnc_dataset_xml.flnc_dataset_xml
    Array[File]? flnc_child_dataset_xmls = populate_flnc_dataset_xml.child_flnc_dataset_xmls
    Array[File]? flnc_dataset_bams = populate_flnc_dataset_xml.packaged_flnc_bams
    Array[File]? flnc_dataset_bam_pbis = populate_flnc_dataset_xml.packaged_flnc_bam_pbis
  }
}
