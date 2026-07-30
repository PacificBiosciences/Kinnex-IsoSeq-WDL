version 1.0

import "tasks/demux_setup.wdl" as DemuxSetupTasks
import "tasks/lima.wdl" as Lima

workflow hifi_demux {
  meta {
    description: "Run upstream HiFi demux for one preprocessing source dataset and return downstream preprocessing datasets. HiFi demux uses barcode-file naming, while the 3-column biosample CSV is used to derive the per-outer cDNA biosample CSVs."
    outputs: {
      normalized_datasets: {
        description: "Normalized outer-barcode HiFi datasets"
      },
      demuxed_hifi_datasets: {
        description: "HiFi-demuxed ConsensusReadSet XMLs"
      },
      demuxed_hifi_bams: {
        description: "HiFi-demuxed BAMs"
      }
    }
  }

  parameter_meta {
    hifi_bam: {
      description: "Input HiFi BAM"
    }
    biosample_csv: {
      description: "Biosample CSV"
    }
    hifi_demux_barcodes: {
      description: "HiFi demux barcode FASTA"
    }
    expected_hifi_barcode_pairs: {
      description: "Expected HiFi barcode pairs"
    }
    hifi_demux_lima_threads: {
      description: "CPU threads for upstream HiFi lima demultiplexing"
    }
    hifi_demux_lima_mem_gb: {
      description: "Memory allocation in GB for upstream HiFi lima demultiplexing"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    File hifi_bam
    File biosample_csv
    File hifi_demux_barcodes
    Array[String] expected_hifi_barcode_pairs
    Int hifi_demux_lima_threads = 16
    Int hifi_demux_lima_mem_gb = 64
    RuntimeAttributes runtime_attributes
  }

  String source_dataset_label = basename(hifi_bam, ".bam")

  call DemuxSetupTasks.derive_hifi_demux_biosample_csv { input:
    three_col_csv = biosample_csv,
    runtime_attributes = runtime_attributes
  }

  call Lima.run_lima { input:
    bam = hifi_bam,
    barcode_file = hifi_demux_barcodes,
    biosample_csv = derive_hifi_demux_biosample_csv.hifi_demux_biosample_csv,
    split_named = true,
    hifi_preset = "SYMMETRIC-ADAPTERS",
    store_unbarcoded = true,
    ignore_xml_biosamples = true,
    output_missing_pairs = true,
    emit_dataset_xml = true,
    output_prefix = source_dataset_label,
    expected_barcode_pairs = expected_hifi_barcode_pairs,
    threads = hifi_demux_lima_threads,
    mem_gb = hifi_demux_lima_mem_gb,
    runtime_attributes = runtime_attributes
  }

  scatter (demuxed_hifi_bam in run_lima.demuxed_bams) {
    String demuxed_hifi_name = basename(demuxed_hifi_bam, ".bam")
    String barcode_pair = sub(demuxed_hifi_name, "^.*\\.", "")
    String normalized_dataset_name = source_dataset_label + "." + barcode_pair

    call DemuxSetupTasks.derive_cdna_biosample_csv as derive_cdna_for_demuxed { input:
      three_col_csv = biosample_csv,
      outer_barcode = barcode_pair,
      runtime_attributes = runtime_attributes
    }

    HiFiDemuxedDataset normalized_dataset = object {
      source_dataset_name: source_dataset_label,
      hifi_barcode: barcode_pair,
      dataset_name: normalized_dataset_name,
      hifi_bam: demuxed_hifi_bam,
      cdna_biosample_csv: derive_cdna_for_demuxed.cdna_biosample_csv,
      cdna_barcode_pairs: derive_cdna_for_demuxed.cdna_barcode_pairs
    }
  }

  output {
    Array[HiFiDemuxedDataset] normalized_datasets = normalized_dataset
    Array[File] demuxed_hifi_datasets = run_lima.demuxed_datasets
    Array[File] demuxed_hifi_bams = run_lima.demuxed_bams
  }
}
