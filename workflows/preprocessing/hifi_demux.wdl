version 1.0

import "lima_tasks.wdl" as Lima
import "tasks.wdl" as Tasks

workflow hifi_demux {
  meta {
    description: "Normalize one preprocessing source dataset by optionally running upstream HiFi demux and returning downstream preprocessing datasets. Upstream HiFi demux uses barcode-file naming, while the 3-column biosample CSV is used to derive the per-outer cDNA biosample CSVs."
    outputs: {
      source_dataset_name: {
        description: "Source dataset name"
      },
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
    hifi_bam_barcode: {
      description: "HiFi BAM barcode pair"
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
    String hifi_bam_barcode
    Boolean needs_hifi_demux
    File biosample_csv
    File hifi_demux_barcodes
    Array[String] expected_hifi_barcode_pairs
    Int hifi_demux_lima_threads = 16
    Int hifi_demux_lima_mem_gb = 64
    RuntimeAttributes runtime_attributes
  }

  String source_dataset_label = basename(hifi_bam, ".bam")

  if (needs_hifi_demux) {
    call Lima.run_lima { input:
      bam = hifi_bam,
      barcode_file = hifi_demux_barcodes,
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

      call Tasks.derive_cdna_biosample_csv as derive_cdna_for_demuxed { input:
        three_col_csv = biosample_csv,
        outer_barcode = barcode_pair,
        runtime_attributes = runtime_attributes
      }

      OuterBarcodeHiFiDataset normalized_dataset = object {
        dataset_name: normalized_dataset_name,
        hifi_bam: demuxed_hifi_bam,
        cdna_biosample_csv: derive_cdna_for_demuxed.cdna_biosample_csv,
        cdna_barcode_pairs: derive_cdna_for_demuxed.cdna_barcode_pairs
      }
    }
  }

  if (!needs_hifi_demux) {
    call Tasks.derive_cdna_biosample_csv as derive_cdna_for_passthrough { input:
      three_col_csv = biosample_csv,
      outer_barcode = hifi_bam_barcode,
      runtime_attributes = runtime_attributes
    }

    OuterBarcodeHiFiDataset passthrough_dataset = object {
      dataset_name: source_dataset_label,
      hifi_bam: hifi_bam,
      cdna_biosample_csv: derive_cdna_for_passthrough.cdna_biosample_csv,
      cdna_barcode_pairs: derive_cdna_for_passthrough.cdna_barcode_pairs
    }
  }

  output {
    String source_dataset_name = source_dataset_label
    Array[OuterBarcodeHiFiDataset] normalized_datasets = flatten([
      select_all([
        passthrough_dataset
      ]),
      flatten(select_all([
        normalized_dataset
      ]))
    ])
    Array[File] demuxed_hifi_datasets = flatten(select_all([
      run_lima.demuxed_datasets
    ]))
    Array[File] demuxed_hifi_bams = flatten(select_all([
      run_lima.demuxed_bams
    ]))
  }
}
