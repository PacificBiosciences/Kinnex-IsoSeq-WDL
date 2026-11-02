version 1.0

import "tasks/isoseq_refine.wdl" as IsoSeqRefineTasks
import "tasks/lima.wdl" as Lima
import "tasks/read_segmentation_report.wdl" as ReadSegmentationReportTasks
import "tasks/skera.wdl" as SkeraTasks

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
      segmentation_summary_csv: {
        description: "Per-HiFi-barcode Skera segmentation metrics"
      },
      segmentation_summary_report: {
        description: "Read Segmentation metrics and per-barcode plots"
      },
      segmentation_summary_report_resources: {
        description: "PNG images and thumbnails referenced by the Read Segmentation report"
      },
      refine_summary_report: {
        description: "Iso-Seq Refine summary metrics, primer data, and read-length plots"
      },
      refine_summary_report_resources: {
        description: "PNG images and thumbnails referenced by the Iso-Seq Refine report"
      }
    }
  }

  parameter_meta {
    datasets: {
      description: "Outer-barcode HiFi datasets"
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
    Array[HiFiDemuxedDataset] datasets
    File segmentation_adapters
    File indexed_primers
    Boolean isoseq_require_polya = true
    Int skera_split_threads = 16
    Int skera_split_mem_gb = 32
    Int cdna_lima_threads = 16
    Int cdna_lima_mem_gb = 64
    Int isoseq_refine_threads = 16
    Int isoseq_refine_mem_gb = 32
    RuntimeAttributes runtime_attributes
  }

  scatter (dataset in datasets) {
    String hifi_barcode = dataset.hifi_barcode
    call SkeraTasks.skera_split_hifi { input:
      dataset_name = dataset.dataset_name,
      hifi_bam = dataset.hifi_bam,
      segmentation_adapters = segmentation_adapters,
      threads = skera_split_threads,
      mem_gb = skera_split_mem_gb,
      runtime_attributes = runtime_attributes
    }

    call Lima.run_lima as cdna_lima { input:
      bam = skera_split_hifi.segmented_bam,
      barcode_file = indexed_primers,
      biosample_csv = dataset.cdna_biosample_csv,
      isoseq = true,
      split_named = true,
      store_unbarcoded = true,
      overwrite_biosample_names = true,
      output_prefix = dataset.dataset_name,
      expected_barcode_pairs = dataset.cdna_barcode_pairs,
      threads = cdna_lima_threads,
      mem_gb = cdna_lima_mem_gb,
      runtime_attributes = runtime_attributes
    }

    Int lima_pbi_count = length(cdna_lima.report_bam_pbis)
    String source_dataset_name = dataset.source_dataset_name
    File cdna_biosample_csv = dataset.cdna_biosample_csv

    scatter (demuxed_index in range(length(cdna_lima.demuxed_bams))) {
      File demuxed_bam = cdna_lima.demuxed_bams[demuxed_index]

      call IsoSeqRefineTasks.isoseq_refine { input:
        dataset_name = dataset.dataset_name,
        source_dataset_name = dataset.source_dataset_name,
        hifi_barcode = dataset.hifi_barcode,
        cdna_barcode = cdna_lima.demuxed_barcode_pairs[demuxed_index],
        bio_sample = cdna_lima.demuxed_bio_samples[demuxed_index],
        demuxed_bam = demuxed_bam,
        indexed_primers = indexed_primers,
        isoseq_require_polya = isoseq_require_polya,
        threads = isoseq_refine_threads,
        mem_gb = isoseq_refine_mem_gb,
        runtime_attributes = runtime_attributes
      }
    }
  }
  call ReadSegmentationReportTasks.generate_read_segmentation_report { input:
    hifi_barcodes = hifi_barcode,
    skera_summary_csvs = skera_split_hifi.skera_summary_csv,
    skera_summary_jsons = skera_split_hifi.skera_summary_json,
    read_lengths_csvs = skera_split_hifi.read_lengths_csv,
    ligations_csvs = skera_split_hifi.ligations_csv,
    segmentation_adapters = segmentation_adapters,
    runtime_attributes = runtime_attributes
  }

  call IsoSeqRefineTasks.gather_isoseq_refine_reports { input:
    refine_report_rows = flatten(isoseq_refine.serialized_refine_report_row),
    flnc_bam_pbis = flatten(isoseq_refine.flnc_bam_pbi),
    lima_bam_pbis = flatten(cdna_lima.report_bam_pbis),
    lima_pbi_counts = lima_pbi_count,
    source_dataset_names = source_dataset_name,
    hifi_barcodes = hifi_barcode,
    cdna_biosample_csvs = cdna_biosample_csv,
    indexed_primers = indexed_primers,
    runtime_attributes = runtime_attributes
  }

  output {
    Array[String] dataset_names = skera_split_hifi.dataset_name_out
    Array[String] flnc_names = flatten(isoseq_refine.flnc_name_out)
    Array[File] flnc_bams = flatten(isoseq_refine.flnc_bam)
    Array[File] flnc_bam_pbis = flatten(isoseq_refine.flnc_bam_pbi)
    File refine_summary_report = gather_isoseq_refine_reports.refine_summary_report
    Array[File] refine_summary_report_resources = gather_isoseq_refine_reports.refine_summary_report_resources
    File segmentation_summary_csv = generate_read_segmentation_report.segmentation_summary_csv
    File segmentation_summary_report = generate_read_segmentation_report.segmentation_summary_report
    Array[File] segmentation_summary_report_resources = generate_read_segmentation_report.segmentation_summary_report_resources
  }
}
