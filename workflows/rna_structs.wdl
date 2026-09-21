version 1.0

struct RuntimeAttributes {
  String backend
  Int max_retries
  Int add_memory_mb
  Int nproc
  String container_registry
}

struct ReferenceOverrides {
  File? genome_fasta
  File? annotation_gtf_gz
  File? pigeon_poly_a
  File? pigeon_cage_peak_bed
  File? pigeon_junction_coverage
  File? hifi_demux_barcodes
  File? indexed_primers
  File? segmentation_adapters
}

struct ResolvedReferenceResources {
  File? genome_fasta
  File? genome_fasta_index
  File? annotation_gtf_gz
  File? pigeon_poly_a
  File? pigeon_cage_peak_bed
  File? pigeon_junction_coverage
  File? hifi_demux_barcodes
  File? indexed_primers
  File? segmentation_adapters
}

struct PigeonResources {
  File annotation_gtf
  File annotation_gtf_pgi
  File genome_fasta
  File genome_fasta_index
  File? pigeon_poly_a
  File? pigeon_cage_peak_bed
  File? pigeon_cage_peak_bed_index
  File? pigeon_junction_coverage
  File? pigeon_junction_coverage_index
}

struct HiFiSource {
  File hifi_bam
  String? hifi_barcode
}

struct HiFiDemuxedDataset {
  String source_dataset_name
  String hifi_barcode
  String dataset_name
  File hifi_bam
  File cdna_biosample_csv
  Array[String] cdna_barcode_pairs
}

struct IsoSeqRefineReportRow {
  String report_json
  String bio_sample
  String source_dataset
  String hifi_barcode
  String cdna_barcode
}
