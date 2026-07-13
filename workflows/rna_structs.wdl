version 1.0

struct RuntimeAttributes {
  String backend
  Int max_retries
  Int add_memory_mb
  Int nproc
  String container_registry
}

struct PigeonResources {
  File annotation_gtf
  File annotation_gtf_pgi
  File genome_fasta
  File genome_fasta_index
  File pigeon_poly_a
  File pigeon_cage_peak_bed
  File pigeon_cage_peak_bed_index
  File pigeon_junction_coverage
  File pigeon_junction_coverage_index
}

struct HiFiSource {
  File hifi_bam
  String? hifi_barcode
}

# Internal-only outer-barcode-scoped HiFi dataset
struct OuterBarcodeHiFiDataset {
  String dataset_name
  File hifi_bam
  File cdna_biosample_csv
  Array[String] cdna_barcode_pairs
}
