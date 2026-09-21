version 1.0

import "../rna_structs.wdl"

task unpack_reference_container {
  meta {
    description: "Unpack the typed Kinnex Iso-Seq static-resource manifest from a reference container."
    outputs: {
      resource_bundle_version: {
        description: "Version of the packaged Kinnex Iso-Seq resource bundle"
      },
      genome_fasta: {
        description: "Reference genome FASTA"
      },
      genome_fasta_index: {
        description: "Reference FASTA index"
      },
      annotation_gtf_gz: {
        description: "Compressed annotation GTF"
      },
      pigeon_poly_a: {
        description: "Pigeon polyA motif list"
      },
      pigeon_cage_peak_bed: {
        description: "Pigeon CAGE peak BED"
      },
      pigeon_junction_coverage: {
        description: "Pigeon junction coverage table"
      },
      hifi_demux_barcodes: {
        description: "Kinnex HiFi demultiplexing barcode FASTA"
      },
      indexed_primers_isoseq_v2: {
        description: "Iso-Seq v2 indexed-primer FASTA"
      },
      indexed_primers_isoseq96: {
        description: "IsoSeq96 indexed-primer FASTA"
      },
      segmentation_adapters_8fold: {
        description: "Eight-fold Kinnex adapter FASTA"
      },
      segmentation_adapters_12fold: {
        description: "Twelve-fold Kinnex adapter FASTA"
      },
      segmentation_adapters_16fold: {
        description: "Sixteen-fold Kinnex adapter FASTA"
      },
      manifest_json: {
        description: "Manifest JSON describing the packaged static workflow inputs"
      }
    }
  }

  parameter_meta {
    reference_container: {
      description: "Complete immutable reference-container URI including an @sha256 digest"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    String reference_container
    RuntimeAttributes runtime_attributes
  }

  Int threads = 1
  Int mem_gb = 4
  Int total_mem_mb = (mem_gb * 1024) + runtime_attributes.add_memory_mb

  String unpack_script = "/opt/scripts/unpack_container.py"
  String container_manifest = "/opt/manifests/manifest.json"

  command <<<
    set -euo pipefail

    python3 "~{unpack_script}" \
      --manifest "~{container_manifest}" \
      --output-dir "."
  >>>

  output {
    String resource_bundle_version = read_string("resource_bundle_version")
    File genome_fasta = glob("genome_fasta/*")[0]
    File genome_fasta_index = glob("genome_fasta_index/*")[0]
    File annotation_gtf_gz = glob("annotation_gtf_gz/*")[0]
    File pigeon_poly_a = glob("pigeon_poly_a/*")[0]
    File pigeon_cage_peak_bed = glob("pigeon_cage_peak_bed/*")[0]
    File pigeon_junction_coverage = glob("pigeon_junction_coverage/*")[0]
    File hifi_demux_barcodes = glob("hifi_demux_barcodes/*")[0]
    File indexed_primers_isoseq_v2 = glob("indexed_primers_isoseq_v2/*")[0]
    File indexed_primers_isoseq96 = glob("indexed_primers_isoseq96/*")[0]
    File segmentation_adapters_8fold = glob("skera_adapters_8fold/*")[0]
    File segmentation_adapters_12fold = glob("skera_adapters_12fold/*")[0]
    File segmentation_adapters_16fold = glob("skera_adapters_16fold/*")[0]
    File manifest_json = "manifest.json"
  }

  runtime {
    cpu: threads
    memory: total_mem_mb + " MB"
    docker: "~{reference_container}"
    maxRetries: runtime_attributes.max_retries
  }
}
