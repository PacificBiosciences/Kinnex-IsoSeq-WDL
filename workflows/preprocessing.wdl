version 1.0

import "backend_configuration.wdl" as BackendConfiguration
import "preprocessing/preprocessing_core.wdl" as PreprocessingCore
import "reference_resources/reference_resources.wdl" as ReferenceResources

workflow preprocessing {
  meta {
    description: "PacBio Kinnex Iso-Seq preprocessing workflow for demuxed-on-instrument or non-demuxed HiFi input: optional HiFi demux, then skera desegmentation, cDNA lima demultiplexing, and isoseq refine to FLNC."
    outputs: {
      workflow_name: {
        description: "Workflow name"
      },
      workflow_version: {
        description: "Workflow version"
      },
      reference_container_uri: {
        description: "Immutable reference-container URI when container defaults were used"
      },
      reference_mode: {
        description: "Reference selection mode: container, hybrid, or custom"
      },
      base_resource_bundle_version: {
        description: "Base container resource-bundle version when container defaults were used"
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
    hifi_sources: {
      description: "Input HiFi source BAMs. Omit hifi_barcode for HiFi demux mode, or provide one hifi_barcode per BAM for cDNA demux-only mode"
    }
    biosample_csv: {
      description: "Shared 3-column biosample CSV with header `HiFi Barcode,cDNA Barcode,Bio Sample`"
    }
    reference_container: {
      description: "Optional immutable reference-container URI used for defaults"
    }
    reference_overrides: {
      description: "Typed optional overrides for reference files"
    }
    kinnex_primers_set: {
      description: "Optional Kinnex primer set used to select the packaged Skera adapter FASTA when no custom override is supplied; omission selects 8fold",
      choices: [
        "8fold",
        "12fold",
        "16fold"
      ]
    }
    consensusreadset_xmls: {
      description: "SMRT Link/internal metadata: optional input ConsensusReadSet XMLs for FLNC dataset XML generation"
    }
    isoseq_require_polya: {
      description: "Pass --require-polya to isoseq refine"
    }
    backend: {
      description: "Backend where the workflow will be executed; only HPC is supported"
    }
    max_retries: {
      description: "Maximum retries for failed task attempts"
    }
    add_memory_mb: {
      name: "Add Task Memory (MB)",
      description: "Increasing this number allocates extra memory per task when submitting jobs to the compute backend."
    }
    nproc: {
      description: "Maximum CPU threads per task for development and backend throttling",
      hidden: true
    }
    container_registry: {
      description: "Optional PacBio registry for registry-relative task images; if omitted, quay.io/pacbio is used"
    }
  }

  input {
    Array[HiFiSource] hifi_sources
    File biosample_csv
    String? reference_container
    ReferenceOverrides reference_overrides = object {
    }
    String? kinnex_primers_set
    Array[File]? consensusreadset_xmls
    Boolean isoseq_require_polya = true

    # Backend configuration
    String backend = "HPC"
    Int max_retries = 2
    Int add_memory_mb = 0
    Int nproc = 32
    String? container_registry
  }

  call BackendConfiguration.backend_configuration { input:
    backend = backend,
    max_retries = max_retries,
    add_memory_mb = add_memory_mb,
    nproc = nproc,
    container_registry = container_registry
  }

  RuntimeAttributes default_runtime_attributes = backend_configuration.runtime_attributes

  call ReferenceResources.resolve_reference_resources { input:
    resolution_profile = "preprocessing",
    reference_container = reference_container,
    reference_overrides = reference_overrides,
    kinnex_primers_set = kinnex_primers_set,
    runtime_attributes = default_runtime_attributes
  }

  ResolvedReferenceResources reference_resources = resolve_reference_resources.resources

  scatter (hifi_source in hifi_sources) {
    File hifi_bam = hifi_source.hifi_bam
    String? optional_hifi_barcode = hifi_source.hifi_barcode
  }

  Array[String] hifi_bam_barcodes = select_all(optional_hifi_barcode)
  Boolean needs_hifi_demux = length(hifi_bam_barcodes) == 0

  call PreprocessingCore.preprocessing_core as preprocessing_core { input:
    hifi_bams = hifi_bam,
    needs_hifi_demux = needs_hifi_demux,
    biosample_csv = biosample_csv,
    hifi_demux_barcodes = select_first([
      reference_resources.hifi_demux_barcodes
    ]),
    hifi_bam_barcodes = hifi_bam_barcodes,
    skera_adapters = select_first([
      reference_resources.skera_adapters
    ]),
    barcoded_primers = select_first([
      reference_resources.barcoded_primers
    ]),
    consensusreadset_xmls = consensusreadset_xmls,
    runtime_attributes = default_runtime_attributes,
    isoseq_require_polya = isoseq_require_polya
  }

  output {
    String workflow_name = "preprocessing"
    String workflow_version = "0.3.0"
    String? reference_container_uri = resolve_reference_resources.reference_container_uri
    String reference_mode = resolve_reference_resources.reference_mode
    String? base_resource_bundle_version = resolve_reference_resources.base_resource_bundle_version
    Array[String] source_dataset_names = preprocessing_core.source_dataset_names
    Array[String] dataset_names = preprocessing_core.dataset_names
    Array[File] hifi_demux_datasets = preprocessing_core.hifi_demux_datasets
    Array[File] hifi_demux_bams = preprocessing_core.hifi_demux_bams
    Array[String] flnc_names = preprocessing_core.flnc_names
    Array[File] flnc_bams = preprocessing_core.flnc_bams
    Array[File] flnc_bam_pbis = preprocessing_core.flnc_bam_pbis
    File refine_summary_report = preprocessing_core.refine_summary_report
    File? flnc_dataset_xml = preprocessing_core.flnc_dataset_xml
    Array[File]? flnc_child_dataset_xmls = preprocessing_core.flnc_child_dataset_xmls
    Array[File]? flnc_dataset_bams = preprocessing_core.flnc_dataset_bams
    Array[File]? flnc_dataset_bam_pbis = preprocessing_core.flnc_dataset_bam_pbis
  }
}
