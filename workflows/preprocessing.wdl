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
    ref_name: {
      description: "Packaged reference to use when reference_container is omitted",
      choices: [
        "GRCh38_gencode49"
      ]
    }
    reference_container: {
      description: "Optional explicit immutable reference-container URI; takes precedence over ref_name and is not rewritten by container_registry"
    }
    reference_overrides: {
      description: "Typed optional overrides for reference files"
    }
    segmentation_adapter_set: {
      description: "Optional segmentation adapter set used to select the packaged segmentation adapter FASTA when no custom override is supplied; omission selects 8-fold",
      choices: [
        "8-fold",
        "12-fold",
        "16-fold"
      ]
    }
    isoseq_primers_set: {
      description: "Optional Iso-Seq primer set used to select the packaged indexed-primer FASTA when no custom override is supplied; omission selects IsoSeq-v2",
      choices: [
        "IsoSeq-v2",
        "IsoSeq96"
      ]
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
      description: "Optional PacBio registry for registry-relative task images and named reference containers; if omitted, quay.io/pacbio is used"
    }
  }

  input {
    Array[HiFiSource] hifi_sources
    File biosample_csv
    String ref_name = "GRCh38_gencode49"
    String? reference_container
    ReferenceOverrides reference_overrides = object {
    }
    String? segmentation_adapter_set
    String? isoseq_primers_set
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
    ref_name = ref_name,
    reference_container = reference_container,
    reference_overrides = reference_overrides,
    segmentation_adapter_set = segmentation_adapter_set,
    isoseq_primers_set = isoseq_primers_set,
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
    segmentation_adapters = select_first([
      reference_resources.segmentation_adapters
    ]),
    indexed_primers = select_first([
      reference_resources.indexed_primers
    ]),
    runtime_attributes = default_runtime_attributes,
    isoseq_require_polya = isoseq_require_polya
  }

  output {
    String workflow_name = "preprocessing"
    String workflow_version = "0.4.0"
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
  }
}
