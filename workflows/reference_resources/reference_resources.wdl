version 1.0

import "reference_container.wdl" as ReferenceContainer

task validate_reference_sources {
  meta {
    description: "Validate reference-container and typed-override selection before resolving workflow resources."
    outputs: {
      needs_genome: {
        description: "Whether the selected profile consumes genome and annotation resources"
      },
      needs_pigeon: {
        description: "Whether the selected profile consumes optional Pigeon support resources"
      },
      needs_preprocessing: {
        description: "Whether the selected profile consumes preprocessing resources"
      },
      needs_container: {
        description: "Whether at least one effective resource must be resolved from the reference container"
      },
      reference_mode: {
        description: "Reference selection mode: container, hybrid, or custom"
      }
    }
  }

  parameter_meta {
    resolution_profile: {
      description: "Internal workflow profile that determines the required reference resources"
    }
    reference_container: {
      description: "Optional immutable reference-container URI"
    }
    kinnex_primers_set: {
      description: "Optional Kinnex primer set used to select packaged Skera adapters when no custom Skera adapter override is supplied"
    }
    present_override_names: {
      description: "Logical names of typed reference overrides supplied by the caller"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    String resolution_profile
    String? reference_container
    String? kinnex_primers_set
    Array[String] present_override_names
    RuntimeAttributes runtime_attributes
  }

  Int threads = 1
  Int mem_gb = 1
  Int total_mem_mb = (mem_gb * 1024) + runtime_attributes.add_memory_mb
  String reference_container_value = select_first([
    reference_container,
    ""
  ])
  String kinnex_primers_set_value = select_first([
    kinnex_primers_set,
    ""
  ])
  Boolean has_kinnex_primers_set = defined(kinnex_primers_set)
  File present_override_names_file = write_lines(present_override_names)

  command <<<
    set -euo pipefail

    python3 - \
      "~{resolution_profile}" \
      "~{reference_container_value}" \
      "~{kinnex_primers_set_value}" \
      "~{true="true" false="false" has_kinnex_primers_set}" \
      "~{present_override_names_file}" <<'PY'
    import pathlib
    import re
    import sys

    profile, container, adapter_set, has_adapter_set_value, override_names_file = sys.argv[1:]
    has_adapter_set = has_adapter_set_value == 'true'
    present_overrides = {name for name in pathlib.Path(override_names_file).read_text().splitlines() if name}

    resource_groups = {
        'genome': {
            'required': {
                'genome_fasta',
                'annotation_gtf_gz',
            },
            'optional': set(),
        },
        'pigeon': {
            'required': set(),
            'optional': {
                'pigeon_poly_a',
                'pigeon_cage_peak_bed',
                'pigeon_junction_coverage',
            },
        },
        'preprocessing': {
            'required': {
                'hifi_demux_barcodes',
                'barcoded_primers',
                'skera_adapters',
            },
            'optional': set(),
        },
    }
    profile_groups = {
        'kinnex_isoseq': {
            'genome',
            'pigeon',
            'preprocessing',
        },
        'preprocessing': {
            'preprocessing',
        },
        'secondary_analysis': {
            'genome',
            'pigeon',
        },
        'isocall': {
            'genome',
        },
        'isoform_classification': {
            'genome',
            'pigeon',
        },
    }

    all_resources = set()
    for resources in resource_groups.values():
        all_resources.update(resources['required'])
        all_resources.update(resources['optional'])

    errors = []
    unknown_overrides = sorted(present_overrides - all_resources)
    if unknown_overrides:
        errors.append('unknown reference override names: ' + ', '.join(unknown_overrides))
    if profile not in profile_groups:
        errors.append(f'unknown resolution profile: {profile}')
    if has_adapter_set and adapter_set not in {'8fold', '12fold', '16fold'}:
        errors.append('kinnex_primers_set must be one of: 8fold, 12fold, 16fold')

    if container and not re.search(r'@sha256:[0-9a-f]{64}$', container):
        errors.append('reference_container must end in a lowercase 64-character @sha256 digest')

    groups = profile_groups.get(profile, set())
    accepted = set()
    required = set()
    for group in groups:
        accepted.update(resource_groups[group]['required'])
        accepted.update(resource_groups[group]['optional'])
        required.update(resource_groups[group]['required'])

    if profile in profile_groups:
        irrelevant = sorted(present_overrides - accepted)
        if irrelevant:
            errors.append(f'overrides are not used by resolution profile {profile}: ' + ', '.join(irrelevant))

        if has_adapter_set and 'preprocessing' not in groups:
            errors.append(f'kinnex_primers_set is not used by resolution profile {profile}')
        if has_adapter_set and 'skera_adapters' in present_overrides:
            errors.append('kinnex_primers_set cannot be supplied with a custom skera_adapters override')

        if 'genome' in groups and 'genome_fasta' in present_overrides and 'annotation_gtf_gz' not in present_overrides:
            errors.append('annotation_gtf_gz override is required when genome_fasta is overridden')

    if not container and profile in profile_groups:
        missing = sorted(required - present_overrides)
        if missing:
            errors.append('reference_container is absent and required overrides are missing: ' + ', '.join(missing))

    if errors:
        raise SystemExit('Invalid reference configuration:\n- ' + '\n- '.join(errors))

    needs_genome = 'genome' in groups
    needs_pigeon = 'pigeon' in groups
    needs_preprocessing = 'preprocessing' in groups
    uses_genome_override = 'genome_fasta' in present_overrides

    container_fallbacks = set()
    if needs_genome and not uses_genome_override:
        container_fallbacks.add('genome_fasta')
        if 'annotation_gtf_gz' not in present_overrides:
            container_fallbacks.add('annotation_gtf_gz')
    if needs_pigeon and not uses_genome_override:
        container_fallbacks.update(resource_groups['pigeon']['optional'] - present_overrides)
    if needs_preprocessing:
        container_fallbacks.update(resource_groups['preprocessing']['required'] - present_overrides)

    needs_container = bool(container and container_fallbacks)
    uses_any_override = bool(present_overrides & accepted)
    reference_mode = 'hybrid' if needs_container and uses_any_override else 'container' if needs_container else 'custom'

    pathlib.Path('needs_genome.txt').write_text(str(needs_genome).lower() + '\n')
    pathlib.Path('needs_pigeon.txt').write_text(str(needs_pigeon).lower() + '\n')
    pathlib.Path('needs_preprocessing.txt').write_text(str(needs_preprocessing).lower() + '\n')
    pathlib.Path('needs_container.txt').write_text(str(needs_container).lower() + '\n')
    pathlib.Path('reference_mode.txt').write_text(reference_mode + '\n')
    PY
  >>>

  output {
    Boolean needs_genome = read_boolean("needs_genome.txt")
    Boolean needs_pigeon = read_boolean("needs_pigeon.txt")
    Boolean needs_preprocessing = read_boolean("needs_preprocessing.txt")
    Boolean needs_container = read_boolean("needs_container.txt")
    String reference_mode = read_string("reference_mode.txt")
  }

  runtime {
    cpu: threads
    memory: total_mem_mb + " MB"
    docker: runtime_attributes.container_registry + "/pb_wdl_base@sha256:03cb3c01937eccc907f8ad71c87b258581504572205fe3f31a657e318f3564ae"
    maxRetries: runtime_attributes.max_retries
  }
}

task index_reference_fasta {
  meta {
    description: "Generate a FASTA index for an overridden reference genome."
    outputs: {
      genome_fasta_index: {
        description: "Generated reference FASTA index"
      }
    }
  }

  parameter_meta {
    genome_fasta: {
      description: "Uncompressed reference FASTA override"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    File genome_fasta
    RuntimeAttributes runtime_attributes
  }

  Int threads = 1
  Int mem_gb = 2
  Int total_mem_mb = (mem_gb * 1024) + runtime_attributes.add_memory_mb
  String genome_fasta_index_basename = basename(genome_fasta) + ".fai"

  command <<<
    set -euo pipefail

    python3 - "~{genome_fasta}" <<'PY'
    import pathlib
    import sys

    genome_fasta = pathlib.Path(sys.argv[1])
    with genome_fasta.open('rb') as handle:
        if handle.read(2) == b'\x1f\x8b':
            raise SystemExit(
                'reference_overrides.genome_fasta must be an uncompressed FASTA; gzip and BGZF inputs are not supported'
            )
    PY

    samtools faidx \
      --fai-idx "~{genome_fasta_index_basename}" \
      "~{genome_fasta}"
  >>>

  output {
    File genome_fasta_index = genome_fasta_index_basename
  }

  runtime {
    cpu: threads
    memory: total_mem_mb + " MB"
    docker: runtime_attributes.container_registry + "/pb_wdl_base@sha256:03cb3c01937eccc907f8ad71c87b258581504572205fe3f31a657e318f3564ae"
    maxRetries: runtime_attributes.max_retries
  }
}

workflow resolve_reference_resources {
  meta {
    description: "Resolve typed container defaults and per-file overrides into one stage-specific effective reference set."
    outputs: {
      resources: {
        description: "Resolved typed reference resources required by the selected workflow profile"
      },
      reference_mode: {
        description: "Reference selection mode: container, hybrid, or custom"
      },
      reference_container_uri: {
        description: "Immutable reference-container URI when container defaults were used"
      },
      base_resource_bundle_version: {
        description: "Base container resource-bundle version when container defaults were used"
      }
    }
  }

  parameter_meta {
    resolution_profile: {
      description: "Internal workflow profile that determines the required reference resources"
    }
    reference_container: {
      description: "Optional immutable reference-container URI used for defaults"
    }
    reference_overrides: {
      description: "Typed optional overrides for packaged reference files; a genome override requires an annotation override and suppresses packaged Pigeon support fallbacks"
    }
    kinnex_primers_set: {
      description: "Optional Kinnex primer set used to select packaged Skera adapters when no custom Skera adapter override is supplied; omission selects 8fold"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    String resolution_profile
    String? reference_container
    ReferenceOverrides reference_overrides = object {
    }
    String? kinnex_primers_set
    RuntimeAttributes runtime_attributes
  }

  Array[String] present_override_names = flatten([
    if defined(reference_overrides.genome_fasta)
      then [
        "genome_fasta"
      ]
      else []
    ,
    if defined(reference_overrides.annotation_gtf_gz)
      then [
        "annotation_gtf_gz"
      ]
      else []
    ,
    if defined(reference_overrides.pigeon_poly_a)
      then [
        "pigeon_poly_a"
      ]
      else []
    ,
    if defined(reference_overrides.pigeon_cage_peak_bed)
      then [
        "pigeon_cage_peak_bed"
      ]
      else []
    ,
    if defined(reference_overrides.pigeon_junction_coverage)
      then [
        "pigeon_junction_coverage"
      ]
      else []
    ,
    if defined(reference_overrides.hifi_demux_barcodes)
      then [
        "hifi_demux_barcodes"
      ]
      else []
    ,
    if defined(reference_overrides.barcoded_primers)
      then [
        "barcoded_primers"
      ]
      else []
    ,
    if defined(reference_overrides.skera_adapters)
      then [
        "skera_adapters"
      ]
      else []
  ])

  call validate_reference_sources { input:
    resolution_profile = resolution_profile,
    reference_container = reference_container,
    kinnex_primers_set = kinnex_primers_set,
    present_override_names = present_override_names,
    runtime_attributes = runtime_attributes
  }

  Boolean needs_preprocessing = validate_reference_sources.needs_preprocessing
  Boolean needs_genome = validate_reference_sources.needs_genome
  Boolean needs_pigeon = validate_reference_sources.needs_pigeon
  Boolean needs_container = validate_reference_sources.needs_container

  Boolean uses_genome_override = defined(reference_overrides.genome_fasta)
  Boolean uses_skera_override = defined(reference_overrides.skera_adapters)
  String effective_kinnex_primers_set = select_first([
    kinnex_primers_set,
    "8fold"
  ])

  if (needs_container) {
    String required_reference_container = select_first([
      reference_container
    ])
    call ReferenceContainer.unpack_reference_container { input:
      reference_container = required_reference_container,
      runtime_attributes = runtime_attributes
    }
  }

  if (needs_genome && uses_genome_override) {
    File overridden_genome_fasta = select_first([
      reference_overrides.genome_fasta
    ])
    call index_reference_fasta { input:
      genome_fasta = overridden_genome_fasta,
      runtime_attributes = runtime_attributes
    }
  }

  if (needs_genome) {
    File selected_genome_fasta = if uses_genome_override
      then select_first([
        reference_overrides.genome_fasta
      ])
      else select_first([
        unpack_reference_container.genome_fasta
      ])
    File selected_genome_fasta_index = if uses_genome_override
      then select_first([
        index_reference_fasta.genome_fasta_index
      ])
      else select_first([
        unpack_reference_container.genome_fasta_index
      ])
  }

  if (needs_genome) {
    File selected_annotation_gtf_gz = if uses_genome_override
      then select_first([
        reference_overrides.annotation_gtf_gz
      ])
      else select_first([
        reference_overrides.annotation_gtf_gz,
        unpack_reference_container.annotation_gtf_gz
      ])
  }

  if (needs_pigeon) {
    File? selected_pigeon_poly_a = if uses_genome_override
      then reference_overrides.pigeon_poly_a
      else select_first([
        reference_overrides.pigeon_poly_a,
        unpack_reference_container.pigeon_poly_a
      ])
    File? selected_pigeon_cage_peak_bed = if uses_genome_override
      then reference_overrides.pigeon_cage_peak_bed
      else select_first([
        reference_overrides.pigeon_cage_peak_bed,
        unpack_reference_container.pigeon_cage_peak_bed
      ])
    File? selected_pigeon_junction_coverage = if uses_genome_override
      then reference_overrides.pigeon_junction_coverage
      else select_first([
        reference_overrides.pigeon_junction_coverage,
        unpack_reference_container.pigeon_junction_coverage
      ])
  }

  if (needs_preprocessing && !uses_skera_override && effective_kinnex_primers_set == "8fold") {
    File selected_packaged_skera_adapters_8fold = select_first([
      unpack_reference_container.skera_adapters_8fold
    ])
  }
  if (needs_preprocessing && !uses_skera_override && effective_kinnex_primers_set == "12fold") {
    File selected_packaged_skera_adapters_12fold = select_first([
      unpack_reference_container.skera_adapters_12fold
    ])
  }
  if (needs_preprocessing && !uses_skera_override && effective_kinnex_primers_set == "16fold") {
    File selected_packaged_skera_adapters_16fold = select_first([
      unpack_reference_container.skera_adapters_16fold
    ])
  }

  if (needs_preprocessing) {
    File selected_hifi_demux_barcodes = select_first([
      reference_overrides.hifi_demux_barcodes,
      unpack_reference_container.hifi_demux_barcodes
    ])
    File selected_barcoded_primers = select_first([
      reference_overrides.barcoded_primers,
      unpack_reference_container.barcoded_primers
    ])
    File selected_skera_adapters = select_first([
      reference_overrides.skera_adapters,
      selected_packaged_skera_adapters_8fold,
      selected_packaged_skera_adapters_12fold,
      selected_packaged_skera_adapters_16fold
    ])
  }

  String selected_reference_mode = validate_reference_sources.reference_mode
  ResolvedReferenceResources resolved_resources = object {
    genome_fasta: selected_genome_fasta,
    genome_fasta_index: selected_genome_fasta_index,
    annotation_gtf_gz: selected_annotation_gtf_gz,
    pigeon_poly_a: selected_pigeon_poly_a,
    pigeon_cage_peak_bed: selected_pigeon_cage_peak_bed,
    pigeon_junction_coverage: selected_pigeon_junction_coverage,
    hifi_demux_barcodes: selected_hifi_demux_barcodes,
    barcoded_primers: selected_barcoded_primers,
    skera_adapters: selected_skera_adapters
  }

  output {
    ResolvedReferenceResources resources = resolved_resources
    String reference_mode = selected_reference_mode
    String? reference_container_uri = required_reference_container
    String? base_resource_bundle_version = unpack_reference_container.resource_bundle_version
  }
}
