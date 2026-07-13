version 1.0

import "../rna_structs.wdl"

task validate_preprocessing_inputs {
  meta {
    description: "Validate preprocessing input shape and barcode/sample-sheet compatibility before any demux or skera work starts."
    outputs: {
      validation_report: {
        description: "Preprocessing input validation report"
      },
      validated_biosample_csv: {
        description: "Validated biosample CSV"
      },
      outer_barcode_pairs: {
        description: "Outer barcode pairs"
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
    source_dataset_names: {
      description: "Source dataset names"
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
    barcoded_primers: {
      description: "Barcoded primer FASTA"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    Array[File] hifi_bams
    Array[String] hifi_bam_barcodes
    Array[String] source_dataset_names
    Boolean needs_hifi_demux
    File biosample_csv
    File hifi_demux_barcodes
    File barcoded_primers
    RuntimeAttributes runtime_attributes
  }

  Int threads = 1
  Int mem_gb = 4
  Int total_mem_mb = (mem_gb * 1024) + runtime_attributes.add_memory_mb

  File hifi_bam_barcodes_file = write_lines(hifi_bam_barcodes)
  File source_dataset_names_file = write_lines(source_dataset_names)

  command <<<
    set -euo pipefail

    cat > hifi_bams.txt <<'EOF'
    ~{sep="\n" hifi_bams}
    EOF

    needs_hifi_demux="~{true="true" false="false" needs_hifi_demux}"

    python3 - \
      "${needs_hifi_demux}" \
      "hifi_bams.txt" \
      "~{source_dataset_names_file}" \
      "~{hifi_bam_barcodes_file}" \
      "~{biosample_csv}" \
      "~{hifi_demux_barcodes}" \
      "~{barcoded_primers}" \
      > preprocessing_input_validation.txt <<'PY'
    import csv
    import re
    import sys
    from collections import Counter, defaultdict

    import pysam

    (
        needs_hifi_demux_raw,
        hifi_bams_path,
        source_names_path,
        hifi_bam_barcodes_path,
        biosample_csv,
        hifi_demux_barcodes,
        barcoded_primers,
    ) = sys.argv[1:]

    needs_hifi_demux = needs_hifi_demux_raw == 'true'
    bio_sample_name_re = re.compile(r'^[A-Za-z0-9_-]{1,40}$')
    pysam.set_verbosity(0)


    def read_lines(path):
        with open(path, 'r', encoding='utf-8') as fh:
            return [line.rstrip('\n') for line in fh if line.rstrip('\n')]


    def fail(message):
        sys.exit(message)


    def duplicates(values):
        counts = Counter(values)
        return sorted(value for value, count in counts.items() if count > 1)


    def split_barcode_pair(value, label):
        parts = value.split('--')
        if len(parts) != 2 or not parts[0] or not parts[1]:
            fail(f'{label} {value!r} must be an exact barcode pair such as left--right')
        return parts


    def read_fasta_names(path, label):
        names = []
        with open(path, 'r', encoding='utf-8') as fh:
            for line_no, line in enumerate(fh, start=1):
                if line.startswith('>'):
                    stripped = line[1:].strip()
                    if not stripped:
                        fail(f'{label}: FASTA record header on line {line_no} is empty')
                    names.append(stripped.split()[0])
        duplicate_names = duplicates(names)
        if duplicate_names:
            fail(f'duplicate FASTA record names in {label}: ' + ', '.join(duplicate_names))
        return set(names)


    def validate_bio_sample_name(sample_name, label):
        if not bio_sample_name_re.fullmatch(sample_name):
            fail(
                f'{label}: Bio Sample value {sample_name!r} '
                'must be 40 characters or fewer and contain only alphanumeric '
                'characters (A-Z, a-z, 0-9), underscores (_), and hyphens (-)'
            )


    def read_single_pu_value(bam):
        try:
            with pysam.AlignmentFile(bam, 'rb', check_sq=False) as bam_fh:
                read_groups = bam_fh.header.to_dict().get('RG', [])
        except (OSError, ValueError) as exc:
            fail(f'failed to read BAM header for {bam}: {exc}')

        observed_pu = set()
        for read_group in read_groups:
            if 'PU' in read_group:
                observed_pu.add(read_group['PU'])

        if not read_groups:
            fail(f'{bam} has no @RG records')
        if not observed_pu:
            fail(f'{bam} has no @RG PU tags')
        if len(observed_pu) != 1:
            fail(f'{bam} has multiple @RG PU values: ' + ', '.join(repr(pu) for pu in sorted(observed_pu)))

        pu = next(iter(observed_pu))
        if any(char in pu for char in '\t\r\n'):
            fail(f'{bam} has a PU value with tab, carriage return, or newline characters: {pu!r}')
        return pu


    hifi_bams = read_lines(hifi_bams_path)
    source_names = read_lines(source_names_path)
    hifi_bam_barcodes = read_lines(hifi_bam_barcodes_path)
    if not hifi_bams:
        fail('hifi_sources must contain at least one source BAM')
    if len(source_names) != len(hifi_bams):
        fail('internal error: source name and hifi_bams arrays differ in length')

    duplicate_source_names = duplicates(source_names)
    if duplicate_source_names:
        fail('duplicate source BAM basenames are not allowed: ' + ', '.join(duplicate_source_names))

    if needs_hifi_demux:
        if len(hifi_bams) != 1:
            fail('hifi_sources must contain exactly one source BAM in HiFi demux mode')
        if hifi_bam_barcodes:
            fail('hifi_sources must omit hifi_barcode in HiFi demux mode')
    else:
        if len(hifi_bam_barcodes) != len(hifi_bams):
            fail(
                'each hifi_sources entry must include hifi_barcode in cDNA demux-only mode: '
                f'got {len(hifi_bam_barcodes)} hifi_barcode values for {len(hifi_bams)} source BAMs'
            )
        bam_pu_values = [(bam, read_single_pu_value(bam)) for bam in hifi_bams]
        distinct_pu_values = sorted({pu for _, pu in bam_pu_values})
        if len(distinct_pu_values) != 1:
            fail(
                'all hifi_sources BAMs must have the same @RG PU value in '
                'cDNA demux-only mode; observed ' + ', '.join(f'{bam}: {pu!r}' for bam, pu in bam_pu_values)
            )

    groups = defaultdict(list)
    outer_barcode_order = []
    expected_header = ('HiFi Barcode', 'cDNA Barcode', 'Bio Sample')
    with open(biosample_csv, 'r', newline='', encoding='utf-8-sig') as fh:
        reader = csv.reader(fh)
        try:
            header = tuple(cell.strip() for cell in next(reader))
        except StopIteration:
            fail(f'{biosample_csv} is empty')
        if header != expected_header:
            fail(f'{biosample_csv} has header {header!r}; expected {expected_header!r}')

        seen_cdna_by_outer = defaultdict(set)
        for line_no, raw_row in enumerate(reader, start=2):
            row = [cell.strip() for cell in raw_row]
            if not row or all(cell == '' for cell in row):
                continue
            if len(row) != 3:
                fail(f'{biosample_csv}:{line_no}: expected 3 columns, got {len(row)}')
            outer, cdna, sample = row
            if not outer or not cdna or not sample:
                fail(f'{biosample_csv}:{line_no}: empty field in row {raw_row!r}')
            outer_parts = split_barcode_pair(outer, f'{biosample_csv}:{line_no}: outer barcode')
            if needs_hifi_demux and outer_parts[0] != outer_parts[1]:
                fail(f'{biosample_csv}:{line_no}: outer barcode {outer!r} must be symmetric in HiFi demux mode')
            split_barcode_pair(cdna, f'{biosample_csv}:{line_no}: cDNA Barcode')
            if cdna in seen_cdna_by_outer[outer]:
                fail(f'{biosample_csv}:{line_no}: duplicate cDNA Barcode {cdna!r} under outer barcode {outer!r}')
            seen_cdna_by_outer[outer].add(cdna)
            validate_bio_sample_name(sample, f'{biosample_csv}:{line_no}')
            if outer not in groups:
                outer_barcode_order.append(outer)
            groups[outer].append((cdna, sample))

    if not groups:
        fail(f'{biosample_csv} has no data rows')

    cdna_names = read_fasta_names(barcoded_primers, 'preprocessing.barcoded_primers')
    hifi_names = read_fasta_names(
        hifi_demux_barcodes,
        'preprocessing.hifi_demux_barcodes',
    )

    for hifi_bam_barcode in hifi_bam_barcodes:
        split_barcode_pair(hifi_bam_barcode, 'hifi_sources.hifi_barcode')
        if hifi_bam_barcode not in groups:
            fail(f'hifi_sources hifi_barcode value {hifi_bam_barcode!r} not found in {biosample_csv}')

    missing_cdna_components = []
    for cdna in sorted({cdna for rows in groups.values() for cdna, _ in rows}):
        for component in split_barcode_pair(cdna, 'cDNA Barcode'):
            if component not in cdna_names:
                missing_cdna_components.append(component)
    if missing_cdna_components:
        fail('cDNA Barcode components are absent from barcoded_primers: ' + ', '.join(sorted(set(missing_cdna_components))))

    hifi_barcodes_to_check = sorted(groups) if needs_hifi_demux else hifi_bam_barcodes
    missing_hifi_components = []
    for outer in hifi_barcodes_to_check:
        for component in split_barcode_pair(outer, 'HiFi Barcode'):
            if component not in hifi_names:
                missing_hifi_components.append(component)
    if missing_hifi_components:
        fail(
            'HiFi Barcode components are absent from hifi_demux_barcodes: '
            + ', '.join(sorted(set(missing_hifi_components)))
        )

    print('preprocessing input validation passed')
    print(f'source_datasets={len(source_names)}')
    print(f'biosample_outer_barcodes={len(groups)}')

    with open('outer_barcode_pairs.txt', 'w', encoding='utf-8') as out_fh:
        for outer in outer_barcode_order:
            out_fh.write(outer + '\n')
    PY
  >>>

  output {
    File validation_report = "preprocessing_input_validation.txt"
    File validated_biosample_csv = biosample_csv
    Array[String] outer_barcode_pairs = read_lines("outer_barcode_pairs.txt")
  }

  runtime {
    cpu: threads
    memory: total_mem_mb + " MB"
    docker: runtime_attributes.container_registry + "/pb_wdl_base@sha256:03cb3c01937eccc907f8ad71c87b258581504572205fe3f31a657e318f3564ae"
    maxRetries: runtime_attributes.max_retries
  }
}

task validate_consensusreadset_xmls {
  meta {
    description: "Validate SMRT Link ConsensusReadSet XML metadata and source BAM acquisition compatibility."
    outputs: {
      validation_report: {
        description: "ConsensusReadSet XML validation report"
      },
      first_consensusreadset_xml: {
        description: "First ConsensusReadSet XML"
      },
      collection_contexts: {
        description: "Collection contexts"
      },
      validated_biosample_csv_out: {
        description: "Validated biosample CSV"
      }
    }
  }

  parameter_meta {
    consensusreadset_xmls: {
      description: "ConsensusReadSet XMLs"
    }
    source_hifi_bams: {
      description: "Source HiFi BAMs"
    }
    validated_biosample_csv: {
      description: "Validated biosample CSV"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    Array[File] consensusreadset_xmls
    Array[File] source_hifi_bams
    File validated_biosample_csv
    RuntimeAttributes runtime_attributes
  }

  Int threads = 1
  Int mem_gb = 4
  Int total_mem_mb = (mem_gb * 1024) + runtime_attributes.add_memory_mb

  command <<<
    set -euo pipefail

    cat > consensusreadset_xmls.txt <<'EOF'
    ~{sep="\n" consensusreadset_xmls}
    EOF

    cat > source_hifi_bams.txt <<'EOF'
    ~{sep="\n" source_hifi_bams}
    EOF

    python3 - \
      "consensusreadset_xmls.txt" \
      "source_hifi_bams.txt" \
      > consensusreadset_xml_validation.txt <<'PY'
    import sys

    from pbcore.io import BamReader, ConsensusReadSet

    consensusreadset_xmls_file, source_hifi_bams_file = sys.argv[1:]


    def fail(message):
        sys.exit(message)


    def read_lines(path):
        with open(path, 'r', encoding='utf-8') as fh:
            return [line.rstrip('\n') for line in fh if line.rstrip('\n')]


    def get_required_nonnegative_metadata(metadata, xml_path, element_name):
        if element_name not in metadata.tags:
            fail(f'{xml_path} is missing {element_name} metadata')
        value = metadata.getMemberV(element_name, default=None, asType=int)
        if value is None:
            fail(f'{xml_path} {element_name} metadata must be a parseable nonnegative integer')
        if value < 0:
            fail(f'{xml_path} {element_name} metadata must be nonnegative; observed {value}')
        return value


    def validate_xml(xml_path):
        try:
            dataset = ConsensusReadSet(xml_path, skipCounts=True, skipMissing=True)
        except Exception as exc:
            fail(f'{xml_path} is not a valid ConsensusReadSet XML: {exc}')

        get_required_nonnegative_metadata(dataset.metadata, xml_path, 'TotalLength')
        get_required_nonnegative_metadata(dataset.metadata, xml_path, 'NumRecords')

        collections = list(dataset.metadata.collections)
        if not collections:
            fail(f'{xml_path} must contain collection metadata')

        collection_ids = set()
        contexts = set()
        consensus_read_set_ref_ids = set()
        for index, collection in enumerate(collections, start=1):
            unique_id = (collection.uniqueId or '').strip()
            context = (collection.context or '').strip()
            consensus_read_set_ref_id = (collection.consensusReadSetRef.uuid or '').strip()
            if not unique_id:
                fail(f'{xml_path} CollectionMetadata[{index}] is missing UniqueId')
            if not context:
                fail(f'{xml_path} CollectionMetadata[{index}] is missing Context')
            if not consensus_read_set_ref_id:
                fail(f'{xml_path} CollectionMetadata[{index}] is missing ConsensusReadSetRef UniqueId')
            collection_ids.add(unique_id)
            contexts.add(context)
            consensus_read_set_ref_ids.add(consensus_read_set_ref_id)

        # TODO: context validation: should all XML collection contexts collapse to one?
        return collection_ids, contexts, consensus_read_set_ref_ids


    def read_bam_acquisition(bam):
        try:
            with BamReader(bam) as reader:
                read_groups = list(reader.readGroupTable)
        except Exception as exc:
            fail(f'pbcore could not read source HiFi BAM header for {bam}: {exc}')

        platform_units = set()
        for index, read_group in enumerate(read_groups, start=1):
            platform_unit = str(read_group.MovieName or '').strip()
            if not platform_unit:
                read_group_id = str(read_group.StringID or f'record {index}')
                fail(f'{bam} @RG {read_group_id!r} is missing PU acquisition/platform-unit')
            platform_units.add(platform_unit)

        if not read_groups:
            fail(f'{bam} has no @RG records')
        if len(platform_units) != 1:
            fail(
                f'{bam} must have exactly one distinct @RG PU acquisition/platform-unit; '
                'observed ' + ', '.join(repr(value) for value in sorted(platform_units))
            )
        return next(iter(platform_units))


    xmls = read_lines(consensusreadset_xmls_file)
    source_bams = read_lines(source_hifi_bams_file)

    if not xmls:
        fail('preprocessing.consensusreadset_xmls is defined but empty')
    if not source_bams:
        fail('hifi_sources must contain at least one source BAM')

    all_collection_ids = set()
    all_contexts = set()
    all_consensus_read_set_ref_ids = set()
    for xml_path in xmls:
        collection_ids, contexts, consensus_read_set_ref_ids = validate_xml(xml_path)
        all_collection_ids.update(collection_ids)
        all_contexts.update(contexts)
        all_consensus_read_set_ref_ids.update(consensus_read_set_ref_ids)

    # TODO: collection metadata selection: should full collection metadata equality be required?
    if len(all_collection_ids) != 1:
        fail(
            'all consensusreadset_xmls must share one CollectionMetadata UniqueId; observed '
            + ', '.join(repr(value) for value in sorted(all_collection_ids))
        )
    if len(all_consensus_read_set_ref_ids) != 1:
        fail(
            'all consensusreadset_xmls must share one CollectionMetadata '
            'ConsensusReadSetRef UniqueId; observed '
            + ', '.join(repr(value) for value in sorted(all_consensus_read_set_ref_ids))
        )

    for bam in source_bams:
        acquisition = read_bam_acquisition(bam)
        if acquisition not in all_contexts:
            fail(
                f'{bam} @RG PU acquisition/platform-unit {acquisition!r} does not match '
                'any consensusreadset_xmls CollectionMetadata Context; observed contexts: '
                + ', '.join(repr(value) for value in sorted(all_contexts))
            )

    print('consensusreadset XML validation passed')
    print(f'consensusreadset_xmls={len(xmls)}')
    print(f'source_bams={len(source_bams)}')
    print(f'collection_unique_id={next(iter(all_collection_ids))}')
    print('consensusreadset_ref_unique_id=' + next(iter(all_consensus_read_set_ref_ids)))
    print('collection_contexts=' + ','.join(sorted(all_contexts)))

    with open('collection_contexts.txt', 'w', encoding='utf-8') as out_fh:
        out_fh.writelines(context + '\n' for context in sorted(all_contexts))

    with open('first.consensusreadset.xml', 'wb') as out_fh, open(xmls[0], 'rb') as in_fh:
        out_fh.write(in_fh.read())
    PY
  >>>

  output {
    File validation_report = "consensusreadset_xml_validation.txt"
    File first_consensusreadset_xml = "first.consensusreadset.xml"
    File collection_contexts = "collection_contexts.txt"
    File validated_biosample_csv_out = validated_biosample_csv
  }

  runtime {
    cpu: threads
    memory: total_mem_mb + " MB"
    docker: runtime_attributes.container_registry + "/pbcore@sha256:47fc6f1174605be9a8a932f17c15d1a29fc5e6b89ea6809f32995318c7afe1e7"  # 2.6.0_build3
    maxRetries: runtime_attributes.max_retries
  }
}

task derive_cdna_biosample_csv {
  meta {
    description: "Derive a lima-compatible 2-column cDNA biosample CSV for a single outer barcode by filtering the 3-column biosample CSV. The output is a Barcodes,Bio Sample mapping from cDNA Barcode to Bio Sample scoped to the given outer_barcode."
    outputs: {
      cdna_biosample_csv: {
        description: "cDNA biosample CSV"
      },
      cdna_barcode_pairs: {
        description: "cDNA barcode pairs"
      }
    }
  }

  parameter_meta {
    three_col_csv: {
      description: "Three-column biosample CSV"
    }
    outer_barcode: {
      description: "Outer barcode pair"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    File three_col_csv
    String outer_barcode
    RuntimeAttributes runtime_attributes
  }

  Int threads = 1
  Int mem_gb = 4
  Int total_mem_mb = (mem_gb * 1024) + runtime_attributes.add_memory_mb

  String output_basename = outer_barcode + ".cdna_biosample.csv"

  command <<<
    set -euo pipefail

    python3 - "~{three_col_csv}" "~{outer_barcode}" "~{output_basename}" <<'PY'
    import csv
    import sys

    three_col_csv, outer_barcode, out_csv = sys.argv[1], sys.argv[2], sys.argv[3]

    expected_header = ('HiFi Barcode', 'cDNA Barcode', 'Bio Sample')
    rows = []
    with open(three_col_csv, 'r', newline='', encoding='utf-8-sig') as fh:
        reader = csv.reader(fh)
        try:
            header = tuple(cell.strip() for cell in next(reader))
        except StopIteration:
            sys.exit(f'{three_col_csv} is empty')
        if header != expected_header:
            sys.exit(f'{three_col_csv} has header {header!r}; expected {expected_header!r}')

        for line_no, raw_row in enumerate(reader, start=2):
            row = [cell.strip() for cell in raw_row]
            if not row or all(cell == '' for cell in row):
                continue
            if len(row) != 3:
                sys.exit(f'{three_col_csv}:{line_no}: expected 3 columns, got {len(row)}')
            outer, cdna, sample = row
            if outer != outer_barcode:
                continue
            rows.append((cdna, sample))

    if not rows:
        sys.exit(f'outer barcode {outer_barcode!r} not found in {three_col_csv}')

    with open(out_csv, 'w', newline='\n', encoding='utf-8') as fh:
        writer = csv.writer(fh, lineterminator='\n')
        writer.writerow(('Barcodes', 'Bio Sample'))
        for cdna, sample in rows:
            writer.writerow((cdna, sample))

    with open('cdna_barcode_pairs.txt', 'w', encoding='utf-8') as out_fh:
        for cdna, _sample in rows:
            out_fh.write(cdna + '\n')
    PY
  >>>

  output {
    File cdna_biosample_csv = "~{output_basename}"
    Array[String] cdna_barcode_pairs = read_lines("cdna_barcode_pairs.txt")
  }

  runtime {
    cpu: threads
    memory: total_mem_mb + " MB"
    docker: runtime_attributes.container_registry + "/pb_wdl_base@sha256:03cb3c01937eccc907f8ad71c87b258581504572205fe3f31a657e318f3564ae"
    maxRetries: runtime_attributes.max_retries
  }
}

task skera_split_hifi {
  meta {
    description: "Split one HiFi BAM into S-reads with skera."
    outputs: {
      dataset_name_out: {
        description: "Dataset name"
      },
      segmented_bam: {
        description: "Segmented BAM"
      }
    }
  }

  parameter_meta {
    dataset_name: {
      description: "Dataset name"
    }
    hifi_bam: {
      description: "HiFi BAM"
    }
    skera_adapters: {
      description: "skera adapter FASTA"
    }
    threads: {
      description: "CPU threads"
    }
    mem_gb: {
      description: "Memory allocation in GB"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    String dataset_name
    File hifi_bam
    File skera_adapters
    Int threads = 8
    Int mem_gb = 32
    RuntimeAttributes runtime_attributes
  }

  Int total_mem_mb = (mem_gb * 1024) + runtime_attributes.add_memory_mb
  Int effective_threads = if (threads > runtime_attributes.nproc)
    then runtime_attributes.nproc
    else threads

  String output_prefix = dataset_name

  command <<<
    set -euo pipefail

    skera split \
      --num-threads ~{effective_threads} \
      --log-level INFO \
      --log-file "~{output_prefix}.skera.log" \
      "~{hifi_bam}" \
      "~{skera_adapters}" \
      "~{output_prefix}.sreads.bam"
  >>>

  output {
    String dataset_name_out = dataset_name
    File segmented_bam = "~{output_prefix}.sreads.bam"
  }

  runtime {
    cpu: effective_threads
    memory: total_mem_mb + " MB"
    docker: runtime_attributes.container_registry + "/skera@sha256:f823bd64f2beec351cae82e67c9c6f257b4896ddd59f8b2050568d59f3a165f6"  # 1.4.0_build3
    maxRetries: runtime_attributes.max_retries
  }
}

task isoseq_refine {
  meta {
    description: "Run isoseq refine on one demultiplexed BAM to remove concatemers and trim polyA tails."
    outputs: {
      dataset_name_out: {
        description: "Dataset name"
      },
      flnc_name_out: {
        description: "FLNC name"
      },
      flnc_bam: {
        description: "FLNC BAM"
      },
      flnc_bam_pbi: {
        description: "FLNC BAM PBI"
      },
      filter_summary_report: {
        description: "Filter summary report"
      },
      refine_report: {
        description: "isoseq refine report"
      }
    }
  }

  parameter_meta {
    dataset_name: {
      description: "Dataset name"
    }
    demuxed_bam: {
      description: "Demuxed BAM"
    }
    barcoded_primers: {
      description: "Barcoded primer FASTA"
    }
    isoseq_require_polya: {
      description: "Pass --require-polya to isoseq refine"
    }
    threads: {
      description: "CPU threads"
    }
    mem_gb: {
      description: "Memory allocation in GB"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    String dataset_name
    File demuxed_bam
    File barcoded_primers
    Boolean isoseq_require_polya = true
    Int threads = 16
    Int mem_gb = 32
    RuntimeAttributes runtime_attributes
  }

  Int total_mem_mb = (mem_gb * 1024) + runtime_attributes.add_memory_mb
  Int effective_threads = if (threads > runtime_attributes.nproc)
    then runtime_attributes.nproc
    else threads

  String demuxed_name = basename(demuxed_bam, ".bam")
  String output_prefix = demuxed_name + ".flnc"

  command <<<
    set -euo pipefail

    isoseq refine \
      --num-threads ~{effective_threads} \
      --log-level INFO \
      --log-file "~{demuxed_name}.isoseq_refine.log" \
      ~{true="--require-polya" false="" isoseq_require_polya} \
      "~{demuxed_bam}" \
      "~{barcoded_primers}" \
      "~{output_prefix}.bam"
  >>>

  output {
    String dataset_name_out = dataset_name
    String flnc_name_out = output_prefix
    File flnc_bam = "~{output_prefix}.bam"
    File flnc_bam_pbi = "~{output_prefix}.bam.pbi"
    File filter_summary_report = "~{output_prefix}.filter_summary.report.json"
    File refine_report = "~{output_prefix}.report.csv"
  }

  runtime {
    cpu: effective_threads
    memory: total_mem_mb + " MB"
    docker: runtime_attributes.container_registry + "/isoseq@sha256:d5e4bfd5dd570b26b9769b7b340beeefd618324685662ce2fc898b9a817e465e"  # 26.2.0_build2
    maxRetries: runtime_attributes.max_retries
  }
}

task populate_flnc_dataset_xml {
  meta {
    description: "Package final FLNC BAM/PBI files and create a top-level FLNC ConsensusReadSet XML from a ConsensusReadSet XML metadata template."
    outputs: {
      flnc_dataset_xml: {
        description: "FLNC ConsensusReadSet XML"
      },
      child_flnc_dataset_xmls: {
        description: "Child FLNC ConsensusReadSet XMLs"
      },
      packaged_flnc_bams: {
        description: "Packaged FLNC BAMs"
      },
      packaged_flnc_bam_pbis: {
        description: "Packaged FLNC BAM PBIs"
      }
    }
  }

  parameter_meta {
    consensusreadset_xml: {
      description: "ConsensusReadSet XML template"
    }
    collection_contexts: {
      description: "Collection contexts"
    }
    validated_biosample_csv: {
      description: "Validated biosample CSV"
    }
    flnc_bams: {
      description: "FLNC BAMs"
    }
    flnc_bam_pbis: {
      description: "FLNC BAM PBIs"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    File consensusreadset_xml
    File collection_contexts
    File validated_biosample_csv
    Array[File] flnc_bams
    Array[File] flnc_bam_pbis
    RuntimeAttributes runtime_attributes
  }

  Int threads = 1
  Int mem_gb = 4
  Int total_mem_mb = (mem_gb * 1024) + runtime_attributes.add_memory_mb

  command <<<
    set -euo pipefail

    cat > flnc_bams.txt <<'EOF'
    ~{sep="\n" flnc_bams}
    EOF

    cat > flnc_bam_pbis.txt <<'EOF'
    ~{sep="\n" flnc_bam_pbis}
    EOF

    python3 - \
      "~{consensusreadset_xml}" \
      "~{collection_contexts}" \
      "~{validated_biosample_csv}" \
      "flnc_bams.txt" \
      "flnc_bam_pbis.txt" <<'PY'
    import copy
    import csv
    import errno
    import os
    import re
    import shutil
    import sys

    from pbcore.io import BamReader, ConsensusReadSet
    from pbcore.io.dataset.DataSetMembers import (
        BioSampleMetadata,
        BioSamplesMetadata,
        DNABarcode,
        DNABarcodes,
        ExternalResource,
    )

    (
        consensusreadset_xml,
        collection_contexts_file,
        validated_biosample_csv,
        flnc_bams_file,
        flnc_bam_pbis_file,
    ) = sys.argv[1:]
    dataset_xml = 'flnc.consensusreadset.xml'

    with open(collection_contexts_file, 'r', encoding='utf-8') as fh:
        collection_contexts = {line.rstrip('\n') for line in fh if line.rstrip('\n')}
    with open(flnc_bams_file, 'r', encoding='utf-8') as fh:
        flnc_bams = [line.rstrip('\n') for line in fh if line.rstrip('\n')]
    with open(flnc_bam_pbis_file, 'r', encoding='utf-8') as fh:
        flnc_bam_pbis = [line.rstrip('\n') for line in fh if line.rstrip('\n')]

    if not collection_contexts:
        sys.exit('validated consensusreadset XML collection context set is empty')
    if not flnc_bams:
        sys.exit('flnc_bams must contain at least one FLNC BAM')
    if len(flnc_bams) != len(flnc_bam_pbis):
        sys.exit(
            f'flnc_bams and flnc_bam_pbis must contain the same number of files: {len(flnc_bams)} != {len(flnc_bam_pbis)}'
        )


    def exactly_one(values, what, bam_path):
        if len(values) != 1:
            observed = ', '.join(repr(value) for value in sorted(values)) or 'none'
            sys.exit(f'{bam_path} must have exactly one distinct {what}; observed {observed}')
        return next(iter(values))


    def read_required_bam_sample_and_platform_unit(bam_path):
        with BamReader(bam_path) as reader:
            read_groups = list(reader.readGroupTable)
        if not read_groups:
            sys.exit(f'{bam_path} has no @RG records')

        sample = exactly_one(
            {str(rg.SampleName).strip() for rg in read_groups if str(rg.SampleName).strip()},
            '@RG SM sample',
            bam_path,
        )
        platform_unit = exactly_one(
            {str(rg.MovieName).strip() for rg in read_groups if str(rg.MovieName).strip()},
            '@RG PU acquisition/platform-unit',
            bam_path,
        )
        if platform_unit not in collection_contexts:
            sys.exit(
                f'{bam_path} @RG PU acquisition/platform-unit {platform_unit!r} does not match '
                'any validated consensusreadset_xmls CollectionMetadata Context; observed contexts: '
                + ', '.join(repr(value) for value in sorted(collection_contexts))
            )
        return sample, platform_unit


    def read_biosample_csv(path):
        expected_header = ('HiFi Barcode', 'cDNA Barcode', 'Bio Sample')
        cdna_barcodes_by_sample = {}
        with open(path, 'r', newline='', encoding='utf-8-sig') as fh:
            reader = csv.reader(fh)
            try:
                header = tuple(cell.strip() for cell in next(reader))
            except StopIteration:
                sys.exit(f'{path} is empty')
            if header != expected_header:
                sys.exit(f'{path} has header {header!r}; expected {expected_header!r}')

            for line_no, raw_row in enumerate(reader, start=2):
                row = [cell.strip() for cell in raw_row]
                if not row or all(cell == '' for cell in row):
                    continue
                if len(row) != 3:
                    sys.exit(f'{path}:{line_no}: expected 3 columns, got {len(row)}')
                _outer_barcode, cdna_barcode, sample = row
                if not cdna_barcode or not sample:
                    sys.exit(f'{path}:{line_no}: empty cDNA barcode or Bio Sample field')
                cdna_barcodes_by_sample.setdefault(sample, [])
                if cdna_barcode not in cdna_barcodes_by_sample[sample]:
                    cdna_barcodes_by_sample[sample].append(cdna_barcode)

        if not cdna_barcodes_by_sample:
            sys.exit(f'{path} has no data rows')
        return cdna_barcodes_by_sample


    def iter_biosamples(dataset):
        for bio_sample in dataset.metadata.bioSamples:
            yield bio_sample
        for collection in dataset.metadata.collections:
            for bio_sample in collection.wellSample.bioSamples:
                yield bio_sample


    def inherited_barcode_uuids(dataset):
        uuids_by_barcode = {}
        for bio_sample in iter_biosamples(dataset):
            for dna_barcode in bio_sample.DNABarcodes:
                barcode_name = (dna_barcode.name or '').strip()
                barcode_uuid = (dna_barcode.uniqueId or '').strip()
                if barcode_name and barcode_uuid:
                    uuids_by_barcode.setdefault(barcode_name, set()).add(barcode_uuid)
        return {
            barcode_name: next(iter(barcode_uuids))
            for barcode_name, barcode_uuids in uuids_by_barcode.items()
            if len(barcode_uuids) == 1
        }


    def required_consensus_read_set_uuid(dataset, xml_path):
        consensus_read_set_uuids = set()
        for index, collection in enumerate(dataset.metadata.collections, start=1):
            consensus_read_set_uuid = (collection.consensusReadSetRef.uuid or '').strip()
            if not consensus_read_set_uuid:
                sys.exit(f'{xml_path} CollectionMetadata[{index}] is missing ConsensusReadSetRef UniqueId')
            consensus_read_set_uuids.add(consensus_read_set_uuid)
        if len(consensus_read_set_uuids) != 1:
            sys.exit(
                f'{xml_path} must contain exactly one ConsensusReadSetRef UniqueId; observed '
                + ', '.join(repr(value) for value in sorted(consensus_read_set_uuids))
            )
        return next(iter(consensus_read_set_uuids))


    def link_or_copy(src, dest):
        try:
            os.link(src, dest)
        except OSError as exc:
            if exc.errno in (errno.EXDEV, errno.EEXIST, errno.EPERM):
                shutil.copy2(src, dest)
            else:
                raise


    def candidate_bam_name(path, index, attempt):
        basename = os.path.basename(path)
        if basename.endswith('.bam'):
            stem = basename.removesuffix('.bam')
            suffix = '.bam'
        else:
            stem, suffix = os.path.splitext(basename)
        marker = f'{index:03d}' if attempt == 1 else f'{index:03d}-{attempt}'
        return f'{stem}.{marker}{suffix}'


    def make_packaged_bam_name(path, index, colliding_basenames, used_names):
        basename = os.path.basename(path)
        if basename not in colliding_basenames and basename not in used_names:
            return basename
        attempt = 1
        while True:
            candidate = candidate_bam_name(path, index, attempt)
            if candidate not in used_names:
                return candidate
            attempt += 1


    def make_biosamples(samples, cdna_barcodes_by_sample, uuid_by_barcode):
        bio_samples = BioSamplesMetadata()
        for sample in samples:
            bio_sample = BioSampleMetadata()
            bio_sample.name = sample
            dna_barcodes = DNABarcodes()
            for cdna_barcode in cdna_barcodes_by_sample[sample]:
                dna_barcode = DNABarcode()
                dna_barcode.name = cdna_barcode
                barcode_uuid = uuid_by_barcode.get(cdna_barcode)
                if barcode_uuid:
                    dna_barcode.uniqueId = barcode_uuid
                dna_barcodes.append(dna_barcode)
            bio_sample.append(dna_barcodes)
            bio_samples.append(bio_sample)
        return bio_samples


    def replace_biosample_metadata(
        dataset,
        samples,
        cdna_barcodes_by_sample,
        uuid_by_barcode,
        consensus_read_set_uuid,
    ):
        # The consensus_read_set_uuid is the instrument ConsensusReadSetRef UUID and
        # must stay distinct from the generated dataset.uuid
        dataset.metadata.removeChildren('BioSamples')
        dataset.metadata.append(make_biosamples(samples, cdna_barcodes_by_sample, uuid_by_barcode))
        for collection in dataset.metadata.collections:
            collection.consensusReadSetRef.uuid = consensus_read_set_uuid
            well_sample = collection.wellSample
            well_sample.removeChildren('BioSamples')
            well_sample.append(make_biosamples(samples, cdna_barcodes_by_sample, uuid_by_barcode))


    def sanitize_xml_name(value):
        cleaned = re.sub(r'[^A-Za-z0-9._-]+', '_', value.strip())
        cleaned = cleaned.strip('._-')
        return cleaned or 'biosample'


    def make_child_xml_name(sample, sample_index):
        # sample_index is unique per produced biosample, so this name is unique.
        return f'{sample_index:03d}.{sanitize_xml_name(sample)}.flnc.consensusreadset.xml'


    def add_child_dataset_resource(parent_dataset, parent_bam, child_dataset, child_xml):
        parent_bam_keys = {
            parent_bam,
            os.path.abspath(parent_bam),
        }
        parent_resource = None
        for ext_resource in parent_dataset.externalResources:
            resource_keys = {
                ext_resource.resourceId,
                os.path.abspath(ext_resource.resourceId),
                ext_resource.attrib.get('ResourceId', ''),
            }
            if parent_bam_keys & resource_keys:
                parent_resource = ext_resource
                break
        if parent_resource is None:
            sys.exit(f'could not find parent FLNC BAM ExternalResource {parent_bam!r} while linking child dataset XML')

        child_resource = ExternalResource()
        child_resource.resourceId = child_xml
        child_resource.uniqueId = child_dataset.uuid
        child_resource.metaType = child_dataset.datasetType
        parent_resource._setSubResByMetaType(child_dataset.datasetType, child_resource)


    basename_counts = {}
    for bam in flnc_bams:
        basename_counts[os.path.basename(bam)] = basename_counts.get(os.path.basename(bam), 0) + 1
    colliding_basenames = {basename for basename, count in basename_counts.items() if count > 1}

    packaged_bams = []
    sample_order = []
    bams_by_sample = {}
    used_packaged_names = set()
    for index, (bam, pbi) in enumerate(zip(flnc_bams, flnc_bam_pbis, strict=True), start=1):
        packaged_bam_name = make_packaged_bam_name(
            bam,
            index,
            colliding_basenames,
            used_packaged_names,
        )
        used_packaged_names.add(packaged_bam_name)

        packaged_bam = packaged_bam_name
        packaged_pbi = packaged_bam + '.pbi'
        link_or_copy(bam, packaged_bam)
        link_or_copy(pbi, packaged_pbi)

        sample, _platform_unit = read_required_bam_sample_and_platform_unit(packaged_bam)

        packaged_bams.append(packaged_bam)
        if sample not in bams_by_sample:
            sample_order.append(sample)
            bams_by_sample[sample] = []
        bams_by_sample[sample].append(packaged_bam)

    ds_template = ConsensusReadSet(consensusreadset_xml, skipCounts=True, skipMissing=True)
    cdna_barcodes_by_sample = read_biosample_csv(validated_biosample_csv)
    missing_samples = [sample for sample in sample_order if sample not in cdna_barcodes_by_sample]
    if missing_samples:
        sys.exit(
            'final FLNC BAM @RG SM sample(s) are absent from the preprocessing biosample CSV: '
            + ', '.join(repr(sample) for sample in missing_samples)
        )
    uuid_by_barcode = inherited_barcode_uuids(ds_template)
    instrument_consensus_read_set_uuid = required_consensus_read_set_uuid(
        ds_template,
        consensusreadset_xml,
    )

    template_tags = [tag.strip() for tag in str(ds_template.tags or '').split(',') if tag.strip()]
    if 'flnc' not in {tag.lower() for tag in template_tags}:
        template_tags.append('flnc')
    flnc_tags = ','.join(template_tags)


    def build_flnc_dataset(dataset, samples, name):
        dataset.updateCounts()
        dataset.metadata.collections = copy.deepcopy(ds_template.metadata.collections)
        dataset.name = name
        dataset.tags = flnc_tags
        dataset.newUuid(random=True)
        replace_biosample_metadata(
            dataset,
            samples,
            cdna_barcodes_by_sample,
            uuid_by_barcode,
            consensus_read_set_uuid=instrument_consensus_read_set_uuid,
        )


    ds_out = ConsensusReadSet(*packaged_bams)
    build_flnc_dataset(ds_out, sample_order, ds_template.name)
    top_level_uuid = ds_out.uuid

    child_xmls = []
    for sample_index, sample in enumerate(sample_order, start=1):
        child_ds = ConsensusReadSet(*bams_by_sample[sample])
        build_flnc_dataset(child_ds, [sample], f'{ds_template.name} {sample} FLNC')
        child_ds.metadata.addParentDataSet(
            top_level_uuid,
            ds_out.datasetType,
            ds_out.timeStampedName or ds_out.name or dataset_xml,
        )
        child_xml = make_child_xml_name(sample, sample_index)
        child_xmls.append(child_xml)
        child_ds.write(child_xml, relPaths=True)
        add_child_dataset_resource(ds_out, bams_by_sample[sample][0], child_ds, child_xml)

    ds_out.write(dataset_xml, relPaths=True)


    def write_manifest(path, values):
        with open(path, 'w', encoding='utf-8') as out_fh:
            out_fh.writelines(value + '\n' for value in values)


    write_manifest('child_flnc_dataset_xmls.txt', child_xmls)
    write_manifest('packaged_flnc_bams.txt', packaged_bams)
    write_manifest('packaged_flnc_bam_pbis.txt', [bam + '.pbi' for bam in packaged_bams])
    PY
  >>>

  output {
    File flnc_dataset_xml = "flnc.consensusreadset.xml"
    Array[File] child_flnc_dataset_xmls = read_lines("child_flnc_dataset_xmls.txt")
    Array[File] packaged_flnc_bams = read_lines("packaged_flnc_bams.txt")
    Array[File] packaged_flnc_bam_pbis = read_lines("packaged_flnc_bam_pbis.txt")
  }

  runtime {
    cpu: threads
    memory: total_mem_mb + " MB"
    docker: runtime_attributes.container_registry + "/pbcore@sha256:47fc6f1174605be9a8a932f17c15d1a29fc5e6b89ea6809f32995318c7afe1e7"  # 2.6.0_build3
    maxRetries: runtime_attributes.max_retries
  }
}
