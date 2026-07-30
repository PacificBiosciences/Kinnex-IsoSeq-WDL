version 1.0

import "../../rna_structs.wdl"

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

    printf '%s\n' "~{sep="\" \"" hifi_bams}" > hifi_bams.txt

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

    printf '%s\n' "~{sep="\" \"" consensusreadset_xmls}" > consensusreadset_xmls.txt
    printf '%s\n' "~{sep="\" \"" source_hifi_bams}" > source_hifi_bams.txt

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
