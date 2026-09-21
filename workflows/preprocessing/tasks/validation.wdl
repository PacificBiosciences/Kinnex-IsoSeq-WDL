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
      },
      movie_name: {
        description: "Acquisition movie name from the BAM header PU value"
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
      description: "Source dataset names used for reporting and filename-barcode validation"
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
    indexed_primers: {
      description: "Iso-Seq indexed-primer FASTA"
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
    File indexed_primers
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
      "~{indexed_primers}" \
      > preprocessing_input_validation.txt <<'PY'
    import csv
    import datetime
    import json
    import re
    import sys
    import uuid
    from collections import Counter, defaultdict

    import pysam

    (
        needs_hifi_demux_raw,
        hifi_bams_path,
        source_names_path,
        hifi_bam_barcodes_path,
        biosample_csv,
        hifi_demux_barcodes,
        indexed_primers,
    ) = sys.argv[1:]

    needs_hifi_demux = needs_hifi_demux_raw == 'true'
    bio_sample_name_re = re.compile(r'^[A-Za-z0-9_-]{1,40}$')
    movie_name_re = re.compile(r'^[A-Za-z0-9_-]+$')
    pysam.set_verbosity(0)


    def read_lines(path):
        with open(path, 'r', encoding='utf-8') as fh:
            return [line.rstrip('\n') for line in fh if line.rstrip('\n')]


    def fail(message):
        # SL server will automatically detect and load alarms.json
        with open('alarms.json', 'wt') as alarms_json:
            alarms_json.write(
                json.dumps(
                    [
                        {
                            'exception': 'ValidationError',
                            'info': message,
                            'message': message,
                            'name': 'Validation Error',
                            'severity': 'ERROR',
                            'owner': 'validate_preprocessing_inputs',
                            'createdAt': datetime.datetime.now().strftime('%Y-%m-%dT%H:%M:%S'),
                            'id': str(uuid.uuid4()),
                        }
                    ]
                )
            )
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


    def filename_hifi_barcode_conflicts(source_name, provided_barcode, hifi_names):
        provided_components = set(split_barcode_pair(provided_barcode, 'hifi_sources.hifi_barcode'))
        conflicts = []
        for field in source_name.split('.'):
            field_parts = field.split('--')
            is_known_pair = (
                len(field_parts) == 2 and all(field_parts) and all(component in hifi_names for component in field_parts)
            )
            if is_known_pair:
                if field != provided_barcode:
                    conflicts.append(field)
            elif field in hifi_names and field not in provided_components:
                conflicts.append(field)
        return list(dict.fromkeys(conflicts))


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
        if not movie_name_re.fullmatch(pu):
            fail(
                f'{bam} has invalid @RG PU movie name {pu!r}; '
                'expected one or more ASCII letters, digits, underscores, or hyphens'
            )
        return pu


    hifi_bams = read_lines(hifi_bams_path)
    source_names = read_lines(source_names_path)
    hifi_bam_barcodes = read_lines(hifi_bam_barcodes_path)
    if not hifi_bams:
        fail('hifi_sources must contain at least one source BAM')
    if len(source_names) != len(hifi_bams):
        fail('internal error: source name and hifi_bams arrays differ in length')

    bam_pu_values = [(bam, read_single_pu_value(bam)) for bam in hifi_bams]
    distinct_pu_values = sorted({pu for _, pu in bam_pu_values})
    if len(distinct_pu_values) != 1:
        fail(
            'all hifi_sources BAMs must have the same @RG PU value; observed '
            + ', '.join(f'{bam}: {pu!r}' for bam, pu in bam_pu_values)
        )
    movie_name = distinct_pu_values[0]

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

    cdna_names = read_fasta_names(indexed_primers, 'preprocessing.indexed_primers')
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
        fail('cDNA Barcode components are absent from indexed_primers: ' + ', '.join(sorted(set(missing_cdna_components))))

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

    if not needs_hifi_demux:
        for source_name, hifi_bam_barcode in zip(source_names, hifi_bam_barcodes):
            conflicts = filename_hifi_barcode_conflicts(source_name, hifi_bam_barcode, hifi_names)
            if conflicts:
                conflict_text = ', '.join(repr(conflict) for conflict in conflicts)
                print(
                    f'WARNING: source BAM basename {source_name!r} contains HiFi barcode '
                    f'value(s) {conflict_text} that do not correspond to provided '
                    f'hifi_barcode {hifi_bam_barcode!r}',
                    file=sys.stderr,
                )

    print('preprocessing input validation passed')
    print(f'source_datasets={len(source_names)}')
    print(f'biosample_outer_barcodes={len(groups)}')
    print(f'movie_name={movie_name}')

    with open('outer_barcode_pairs.txt', 'w', encoding='utf-8') as out_fh:
        for outer in outer_barcode_order:
            out_fh.write(outer + '\n')

    with open('movie_name.txt', 'w', encoding='utf-8') as out_fh:
        out_fh.write(movie_name + '\n')
    PY
  >>>

  output {
    File validation_report = "preprocessing_input_validation.txt"
    File validated_biosample_csv = biosample_csv
    Array[String] outer_barcode_pairs = read_lines("outer_barcode_pairs.txt")
    String movie_name = read_string("movie_name.txt")
  }

  runtime {
    cpu: threads
    memory: total_mem_mb + " MB"
    docker: runtime_attributes.container_registry + "/pb_wdl_base@sha256:03cb3c01937eccc907f8ad71c87b258581504572205fe3f31a657e318f3564ae"
    maxRetries: runtime_attributes.max_retries
  }
}
