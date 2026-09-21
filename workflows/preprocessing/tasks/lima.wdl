version 1.0

import "../../rna_structs.wdl"

task run_lima {
  meta {
    description: "Run lima with explicit parameterization for either HiFi demux or cDNA demultiplexing."
    outputs: {
      demuxed_bams: {
        description: "Demuxed BAMs"
      },
      demuxed_datasets: {
        description: "Demuxed ConsensusReadSet XMLs"
      },
      demuxed_barcode_pairs: {
        description: "Barcode pairs for demuxed BAMs"
      },
      demuxed_bio_samples: {
        description: "Bio Sample names for demuxed BAMs"
      },
      unbarcoded_bams: {
        description: "Unbarcoded BAMs"
      }
    }
  }

  parameter_meta {
    bam: {
      description: "Input BAM"
    }
    barcode_file: {
      description: "Barcode FASTA"
    }
    biosample_csv: {
      description: "Biosample CSV"
    }
    isoseq: {
      description: "Run lima in isoseq mode"
    }
    split_named: {
      description: "Emit split-named lima outputs"
    }
    hifi_preset: {
      description: "lima HiFi preset"
    }
    store_unbarcoded: {
      description: "Emit unbarcoded reads"
    }
    ignore_xml_biosamples: {
      description: "Ignore XML biosample metadata"
    }
    overwrite_biosample_names: {
      description: "Overwrite biosample names"
    }
    output_missing_pairs: {
      description: "Emit missing barcode pairs"
    }
    emit_dataset_xml: {
      description: "Emit ConsensusReadSet XMLs"
    }
    output_prefix: {
      description: "Output prefix"
    }
    expected_barcode_pairs: {
      description: "Expected barcode pairs"
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
    File bam
    File barcode_file
    File? biosample_csv
    Boolean isoseq = false
    Boolean split_named = true
    String? hifi_preset
    Boolean store_unbarcoded = false
    Boolean ignore_xml_biosamples = false
    Boolean overwrite_biosample_names = false
    Boolean output_missing_pairs = false
    Boolean emit_dataset_xml = false
    String output_prefix
    Array[String] expected_barcode_pairs
    Int threads = 16
    Int mem_gb = 64
    RuntimeAttributes runtime_attributes
  }

  Int total_mem_mb = (mem_gb * 1024) + runtime_attributes.add_memory_mb
  Int effective_threads = if (threads > runtime_attributes.nproc)
    then runtime_attributes.nproc
    else threads
  String output_target = if emit_dataset_xml
    then "demux/" + output_prefix + ".consensusreadset.xml"
    else "demux/" + output_prefix + ".bam"
  File expected_barcode_pairs_file = write_lines(expected_barcode_pairs)

  command <<<
    set -euo pipefail

    if [ "~{emit_dataset_xml}" == "true" ] && [ "~{isoseq}" == "true" ]; then
      echo "emit_dataset_xml is incompatible with --isoseq output mode" >&2
      exit 1
    fi

    mkdir -p demux

    biosample_csv="~{biosample_csv}"
    lima_biosample_args=()
    if [[ "~{defined(biosample_csv)}" == "true" ]]; then
      lima_biosample_args+=(--biosample-csv "${biosample_csv}")
    fi

    lima \
      ~{true="--isoseq" false="" isoseq} \
      ~{true="--split-named" false="" split_named} \
      ~{if defined(hifi_preset)
        then "--hifi-preset '" + hifi_preset + "'"
        else ""} \
      ~{true="--store-unbarcoded" false="" store_unbarcoded} \
      ~{true="--ignore-xml-biosamples" false="" ignore_xml_biosamples} \
      ~{true="--overwrite-biosample-names" false="" overwrite_biosample_names} \
      ~{true="--output-missing-pairs" false="" output_missing_pairs} \
      --num-threads ~{effective_threads} \
      --log-level INFO \
      --log-file "~{output_prefix}.lima.log" \
      "${lima_biosample_args[@]}" \
      "~{bam}" \
      "~{barcode_file}" \
      "~{output_target}"

    python3 - \
      "~{expected_barcode_pairs_file}" \
      "~{output_prefix}" \
      "~{emit_dataset_xml}" \
      "~{store_unbarcoded}" \
      "${biosample_csv}" <<'PY'
    import csv
    import os
    import sys

    (
        expected_barcode_pairs_file,
        output_prefix,
        emit_dataset_xml_raw,
        store_unbarcoded_raw,
        biosample_csv,
    ) = sys.argv[1:]
    emit_dataset_xml = emit_dataset_xml_raw == 'true'
    store_unbarcoded = store_unbarcoded_raw == 'true'


    def read_lines(path):
        with open(path, 'r', encoding='utf-8') as fh:
            return [line.rstrip('\n') for line in fh if line.rstrip('\n')]


    def fail(message):
        sys.exit(message)


    def write_manifest(path, values):
        with open(path, 'w', encoding='utf-8') as out_fh:
            out_fh.writelines(value + '\n' for value in values)


    def read_biosample_mapping(path):
        if not path:
            return {}

        expected_header = ('Barcodes', 'Bio Sample')
        mapping = {}
        with open(path, 'r', newline='', encoding='utf-8-sig') as fh:
            reader = csv.reader(fh)
            try:
                header = tuple(cell.strip() for cell in next(reader))
            except StopIteration:
                fail(f'{path} is empty')
            if header != expected_header:
                fail(f'{path} has header {header!r}; expected {expected_header!r}')

            for line_no, raw_row in enumerate(reader, start=2):
                row = [cell.strip() for cell in raw_row]
                if not row or all(cell == '' for cell in row):
                    continue
                if len(row) != 2:
                    fail(f'{path}:{line_no}: expected 2 columns, got {len(row)}')
                barcode_pair, bio_sample = row
                if not barcode_pair or not bio_sample:
                    fail(f'{path}:{line_no}: empty field in row {raw_row!r}')
                if barcode_pair in mapping:
                    fail(f'{path}:{line_no}: duplicate barcode pair {barcode_pair!r}')
                mapping[barcode_pair] = bio_sample

        if not mapping:
            fail(f'{path} has no data rows')
        return mapping


    expected_barcode_pairs = read_lines(expected_barcode_pairs_file)
    if not expected_barcode_pairs:
        fail('expected_barcode_pairs must contain at least one barcode pair')

    duplicates = sorted(
        {barcode_pair for barcode_pair in expected_barcode_pairs if expected_barcode_pairs.count(barcode_pair) > 1}
    )
    if duplicates:
        fail('expected_barcode_pairs contains duplicates: ' + ', '.join(duplicates))

    bio_sample_by_barcode = read_biosample_mapping(biosample_csv)
    if bio_sample_by_barcode:
        missing_bio_samples = sorted(set(expected_barcode_pairs) - set(bio_sample_by_barcode))
        unexpected_bio_samples = sorted(set(bio_sample_by_barcode) - set(expected_barcode_pairs))
        if missing_bio_samples or unexpected_bio_samples:
            fail(
                'biosample CSV barcode pairs do not match expected_barcode_pairs: '
                f'missing={missing_bio_samples!r}, unexpected={unexpected_bio_samples!r}'
            )

    observed_files = {
        os.path.join('demux', path) for path in os.listdir('demux') if os.path.isfile(os.path.join('demux', path))
    }
    allowed_files = set()
    demuxed_bams = []
    demuxed_datasets = []
    demuxed_barcode_pairs = []
    demuxed_bio_samples = []

    aggregate_dataset_xml = os.path.join(
        'demux',
        output_prefix + '.consensusreadset.xml',
    )
    allowed_files.add(aggregate_dataset_xml)
    allowed_files.update(
        {
            os.path.join('demux', output_prefix + '.lima.clips'),
            os.path.join('demux', output_prefix + '.lima.counts'),
            os.path.join('demux', output_prefix + '.lima.guess'),
            os.path.join('demux', output_prefix + '.lima.guess.json'),
            os.path.join('demux', output_prefix + '.lima.report'),
            os.path.join('demux', output_prefix + '.lima.summary'),
        }
    )

    for barcode_pair in expected_barcode_pairs:
        demuxed_bam = os.path.join(
            'demux',
            output_prefix + '.' + barcode_pair + '.bam',
        )
        demuxed_bam_pbi = demuxed_bam + '.pbi'
        demuxed_dataset = os.path.join(
            'demux',
            output_prefix + '.' + barcode_pair + '.consensusreadset.xml',
        )
        allowed_files.update({demuxed_bam, demuxed_bam_pbi, demuxed_dataset})

        has_bam = demuxed_bam in observed_files
        has_dataset = demuxed_dataset in observed_files
        if emit_dataset_xml and has_bam != has_dataset:
            fail(
                'lima produced incomplete BAM/XML outputs for barcode pair '
                f'{barcode_pair!r}: bam={has_bam}, dataset_xml={has_dataset}'
            )
        if has_bam:
            demuxed_bams.append(demuxed_bam)
            demuxed_barcode_pairs.append(barcode_pair)
            if bio_sample_by_barcode:
                demuxed_bio_samples.append(bio_sample_by_barcode[barcode_pair])
        if has_dataset:
            demuxed_datasets.append(demuxed_dataset)

    unbarcoded_bams = []
    unbarcoded_bam = os.path.join('demux', output_prefix + '.unbarcoded.bam')
    unbarcoded_sidecars = {
        unbarcoded_bam,
        unbarcoded_bam + '.pbi',
        os.path.join('demux', output_prefix + '.unbarcoded.consensusreadset.xml'),
        os.path.join('demux', output_prefix + '.unbarcoded.json'),
    }
    if store_unbarcoded:
        allowed_files.update(unbarcoded_sidecars)
        if unbarcoded_bam in observed_files:
            unbarcoded_bams.append(unbarcoded_bam)
    else:
        disallowed_unbarcoded = sorted(observed_files & unbarcoded_sidecars)
        if disallowed_unbarcoded:
            fail('lima produced unbarcoded outputs while store_unbarcoded=false: ' + ', '.join(disallowed_unbarcoded))

    unexpected_files = sorted(
        path for path in observed_files - allowed_files if path.endswith('.bam') or path.endswith('.consensusreadset.xml')
    )
    if unexpected_files:
        print(
            'WARNING: lima produced unexpected demux outputs; '
            'excluding them from downstream manifests: ' + ', '.join(unexpected_files),
            file=sys.stderr,
        )

    write_manifest('demuxed_bams.txt', demuxed_bams)
    write_manifest('demuxed_datasets.txt', demuxed_datasets)
    write_manifest('demuxed_barcode_pairs.txt', demuxed_barcode_pairs)
    write_manifest('demuxed_bio_samples.txt', demuxed_bio_samples)
    write_manifest('unbarcoded_bams.txt', unbarcoded_bams)
    PY
  >>>

  output {
    Array[File] demuxed_bams = read_lines("demuxed_bams.txt")
    Array[File] demuxed_datasets = read_lines("demuxed_datasets.txt")
    Array[String] demuxed_barcode_pairs = read_lines("demuxed_barcode_pairs.txt")
    Array[String] demuxed_bio_samples = read_lines("demuxed_bio_samples.txt")
    Array[File] unbarcoded_bams = read_lines("unbarcoded_bams.txt")
  }

  runtime {
    cpu: effective_threads
    memory: total_mem_mb + " MB"
    docker: runtime_attributes.container_registry + "/lima@sha256:c49a067e447fdd2a13b4ebb41cfd66d1242032233671dbe4a4fdbe3e62ac82ef"  # 26.2.1_build3
    maxRetries: runtime_attributes.max_retries
  }
}
