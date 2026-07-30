version 1.0

import "../../rna_structs.wdl"

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

task derive_hifi_demux_biosample_csv {
  meta {
    description: "Derive a lima-compatible 2-column HiFi demux biosample CSV from the 3-column biosample CSV. The output is a Barcodes,Bio Sample mapping from HiFi Barcode to Bio Sample."
    outputs: {
      hifi_demux_biosample_csv: {
        description: "HiFi demux biosample CSV"
      }
    }
  }

  parameter_meta {
    three_col_csv: {
      description: "Three-column biosample CSV"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    File three_col_csv
    RuntimeAttributes runtime_attributes
  }

  Int threads = 1
  Int mem_gb = 4
  Int total_mem_mb = (mem_gb * 1024) + runtime_attributes.add_memory_mb

  String output_basename = "hifi_demux_biosample.csv"

  command <<<
    set -euo pipefail

    python3 - "~{three_col_csv}" "~{output_basename}" <<'PY'
    import csv
    import sys

    three_col_csv, out_csv = sys.argv[1], sys.argv[2]

    expected_header = ('HiFi Barcode', 'cDNA Barcode', 'Bio Sample')
    hifi_barcodes = set()
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
            hifi_bc, _cdna_bc, _sample = row
            hifi_barcodes.add(hifi_bc)

    if not hifi_barcodes:
        sys.exit(f'{three_col_csv} has no data rows')

    with open(out_csv, 'w', newline='\n', encoding='utf-8') as fh:
        writer = csv.writer(fh, lineterminator='\n')
        writer.writerow(('Barcodes', 'Bio Sample'))
        for hifi_bc in sorted(hifi_barcodes):
            writer.writerow((hifi_bc, f'hifi_biosample_{hifi_bc}'))
    PY
  >>>

  output {
    File hifi_demux_biosample_csv = "~{output_basename}"
  }

  runtime {
    cpu: threads
    memory: total_mem_mb + " MB"
    docker: runtime_attributes.container_registry + "/pb_wdl_base@sha256:03cb3c01937eccc907f8ad71c87b258581504572205fe3f31a657e318f3564ae"
    maxRetries: runtime_attributes.max_retries
  }
}
