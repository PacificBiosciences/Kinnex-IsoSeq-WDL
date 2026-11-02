version 1.0

import "../../rna_structs.wdl"

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
      },
      refine_report_row: {
        description: "Typed isoseq refine report row"
      },
      serialized_refine_report_row: {
        description: "Serialized isoseq refine report row"
      }
    }
  }

  parameter_meta {
    dataset_name: {
      description: "Dataset name"
    }
    source_dataset_name: {
      description: "Source dataset name"
    }
    hifi_barcode: {
      description: "HiFi barcode pair"
    }
    cdna_barcode: {
      description: "cDNA barcode pair"
    }
    bio_sample: {
      description: "Bio Sample name"
    }
    demuxed_bam: {
      description: "Demuxed BAM"
    }
    indexed_primers: {
      description: "Iso-Seq indexed-primer FASTA"
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
    String source_dataset_name
    String hifi_barcode
    String cdna_barcode
    String bio_sample
    File demuxed_bam
    File indexed_primers
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
      "~{indexed_primers}" \
      "~{output_prefix}.bam"
  >>>

  output {
    String dataset_name_out = dataset_name
    String flnc_name_out = output_prefix
    File flnc_bam = "~{output_prefix}.bam"
    File flnc_bam_pbi = "~{output_prefix}.bam.pbi"
    File filter_summary_report = "~{output_prefix}.filter_summary.report.json"
    File refine_report = "~{output_prefix}.report.csv"
    IsoSeqRefineReportRow refine_report_row = object {
      report_json: read_string(filter_summary_report),
      bio_sample: bio_sample,
      source_dataset: source_dataset_name,
      hifi_barcode: hifi_barcode,
      cdna_barcode: cdna_barcode
    }
    File serialized_refine_report_row = write_json(refine_report_row)
  }

  runtime {
    cpu: effective_threads
    memory: total_mem_mb + " MB"
    docker: runtime_attributes.container_registry + "/isoseq@sha256:d5e4bfd5dd570b26b9769b7b340beeefd618324685662ce2fc898b9a817e465e"  # 26.2.0_build2
    maxRetries: runtime_attributes.max_retries
  }
}

task gather_isoseq_refine_reports {
  meta {
    description: "Gather partition-aware Iso-Seq Refine metrics, primer data, and read-length plots."
    outputs: {
      refine_summary_report: {
        description: "Iso-Seq Refine summary metrics, primer data, and read-length plots"
      },
      refine_summary_report_resources: {
        description: "PNG images and thumbnails referenced by the Iso-Seq Refine report"
      }
    }
  }

  parameter_meta {
    refine_report_rows: {
      description: "Serialized per-partition isoseq refine report rows"
    }
    flnc_bam_pbis: {
      description: "FLNC indexes in matching refine-row order"
    }
    lima_bam_pbis: {
      description: "All assigned and unbarcoded cDNA Lima indexes"
    }
    lima_pbi_counts: {
      description: "Number of Lima indexes per outer dataset"
    }
    source_dataset_names: {
      description: "Source dataset identity per outer dataset"
    }
    hifi_barcodes: {
      description: "HiFi barcode pair per outer dataset"
    }
    cdna_biosample_csvs: {
      description: "cDNA barcode-to-Bio Sample mappings per outer dataset"
    }
    indexed_primers: {
      description: "Iso-Seq indexed-primer FASTA in PBI index order"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    Array[File] refine_report_rows
    Array[File] flnc_bam_pbis
    Array[File] lima_bam_pbis
    Array[Int] lima_pbi_counts
    Array[String] source_dataset_names
    Array[String] hifi_barcodes
    Array[File] cdna_biosample_csvs
    File indexed_primers
    RuntimeAttributes runtime_attributes
  }

  Int threads = 1
  # One PBI is resident at a time, peaking near 3x its compressed size. Each
  # FLNC PBI holds a subset of one lima PBI's reads, so the lima total bounds it.
  Int mem_gb = 4 + ceil(4 * size(lima_bam_pbis, "GB"))
  Int total_mem_mb = (mem_gb * 1024) + runtime_attributes.add_memory_mb
  File lima_metadata = write_json(object {
    source_dataset_names: source_dataset_names,
    hifi_barcodes: hifi_barcodes,
    lima_pbi_counts: lima_pbi_counts
  })

  command <<<
    set -euo pipefail

    printf '%s\n' "~{sep="\" \"" refine_report_rows}" > refine_rows.txt
    printf '%s\n' "~{sep="\" \"" flnc_bam_pbis}" > flnc_pbis.txt
    printf '%s\n' "~{sep="\" \"" lima_bam_pbis}" > lima_pbis.txt
    printf '%s\n' "~{sep="\" \"" cdna_biosample_csvs}" > biosample_csvs.txt

    python3 - refine_rows.txt flnc_pbis.txt lima_pbis.txt biosample_csvs.txt \
      "~{lima_metadata}" "~{indexed_primers}" isoseq_refine.report.json <<'PY'
    import csv
    import gzip
    import io
    import json
    import os
    import struct
    import sys
    import uuid

    import matplotlib
    import numpy as np

    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import PowerNorm
    from matplotlib.ticker import FuncFormatter, MaxNLocator

    (
        row_manifest,
        flnc_manifest,
        lima_manifest,
        csv_manifest,
        metadata_path,
        primers_path,
        output_json,
    ) = sys.argv[1:]
    report_id = 'isoseq_refine'
    PBI_FLAGS_BARCODE = 4
    metric_names = (
        ('num_reads_fl', "Reads with 5' and 3' Primers"),
        ('num_reads_flnc', "Non-Concatamer Reads with 5' and 3' Primers"),
        ('num_reads_flnc_polya', "Non-Concatamer Reads with 5' and 3' Primers and Poly-A Tail"),
    )


    def fail(message):
        sys.exit(message)


    def read_json(path, label):
        try:
            with open(path, encoding='utf-8') as fh:
                return json.load(fh)
        except (OSError, ValueError) as exc:
            fail(f'{label}: cannot read JSON: {exc}')


    def read_paths(path):
        with open(path, encoding='utf-8') as fh:
            return [line.rstrip('\n') for line in fh if line.rstrip('\n')]


    def require_string(value, label):
        if not isinstance(value, str) or not value:
            fail(f'{label} must be a non-empty string')
        return value


    def require_uuid(value, label):
        require_string(value, label)
        try:
            uuid.UUID(value)
        except ValueError:
            fail(f'{label} is not a valid UUID: {value!r}')
        return value


    def empty_stats():
        return {'count': 0, 'bases': 0, 'quality': [0] * 101, 'histogram': {}, 'max_length': 0}


    def add_stats(target, source):
        target['count'] += source['count']
        target['bases'] += source['bases']
        target['max_length'] = max(target['max_length'], source['max_length'])
        for score, count in enumerate(source['quality']):
            target['quality'][score] += count
        for bin_index, count in source['histogram'].items():
            target['histogram'][bin_index] = target['histogram'].get(bin_index, 0) + count


    def reduce_reads(lengths, qualities=None, histogram=False):
        stats = empty_stats()
        stats['count'] = len(lengths)
        stats['bases'] = int(lengths.sum(dtype=np.int64))
        stats['max_length'] = int(lengths.max()) if len(lengths) else 0
        if qualities is not None:
            stats['quality'] = np.bincount(qualities.astype(np.int64), minlength=101).tolist()
        if histogram:
            counts = np.bincount(lengths // 100)
            stats['histogram'] = {int(b): int(counts[b]) for b in np.flatnonzero(counts)}
        return stats


    def read_column(fh, count, dtype):
        values = np.empty(count, dtype)
        view = memoryview(values).cast('B')
        filled = 0
        while filled < len(view):
            read = fh.readinto(view[filled:])
            if not read:
                raise ValueError('truncated PBI column')
            filled += read
        return values


    def load_index(path, label):
        # pbcore's PacBioBamIndex grows each column with bytes +=, which is
        # quadratic in read count. BGZF is multi-member gzip, so stream only
        # the unaligned basic and barcode columns this report reads.
        try:
            with gzip.open(path, 'rb') as fh:
                magic, patch, minor, major, flags, count = struct.unpack('<4sBBBxHI18x', fh.read(32))
                if magic != b'PBI\x01':
                    fail(f'{label}: corrupt PBI magic')
                if (major, minor, patch) < (3, 0, 1) or flags & ~PBI_FLAGS_BARCODE:
                    fail(f'{label}: PBI {os.path.basename(path)!r} must be an unaligned version 3.0.1+ index')
                fh.seek(4 * count, io.SEEK_CUR)
                q_start = read_column(fh, count, '<i4')
                lengths = read_column(fh, count, '<i4').astype(np.int64) - q_start
                if np.any(q_start < 0) or np.any(lengths < 0):
                    fail(f'{label}: incompatible PBI read lengths')
                del q_start
                barcodes = None
                if flags & PBI_FLAGS_BARCODE:
                    fh.seek(17 * count, io.SEEK_CUR)
                    barcodes = tuple(read_column(fh, count, dtype) for dtype in ('<i2', '<i2', '<i1'))
            return lengths, barcodes
        except Exception as exc:
            fail(f'{label}: cannot read PBI {os.path.basename(path)!r}: {exc}')


    row_paths, flnc_paths = read_paths(row_manifest), read_paths(flnc_manifest)
    lima_paths, csv_paths = read_paths(lima_manifest), read_paths(csv_manifest)
    metadata = read_json(metadata_path, 'Lima metadata')
    sources = metadata.get('source_dataset_names', [])
    hifis = metadata.get('hifi_barcodes', [])
    pbi_counts = metadata.get('lima_pbi_counts', [])
    outer_label = f'outer identities {list(zip(sources, hifis))!r}'
    if not sources or not (len(sources) == len(hifis) == len(pbi_counts) == len(csv_paths)):
        fail(
            f'{outer_label}: source_dataset_names, hifi_barcodes, lima_pbi_counts and CSV arrays must be nonempty and equally sized'
        )
    if any(type(n) is not int or n < 0 for n in pbi_counts) or sum(pbi_counts) != len(lima_paths):
        fail(f'{outer_label}: lima_pbi_counts must be nonnegative and sum to lima_bam_pbis length')
    if len(row_paths) != len(flnc_paths):
        fail(f'{outer_label}: refine_report_rows and flnc_bam_pbis lengths must match')

    with open(primers_path, encoding='utf-8') as fh:
        headers = [line[1:].split() for line in fh if line.startswith('>')]
    if not headers or any(not fields for fields in headers):
        fail('indexed_primers must contain unique nonempty FASTA IDs')
    primers = [fields[0] for fields in headers]
    if len(set(primers)) != len(primers):
        fail('indexed_primers must contain unique nonempty FASTA IDs')
    primer_indices = {name: i for i, name in enumerate(primers)}


    def normalize_pair(name, label):
        parts = name.split('--')
        if len(parts) != 2 or any(part not in primer_indices for part in parts):
            fail(f'{label}: unknown cDNA primer pair {name!r}')
        return tuple(sorted(primer_indices[part] for part in parts))


    # Orientation equivalence follows pbreports, but each PBI uses masks rather
    # than contiguous ranges. Adaptation sources: internal/docs/upstream-sources.md.
    outers = {}
    hifi_order = []
    selected_order = {}
    canonical = {}
    pbi_offset = 0
    assigned_total, unbarcoded_total, flnc_total = empty_stats(), empty_stats(), empty_stats()
    for source, hifi, count, csv_path in zip(sources, hifis, pbi_counts, csv_paths):
        require_string(source, 'Source Dataset')
        require_string(hifi, f'{source}: HiFi Barcode')
        key = (source, hifi)
        label = '/'.join(key)
        if key in outers:
            fail(f'{label}: duplicate outer dataset identity')
        if hifi not in hifi_order:
            hifi_order.append(hifi)
            selected_order[hifi] = []
        selected = {}
        with open(csv_path, newline='', encoding='utf-8-sig') as fh:
            reader = csv.reader(fh)
            header = tuple(cell.strip() for cell in next(reader, []))
            if header != ('Barcodes', 'Bio Sample'):
                fail(f'{label}: biosample CSV has header {header!r}; expected Barcodes,Bio Sample')
            for line_no, raw in enumerate(reader, 2):
                fields = [cell.strip() for cell in raw]
                if not fields or all(not cell for cell in fields):
                    continue
                if len(fields) != 2 or not all(fields):
                    fail(f'{label}: biosample CSV row {line_no}: expected two nonempty fields')
                pair_name, sample = fields
                pair = normalize_pair(pair_name, label)
                if pair in selected:
                    fail(f'{label}/{pair_name}: duplicate barcode pair in biosample CSV')
                selected[pair] = (pair_name, sample)
                group_key = (hifi, pair)
                if group_key in canonical and canonical[group_key][1] != sample:
                    fail(f'{label}/{pair_name}: contradictory Bio Sample mapping')
                if group_key not in canonical:
                    canonical[group_key] = (pair_name, sample)
                    selected_order[hifi].append(pair)
        if not selected:
            fail(f'{label}: biosample CSV has no data rows')
        paths = lima_paths[pbi_offset : pbi_offset + count]
        pbi_offset += count
        unbarcoded = [path for path in paths if os.path.basename(path).endswith('.unbarcoded.bam.pbi')]
        if len(unbarcoded) != 1:
            fail(f'{label}/No Primer: expected one unbarcoded PBI; observed {len(unbarcoded)}')
        outer = {'selected': selected, 'assigned': {}, 'refine': {}, 'unbarcoded': empty_stats()}
        outers[key] = outer
        for path in paths:
            index_label = f'{label}/{os.path.basename(path)}'
            lengths, barcodes = load_index(path, index_label)
            if path == unbarcoded[0]:
                outer['unbarcoded'] = reduce_reads(lengths)
                add_stats(unbarcoded_total, outer['unbarcoded'])
            elif len(lengths):
                if barcodes is None:
                    fail(f'{index_label}: assigned PBI requires barcode/quality columns')
                forward, reverse, quality = barcodes
                if (
                    np.any(forward < 0)
                    or np.any(reverse < 0)
                    or np.any(forward >= len(primers))
                    or np.any(reverse >= len(primers))
                ):
                    fail(f'{index_label}: assigned primer index outside indexed_primers')
                if np.any(quality < 0) or np.any(quality > 100):
                    bad = int(np.flatnonzero((quality < 0) | (quality > 100))[0])
                    pair_name = '--'.join(primers[i] for i in sorted((int(forward[bad]), int(reverse[bad]))))
                    fail(f'{label}/{pair_name}: assigned primer quality must be in [0,100]')
                codes = np.minimum(forward, reverse).astype(np.int64) * len(primers) + np.maximum(forward, reverse)
                for code in np.unique(codes):
                    first, second = divmod(int(code), len(primers))
                    pair = (first, second)
                    mask = ((forward == first) & (reverse == second)) | ((forward == second) & (reverse == first))
                    stats = reduce_reads(lengths[mask], quality[mask])
                    add_stats(outer['assigned'].setdefault(pair, empty_stats()), stats)
                    add_stats(assigned_total, stats)
                del forward, reverse, quality, codes, mask
            del lengths, barcodes

    row_fields = {'report_json', 'bio_sample', 'source_dataset', 'hifi_barcode', 'cdna_barcode'}
    seen_keys, dataset_uuids, report_tags = set(), set(), set()
    uuid_reports = []
    for row_number, (row_path, flnc_path) in enumerate(zip(row_paths, flnc_paths), 1):
        row = read_json(row_path, f'refine report row {row_number}')
        if not isinstance(row, dict) or set(row) != row_fields:
            fail(f'refine report row {row_number}: expected fields {sorted(row_fields)!r}')
        raw_identity = '/'.join(str(row[field]) for field in ('source_dataset', 'hifi_barcode', 'cdna_barcode'))
        for field in row_fields:
            require_string(row[field], f'{raw_identity}: refine report row {row_number} {field}')
        key = (row['source_dataset'], row['hifi_barcode'], row['cdna_barcode'])
        label = '/'.join(key)
        if key in seen_keys:
            fail(f'duplicate refine report row key (Source Dataset, HiFi Barcode, cDNA Barcode)={key!r}')
        seen_keys.add(key)
        outer = outers.get(key[:2])
        pair = normalize_pair(key[2], label)
        if outer is None or outer['selected'].get(pair) != (key[2], row['bio_sample']):
            fail(f'{label}: refine identity does not match outer dataset biosample CSV')
        if pair in outer['refine']:
            fail(f'duplicate refine report row key (Source Dataset, HiFi Barcode, cDNA Barcode)={key!r}')
        try:
            report = json.loads(row['report_json'])
        except ValueError as exc:
            fail(f'{label}: report_json is invalid JSON: {exc}')
        if not isinstance(report, dict) or report.get('id') != report_id:
            fail(f'{label}: report id must be {report_id!r}')
        report_uuid = require_uuid(report.get('uuid'), f'{label}: report uuid')
        attributes = report.get('attributes')
        if not isinstance(attributes, list):
            fail(f'{label}: report attributes must be an array')
        metrics = {}
        for attribute in attributes:
            if not isinstance(attribute, dict):
                fail(f'{label}: report attribute must be an object')
            leaf = require_string(attribute.get('id'), f'{label}: attribute id').rsplit('.', 1)[-1]
            if leaf in metrics:
                fail(f'{label}: duplicate metric id leaf {leaf!r}')
            metrics[leaf] = attribute.get('value')
        if metrics.get('sample_name') != row['bio_sample']:
            fail(
                f'{label}: refine Sample Name {metrics.get("sample_name")!r} does not match Bio Sample {row["bio_sample"]!r}'
            )
        counts = {}
        for leaf, _ in metric_names:
            if leaf not in metrics:
                fail(f'{label}: missing required refine metric {leaf}')
            value = metrics[leaf]
            if type(value) is not int or value < 0:
                fail(f'{label}: required refine metric {leaf} must be a nonnegative integer')
            counts[leaf] = value
        uuids, tags = report.get('dataset_uuids', []), report.get('tags', [])
        if not isinstance(uuids, list) or not isinstance(tags, list) or any(not isinstance(tag, str) for tag in tags):
            fail(f'{label}: dataset_uuids and tags must be arrays; tags must be strings')
        dataset_uuids.update(require_uuid(value, f'{label}: dataset UUID') for value in uuids)
        report_tags.update(tags)
        lengths, _ = load_index(flnc_path, label)
        flnc = reduce_reads(lengths, histogram=True)
        del lengths
        add_stats(flnc_total, flnc)
        outer['refine'][pair] = {'counts': counts, 'flnc': flnc}
        uuid_reports.append({'identity': key, 'uuid': report_uuid, 'counts': counts, 'flnc': flnc})

    # Table identity stays technical; plot identity pools only exact
    # HiFi/cDNA/Bio Sample matches, never the Bio Sample string alone.
    table_rows, pooled = [], {}
    selected_groups = set()
    for (source, hifi), outer in outers.items():
        unexpected = sorted(
            set(outer['assigned']) - set(outer['selected']),
            key=lambda pair: '--'.join(primers[i] for i in pair),
        )
        for pair in list(outer['selected']) + unexpected:
            assigned = outer['assigned'].get(pair, empty_stats())
            refine = outer['refine'].get(pair)
            pair_name, sample = outer['selected'].get(pair, ('--'.join(primers[i] for i in pair), 'UNASSIGNED'))
            label = f'{source}/{hifi}/{pair_name}'
            if pair in outer['selected'] and assigned['count'] and refine is None:
                fail(f'{label}: observed selected pair is missing its refine report row')
            if not assigned['count'] and refine is None:
                continue
            if pair in outer['selected']:
                canonical_name, canonical_sample = canonical[(hifi, pair)]
            else:
                canonical_name, canonical_sample = pair_name, sample
            group_key = (hifi, canonical_name, canonical_sample)
            if pair in outer['selected']:
                selected_groups.add(group_key)
            group = pooled.setdefault(group_key, {'assigned': empty_stats(), 'flnc': empty_stats(), 'pair': pair})
            add_stats(group['assigned'], assigned)
            if refine is not None:
                add_stats(group['flnc'], refine['flnc'])
            quality_sum = sum(score * n for score, n in enumerate(assigned['quality']))
            table_rows.append(
                {
                    'biosample': sample,
                    'source_dataset': source,
                    'hifi_barcode': hifi,
                    'barcode': pair_name,
                    'number_of_reads': assigned['count'],
                    'mean_bcqual': quality_sum / assigned['count'] if assigned['count'] else 0.0,
                    **(refine['counts'] if refine is not None else {leaf: None for leaf, _ in metric_names}),
                }
            )
        table_rows.append(
            {
                'biosample': 'UNASSIGNED',
                'source_dataset': source,
                'hifi_barcode': hifi,
                'barcode': 'No Primer',
                'number_of_reads': outer['unbarcoded']['count'],
                'mean_bcqual': 0.0,
                **{leaf: 0 for leaf, _ in metric_names},
            }
        )

    groups = []
    for hifi in hifi_order:
        selected_ranks = {pair: i for i, pair in enumerate(selected_order[hifi])}
        keys = [key for key in pooled if key[0] == hifi]
        keys.sort(
            key=lambda key: (
                selected_ranks[pooled[key]['pair']] if key in selected_groups else len(selected_ranks),
                key[1],
                key[2],
            )
        )
        groups.extend((key, pooled[key]) for key in keys)
    observed = [(key, group) for key, group in groups if group['assigned']['count']]
    group_counts = [group['assigned']['count'] for _, group in observed]
    total_count = assigned_total['count'] + unbarcoded_total['count']
    total_bases = assigned_total['bases'] + unbarcoded_total['bases']
    summary = [
        ('n_reads_total', 'Reads', total_count),
        *[
            (
                leaf,
                name + (' (FLNC Reads)' if leaf == 'num_reads_flnc_polya' else ''),
                sum(row[leaf] or 0 for row in table_rows),
            )
            for leaf, name in metric_names
        ],
        (
            'mean_flnc_read_length',
            'Mean Length of FLNC Reads',
            flnc_total['bases'] // flnc_total['count'] if flnc_total['count'] else 0,
        ),
        ('n_barcodes', 'Unique HiFi/cDNA Barcode Combinations', len(observed)),
        (
            'mean_reads',
            'Mean Reads per HiFi/cDNA Barcode Combination',
            assigned_total['count'] // len(observed) if observed else 0,
        ),
        ('max_reads', 'Max. Reads per HiFi/cDNA Barcode Combination', max(group_counts, default=0)),
        ('min_reads', 'Min. Reads per HiFi/cDNA Barcode Combination', min(group_counts, default=0)),
        ('n_unbarcoded_reads', 'Reads without Primers', unbarcoded_total['count']),
        (
            'percent_barcoded_bases',
            'Percent Bases in Reads with Primers',
            assigned_total['bases'] / total_bases if total_bases else 0,
        ),
        (
            'percent_barcoded_reads',
            'Percent Reads with Primers',
            assigned_total['count'] / total_count if total_count else 0,
        ),
    ]
    table_columns = (
        ('biosample', 'Bio Sample Name'),
        ('source_dataset', 'Source Dataset'),
        ('hifi_barcode', 'HiFi Barcode'),
        ('barcode', 'Primer Name (cDNA Barcode)'),
        ('number_of_reads', 'CCS Reads'),
        ('mean_bcqual', 'Mean Primer Quality'),
        *metric_names,
    )
    uuid_seed = json.dumps(
        {
            'outer_identities': list(outers),
            'source_reports': uuid_reports,
            'table_rows': table_rows,
            'assigned': assigned_total,
            'unbarcoded': unbarcoded_total,
            'flnc': flnc_total,
            'groups': [(key, group) for key, group in groups],
        },
        sort_keys=True,
        separators=(',', ':'),
    )
    # Plot algorithms/styles adapted locally; no runtime pbreports dependency.
    # Source mappings are recorded in internal/docs/upstream-sources.md.
    plt.rcParams.update({'font.size': 10, 'axes.titlesize': 12, 'axes.labelsize': 11})
    count_formatter = FuncFormatter(lambda value, _position: f'{int(value):,}')
    plot_groups = []


    def new_axes(title, xlabel, ylabel):
        fig, ax = plt.subplots(figsize=(10, 6), constrained_layout=True)
        fig.set_facecolor('white')
        ax.set_facecolor('white')
        ax.set_axisbelow(True)
        ax.grid(color='#e0e0e0', linewidth=0.5)
        ax.set(title=title, xlabel=xlabel, ylabel=ylabel)
        ax.xaxis.set_major_formatter(count_formatter)
        ax.yaxis.set_major_formatter(count_formatter)
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        return fig, ax


    def histogram(ax, values, bins, weights=None, color='#009ca2'):
        return ax.hist(values, bins=bins, weights=weights, color=color, edgecolor='none', rwidth=0.8)


    def save_plot(fig, group_leaf, plot_leaf, title, caption, filename=None):
        image = filename or f'isoseq_refine_{plot_leaf}.png'
        thumbnail = image[:-4] + '_thumb.png'
        fig.savefig(image, dpi=180)
        fig.savefig(thumbnail, dpi=20)
        plt.close(fig)
        return {
            'id': f'{report_id}.{group_leaf}.{plot_leaf}',
            'title': title,
            'caption': caption,
            'plotType': 'image',
            'image': image,
            'thumbnail': thumbnail,
        }


    def add_plot_group(leaf, title, plots):
        if plots:
            plot_groups.append(
                {
                    'id': f'{report_id}.{leaf}',
                    'title': title,
                    'plots': plots,
                    'thumbnail': plots[0]['thumbnail'],
                }
            )


    if assigned_total['count']:
        read_stats = []
        title = 'Number of Reads per HiFi/cDNA Barcode Combination'
        fig, ax = new_axes(title, 'HiFi/cDNA Barcode Combination Rank', 'Count of Reads')
        ax.plot(np.arange(len(group_counts)) + 1, sorted(group_counts, reverse=True), color='#008eb9')
        ax.axhline(assigned_total['count'] / len(observed), color='red', label='Mean Number of Reads')
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.legend(loc='best')
        read_stats.append(save_plot(fig, 'read_stats', 'nreads', title, title))

        title = 'HiFi/cDNA Barcode Combination Frequency Distribution'
        fig, ax = new_axes(title, 'Number of Reads with Primers', 'Number of HiFi/cDNA Barcode Combinations')
        histogram(ax, group_counts, max(1, min(len(observed), 40)))
        read_stats.append(save_plot(fig, 'read_stats', 'nreads_histogram', title, title))

        mean_lengths = [group['assigned']['bases'] // group['assigned']['count'] for _, group in observed]
        max_mean_length = max(mean_lengths)
        bin_width = next((width for width in (100, 200, 500, 1000, 2000) if max_mean_length / width < 100), 2000)
        mean_edges = np.arange(0, (max_mean_length // bin_width + 1) * bin_width + 1, bin_width)
        title = 'Mean Read Length Distribution'
        fig, ax = new_axes(title, 'Mean Read Length (bp)', 'Number of HiFi/cDNA Barcode Combinations')
        histogram(ax, mean_lengths, mean_edges)
        read_stats.append(
            save_plot(
                fig,
                'read_stats',
                'readlength_histogram',
                title,
                'Read-weighted mean assigned read length per HiFi/cDNA Barcode Combination',
            )
        )
        add_plot_group('read_stats', 'Primer Read Statistics', read_stats)

        quality_counts = np.asarray(assigned_total['quality'], dtype=np.int64)
        lower = min(int(np.flatnonzero(quality_counts)[0]), 80)
        lower -= lower % 2
        bin_width = 1 if lower > 50 else 2
        title = 'Primer Quality Score Distribution'
        fig, ax = new_axes(title, 'Primer Quality Score', 'Number of Reads with Primers')
        histogram(ax, np.arange(101), np.arange(lower, 102, bin_width), weights=quality_counts)
        ax.axvline(lower, color='red')
        quality_plot = save_plot(
            fig, 'bq_plots', 'bq_histogram', title, 'Assigned primer quality scores from all cDNA Lima reads with primers'
        )
        add_plot_group('bq_plots', 'Primer Quality Scores', [quality_plot])

    if flnc_total['count']:
        flnc_group_leaf = 'fulllength_nonconcatemer_readlength_group'
        flnc_title = 'Length of Full-length Non-Concatemer Reads'
        flnc_plots = []
        final_edge = (flnc_total['max_length'] // 100 + 1) * 100
        full_edges = np.arange(0, final_edge + 1, 100)
        # Sigma is one 100-bp bin; truncate the Gaussian at four standard deviations.
        kernel = np.exp(-0.5 * np.arange(-4, 5, dtype=float) ** 2)
        kernel /= kernel.sum()
        display_bin_count = 1
        for _key, group in groups:
            stats = group['flnc']
            if stats['count']:
                # ceil(99.5% * count), without floating-point rank rounding.
                required = (199 * stats['count'] + 199) // 200
                covered = 0
                for bin_index, count in sorted(stats['histogram'].items()):
                    covered += count
                    if covered >= required:
                        display_bin_count = max(display_bin_count, bin_index + 1)
                        break
        heatmap_edges = np.arange(0, display_bin_count * 100 + 1, 100)
        cmap = plt.colormaps['YlOrBr'].copy()
        # Transparent zero-intensity cells show white; empty groups have a gray underlay.
        cmap.set_bad('white', alpha=0)
        heatmap_columns = [(key, group) for key, group in groups if group['assigned']['count'] or group['flnc']['count']]
        page_size = 96
        page_count = (len(heatmap_columns) + page_size - 1) // page_size
        for page_index, start in enumerate(range(0, len(heatmap_columns), page_size)):
            columns = heatmap_columns[start : start + page_size]
            density = np.zeros((display_bin_count, len(columns)), dtype=float)
            mode_positions, mode_lengths = [], []
            for column, (_key, group) in enumerate(columns):
                stats = group['flnc']
                if stats['count']:
                    counts = np.zeros(len(full_edges) - 1, dtype=float)
                    for bin_index, count in stats['histogram'].items():
                        counts[bin_index] = count
                    smoothed = np.convolve(np.pad(counts, 4), kernel, mode='valid')
                    density[:, column] = smoothed[:display_bin_count] / smoothed.max()
                    mode_positions.append(column + 0.5)
                    mode_lengths.append((int(np.argmax(smoothed)) + 0.5) * 100)
            fig, ax = plt.subplots(figsize=(max(12, 0.38 * len(columns)), 9), constrained_layout=True)
            fig.set_facecolor('white')
            ax.set_facecolor('white')
            for column, (_key, group) in enumerate(columns):
                if not group['flnc']['count']:
                    ax.axvspan(column, column + 1, facecolor='#d3d3d3', edgecolor='none', zorder=0)
            title = 'FLNC read-length density'
            if page_count > 1:
                title += f' (page {page_index + 1} of {page_count})'
            mesh = ax.pcolormesh(
                np.arange(len(columns) + 1),
                heatmap_edges,
                np.ma.masked_equal(density, 0, copy=False),
                shading='flat',
                cmap=cmap,
                norm=PowerNorm(gamma=1.5, vmin=0, vmax=1),
            )
            if mode_positions:
                ax.scatter(
                    mode_positions,
                    mode_lengths,
                    color='red',
                    edgecolor='white',
                    linewidth=0.5,
                    s=15,
                    zorder=3,
                    label='Smoothed modal bin',
                )
                fig.legend(loc='outside upper right', fontsize=8)
            ax.set(
                title=title,
                xlabel='HiFi Barcode / cDNA Barcode / Bio Sample',
                ylabel='FLNC read length (bp)',
                ylim=(0, display_bin_count * 100),
            )
            ax.set_xticks(np.arange(len(columns)) + 0.5)
            ax.set_xticklabels(['\n'.join(key) for key, _ in columns], rotation=90, fontsize=8, linespacing=1.0)
            for column in range(1, len(columns)):
                if columns[column - 1][0][0] != columns[column][0][0]:
                    ax.axvline(column, color='#777777', linewidth=0.8)
            ax.yaxis.set_major_formatter(count_formatter)
            colorbar = fig.colorbar(mesh, ax=ax, ticks=[0, 0.2, 0.4, 0.6, 0.8, 1])
            colorbar.set_label('Density relative to peak')
            flnc_plots.append(
                save_plot(
                    fig,
                    flnc_group_leaf,
                    f'flnc_density_heatmap_page_{page_index}',
                    title,
                    f'Pooled groups {start + 1}–{start + len(columns)} of {len(heatmap_columns)}; '
                    'columns identify HiFi barcode, cDNA barcode and Bio Sample; lines separate HiFi barcodes. Gray means no FLNC reads.',
                    f'isoseq_refine_flnc_density_page_{page_index}.png',
                )
            )

        for group_index, (key, group) in enumerate(groups):
            stats = group['flnc']
            if not stats['count']:
                continue
            identity = ' | '.join(key)
            title = f'{flnc_title}\n{identity}'
            fig, ax = new_axes(title, 'Read Length (bp)', 'Number of Reads')
            positions = np.asarray(list(stats['histogram']), dtype=np.int64) * 100
            frequencies = np.asarray(list(stats['histogram'].values()), dtype=np.int64)
            histogram(ax, positions, full_edges, weights=frequencies, color='#008eb9')
            ax.set_xlim(0, final_edge)
            flnc_plots.append(
                save_plot(
                    fig,
                    flnc_group_leaf,
                    f'fulllength_nonconcatemer_readlength_hist_group_{group_index}',
                    title,
                    f'FLNC read lengths pooled across source datasets for {identity}',
                    f'isoseq_refine_flnc_lengths_group_{group_index}.png',
                )
            )
        add_plot_group(flnc_group_leaf, flnc_title, flnc_plots)


    def format_summary_value(leaf, value):
        if leaf.startswith('percent_'):
            return f'{100 * value:.2f}'.rstrip('0').rstrip('.') + '%'
        return f'{value:,d}'


    def format_table_value(leaf, value):
        if value is None or isinstance(value, str):
            return value
        if leaf == 'mean_bcqual':
            return f'{value:.1f}'
        return f'{value:,d}'


    combined_report = {
        'id': report_id,
        'version': '1.0.1',
        'title': 'Iso-Seq Refine',
        'uuid': str(uuid.uuid5(uuid.NAMESPACE_URL, f'kinnex-isoseq-wdl/isoseq_refine/{uuid_seed}')),
        'attributes': [
            {'id': f'{report_id}.{leaf}', 'name': name, 'value': format_summary_value(leaf, value)}
            for leaf, name, value in summary
        ],
        'tables': [
            {
                'id': f'{report_id}.primer_table',
                'title': 'Primer Data',
                'columns': [
                    {
                        'id': f'{report_id}.primer_table.{leaf}',
                        'header': header,
                        'values': [format_table_value(leaf, row[leaf]) for row in table_rows],
                    }
                    for leaf, header in table_columns
                ],
            }
        ],
        'plotGroups': plot_groups,
        'dataset_uuids': sorted(dataset_uuids),
        'tags': sorted(report_tags),
    }
    with open(output_json, 'w', encoding='utf-8') as fh:
        json.dump(combined_report, fh, indent=2, sort_keys=True, allow_nan=False)
        fh.write('\n')
    PY
  >>>

  output {
    File refine_summary_report = "isoseq_refine.report.json"
    Array[File] refine_summary_report_resources = glob("*.png")
  }

  runtime {
    cpu: threads
    memory: total_mem_mb + " MB"
    docker: runtime_attributes.container_registry + "/pbcore@sha256:47fc6f1174605be9a8a932f17c15d1a29fc5e6b89ea6809f32995318c7afe1e7"
    maxRetries: runtime_attributes.max_retries
  }
}
