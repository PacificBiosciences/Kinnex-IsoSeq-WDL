version 1.0

import "../../rna_structs.wdl"

task generate_read_segmentation_report {
  meta {
    description: "Aggregate Skera segmentation metrics and generate the Read Segmentation report."
    outputs: {
      segmentation_summary_csv: {
        description: "Per-HiFi-barcode Skera segmentation metrics"
      },
      segmentation_summary_report: {
        description: "Read Segmentation metrics and per-barcode plots"
      },
      segmentation_summary_report_resources: {
        description: "PNG images and thumbnails referenced by the Read Segmentation report"
      }
    }
  }

  parameter_meta {
    hifi_barcodes: {
      description: "HiFi barcode for each normalized input"
    }
    skera_summary_csvs: {
      description: "Skera summary CSV for each normalized input"
    }
    skera_summary_jsons: {
      description: "Skera summary report JSON for each normalized input"
    }
    read_lengths_csvs: {
      description: "Skera read-length CSV for each normalized input"
    }
    ligations_csvs: {
      description: "Skera ligation CSV for each normalized input"
    }
    segmentation_adapters: {
      description: "Segmentation adapter FASTA"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    Array[String] hifi_barcodes
    Array[File] skera_summary_csvs
    Array[File] skera_summary_jsons
    Array[File] read_lengths_csvs
    Array[File] ligations_csvs
    File segmentation_adapters
    RuntimeAttributes runtime_attributes
  }

  Int threads = 1
  # Each read-length CSV loads at about 2x its size, then only compact
  # per-read arrays (about 0.25x every CSV) stay resident until plotting.
  Int mem_gb = 4 + ceil(2 * size(read_lengths_csvs, "GB"))
  Int total_mem_mb = (mem_gb * 1024) + runtime_attributes.add_memory_mb

  command <<<
    set -euo pipefail

    printf '%s\n' "~{sep="\" \"" hifi_barcodes}" > hifi_barcodes.txt
    printf '%s\n' "~{sep="\" \"" skera_summary_csvs}" > skera_summary_csvs.txt
    printf '%s\n' "~{sep="\" \"" skera_summary_jsons}" > skera_summary_jsons.txt
    printf '%s\n' "~{sep="\" \"" read_lengths_csvs}" > read_lengths_csvs.txt
    printf '%s\n' "~{sep="\" \"" ligations_csvs}" > ligations_csvs.txt

    python3 - \
      hifi_barcodes.txt \
      skera_summary_csvs.txt \
      skera_summary_jsons.txt \
      read_lengths_csvs.txt \
      ligations_csvs.txt \
      "~{segmentation_adapters}" \
      ~{length(hifi_barcodes)} \
      ~{length(skera_summary_csvs)} \
      ~{length(skera_summary_jsons)} \
      ~{length(read_lengths_csvs)} \
      ~{length(ligations_csvs)} <<'PY'
    # Upstream sources: internal/docs/upstream-sources.md#read-segmentation-report
    # Copyright (c) 2016, Pacific Biosciences. All rights reserved.
    # The Clear BSD terms are provided in this repository's LICENSE.
    import csv
    import json
    import math
    import os.path as op
    import sys
    import uuid
    from collections import OrderedDict

    import matplotlib
    import numpy as np
    import pandas as pd

    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib import ticker
    from matplotlib.patches import Rectangle as Rect

    (
        hifi_barcodes_manifest,
        skera_summary_csvs_manifest,
        skera_summary_jsons_manifest,
        read_lengths_csvs_manifest,
        ligations_csvs_manifest,
        segmentation_adapters,
        *expected_counts_text,
    ) = sys.argv[1:]

    output_csv = 'segmentation.summary.csv'
    output_json = 'segmentation.summary.report.json'
    report_id = 'read_segmentation'
    report_version = '1.0.1'
    report_title = 'Read Segmentation'
    metric_rows = (
        ('Input Reads', 'reads'),
        ('Segmented Reads (S-Reads)', 's_reads'),
        ('Mean Length of S-Reads', 'mean_len_s_reads'),
        ('Percentage of Reads with Full Array', 'percent_full_array'),
        ('Mean Array Size (Concatenation Factor)', 'mean_array_size'),
    )
    csv_header = ('HiFi barcode', *(name for name, _metric_id in metric_rows))
    read_length_header = (
        'zmw',
        'original_read_length',
        'm_read_length',
        'concatenation_factor',
    )
    ligation_header = ('adapter_1', 'adapter_2', 'ligations')
    plot_specs = (
        (
            'segmentation_statistics',
            'Segmentation Statistics',
            (
                (
                    'concatemer_histogram',
                    'Histogram of concatemers',
                    'Histogram of concatemers',
                ),
                (
                    'ligations_heatmap',
                    'Heatmap of concatemer ligations',
                    'Heatmap of concatemer ligations',
                ),
            ),
        ),
        (
            'read_length_group',
            'Length of Reads',
            (('hifi_read_length', 'Length of Reads', 'Length of Reads'),),
        ),
        (
            'sread_length_group',
            'S-read Length Distribution',
            (('sread_length', 'Length of S-read', 'Length of S-reads'),),
        ),
    )


    def fail(message):
        sys.exit(message)


    def read_manifest(path, expected_count, label):
        with open(path, 'r', encoding='utf-8') as manifest_fh:
            text = manifest_fh.read()
        lines = text.splitlines()
        if expected_count == 0 and lines == ['']:
            lines = []
        if len(lines) != expected_count:
            fail(f'{label} manifest contains {len(lines)} entries; WDL array contains {expected_count}')
        return lines


    def require_uuid(value, label):
        if not isinstance(value, str) or not value:
            fail(f'{label} must be a non-empty UUID string')
        try:
            return str(uuid.UUID(value))
        except ValueError:
            fail(f'{label} is not a valid UUID: {value!r}')


    def parse_number(value, label, *, integer=False, nonnegative=True):
        if not isinstance(value, str) or not value:
            fail(f'{label} must be a non-empty numeric value')
        try:
            number = float(value)
        except ValueError:
            fail(f'{label} is not numeric: {value!r}')
        if not math.isfinite(number):
            fail(f'{label} must be finite: {value!r}')
        if nonnegative and number < 0:
            fail(f'{label} must be non-negative: {value!r}')
        if integer:
            if not number.is_integer():
                fail(f'{label} must be an integer: {value!r}')
            return int(number)
        return number


    def parse_json_number(value, label):
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            fail(f'{label} must be numeric')
        number = float(value)
        if not math.isfinite(number) or number < 0:
            fail(f'{label} must be a finite non-negative number')
        return value


    def stable_calculation(numerator, denominator):
        if denominator == 0:
            return '0'
        return format(numerator / denominator, '.15g')


    def numeric_csv_value(value):
        number = float(value)
        if number.is_integer() and all(marker not in value.lower() for marker in ('.', 'e')):
            return int(number)
        return number


    def load_summary_csv(path, input_label):
        try:
            with open(path, 'r', encoding='utf-8', newline='') as csv_fh:
                rows = list(csv.reader(csv_fh))
        except OSError as exc:
            fail(f'{input_label}: could not read Skera summary CSV: {exc}')
        expected_names = [name for name, _metric_id in metric_rows]
        if len(rows) != len(metric_rows) or any(len(row) != 2 for row in rows):
            fail(f'{input_label}: Skera summary CSV must contain exactly five two-column rows')
        observed_names = [row[0] for row in rows]
        if observed_names != expected_names:
            fail(f'{input_label}: Skera summary CSV metric rows are {observed_names!r}; expected {expected_names!r}')
        raw_values = [row[1] for row in rows]
        return {
            'raw_values': raw_values,
            'reads': parse_number(raw_values[0], f'{input_label}: Input Reads', integer=True),
            's_reads': parse_number(
                raw_values[1],
                f'{input_label}: Segmented Reads (S-Reads)',
                integer=True,
            ),
            'mean_len_s_reads': parse_number(
                raw_values[2],
                f'{input_label}: Mean Length of S-Reads',
            ),
            'percent_full_array': parse_number(
                raw_values[3],
                f'{input_label}: Percentage of Reads with Full Array',
            ),
            'mean_array_size': parse_number(
                raw_values[4],
                f'{input_label}: Mean Array Size (Concatenation Factor)',
            ),
        }


    def load_skera_json(path, input_label):
        try:
            with open(path, 'r', encoding='utf-8') as json_fh:
                report = json.load(json_fh)
        except (OSError, json.JSONDecodeError) as exc:
            fail(f'{input_label}: malformed Skera summary JSON: {exc}')
        if not isinstance(report, dict):
            fail(f'{input_label}: Skera summary JSON must contain an object')
        if report.get('id') != 'skera_summary':
            fail(f"{input_label}: Skera summary JSON id must be 'skera_summary'; observed {report.get('id')!r}")
        if report.get('version') != report_version:
            fail(
                f'{input_label}: Skera summary JSON version must be {report_version!r}; observed {report.get("version")!r}'
            )
        report_uuid = require_uuid(report.get('uuid'), f'{input_label}: report uuid')
        for field in ('tables', 'plotGroups'):
            value = report.get(field)
            if not isinstance(value, list) or value:
                fail(f'{input_label}: Skera summary JSON {field} must be an empty array')
        attributes = report.get('attributes')
        if not isinstance(attributes, list):
            fail(f'{input_label}: Skera summary JSON attributes must be an array')
        expected_ids = [metric_id for _name, metric_id in metric_rows]
        observed = {}
        for attribute_index, attribute in enumerate(attributes, start=1):
            attribute_label = f'{input_label}: attribute {attribute_index}'
            if not isinstance(attribute, dict):
                fail(f'{attribute_label} must be an object')
            raw_id = attribute.get('id')
            if not isinstance(raw_id, str) or not raw_id:
                fail(f'{attribute_label} id must be a non-empty string')
            metric_id = raw_id.rsplit('.', 1)[-1]
            if metric_id in observed:
                fail(f'{input_label}: duplicate Skera attribute leaf {metric_id!r}')
            if 'value' not in attribute:
                fail(f'{attribute_label} is missing value')
            parse_json_number(attribute['value'], f'{attribute_label} value')
            observed[metric_id] = attribute
        if list(observed) != expected_ids:
            fail(f'{input_label}: Skera attribute leaves are {list(observed)!r}; expected {expected_ids!r}')
        json_s_reads_value = float(observed['s_reads']['value'])
        if not json_s_reads_value.is_integer():
            fail(f'{input_label}: Skera summary JSON s_reads must be an integer')
        json_s_reads = int(json_s_reads_value)
        dataset_uuids = report.get('dataset_uuids', [])
        if not isinstance(dataset_uuids, list):
            fail(f'{input_label}: dataset_uuids must be an array')
        validated_dataset_uuids = [
            require_uuid(value, f'{input_label}: dataset_uuids[{index}]')
            for index, value in enumerate(dataset_uuids, start=1)
        ]
        tags = report.get('tags', [])
        if not isinstance(tags, list) or any(not isinstance(tag, str) for tag in tags):
            fail(f'{input_label}: tags must be an array of strings')
        return report_uuid, validated_dataset_uuids, tags, json_s_reads


    def load_integer_csv(path, columns, label):
        try:
            data = pd.read_csv(path, dtype=np.int64)
        except (OSError, ValueError) as exc:
            fail(f'{label}: could not read integer CSV: {exc}')
        if tuple(data.columns) != columns:
            fail(f'{label}: CSV columns are {tuple(data.columns)!r}; expected {columns!r}')
        if data.min().lt(0).any():
            fail(f'{label}: CSV values must be non-negative')
        return data


    def load_read_lengths(path, expected_rows, input_label):
        label = f'{input_label}: read-length'
        int32_max = np.iinfo(np.int32).max
        zmws = np.empty(expected_rows, dtype=np.int64)
        values = {column: np.empty(expected_rows, dtype=np.int32) for column in read_length_header[1:]}
        rows = 0
        negative = out_of_range = False
        try:
            with pd.read_csv(path, dtype=np.int64, chunksize=1_000_000) as chunks:
                for chunk in chunks:
                    if tuple(chunk.columns) != read_length_header:
                        fail(f'{label}: CSV columns are {tuple(chunk.columns)!r}; expected {read_length_header!r}')
                    start, rows = rows, rows + len(chunk)
                    negative = negative or bool(chunk.min().lt(0).any())
                    out_of_range = out_of_range or any(
                        chunk[column].min() <= 0 or chunk[column].max() > int32_max for column in values
                    )
                    if rows <= expected_rows and not negative and not out_of_range:
                        zmws[start:rows] = chunk['zmw'].to_numpy()
                        for column, array in values.items():
                            array[start:rows] = chunk[column].to_numpy()
        except (OSError, ValueError) as exc:
            fail(f'{label}: could not read integer CSV: {exc}')
        if negative:
            fail(f'{label}: CSV values must be non-negative')
        if rows != expected_rows:
            fail(f'{input_label}: read-length CSV contains {rows} rows; Skera summary reports {expected_rows} S-reads')
        if out_of_range:
            fail(f'{input_label}: read lengths and concatenation factors must be between 1 and {int32_max}')
        first_reads = ~pd.Series(zmws).duplicated().to_numpy()
        del zmws
        return (
            values['m_read_length'],
            values['original_read_length'][first_reads],
            values['concatenation_factor'][first_reads],
        )


    def read_adapter_names(path):
        names = []
        try:
            with open(path, 'r', encoding='utf-8') as fasta_fh:
                for line_number, line in enumerate(fasta_fh, start=1):
                    if not line.startswith('>'):
                        continue
                    name = line[1:].strip().split()[0] if line[1:].strip() else ''
                    if not name:
                        fail(f'segmentation adapter FASTA line {line_number} has an empty name')
                    names.append(name)
        except OSError as exc:
            fail(f'could not read segmentation adapter FASTA: {exc}')
        if not names:
            fail('segmentation adapter FASTA contains no records')
        if len(names) != len(set(names)):
            fail('segmentation adapter FASTA contains duplicate record names')
        return names


    def load_ligations(path, adapter_count, input_label):
        data = load_integer_csv(path, ligation_header, f'{input_label}: ligation')
        if data[['adapter_1', 'adapter_2']].max().ge(adapter_count).any():
            fail(f'{input_label}: ligation adapter pair is outside the adapter FASTA')
        if data.duplicated(subset=['adapter_1', 'adapter_2']).any():
            fail(f'{input_label}: duplicate ligation adapter pair')
        if len(data) != adapter_count * adapter_count:
            fail(f'{input_label}: ligation CSV is missing adapter pairs')
        counts = np.zeros((adapter_count, adapter_count), dtype=np.int64)
        counts[data['adapter_2'].to_numpy(), data['adapter_1'].to_numpy()] = data['ligations'].to_numpy()
        return counts


    class Palette:
        white = (1, 1, 1)
        orange = (225 / 255, 106 / 255, 44 / 255)
        pink = (223 / 255, 25 / 255, 149 / 255)
        purple = (95 / 255, 36 / 255, 159 / 255)
        blue = (19 / 255, 131 / 255, 198 / 255)
        teal = (0 / 255, 156 / 255, 162 / 255)
        green = (0 / 255, 157 / 255, 78 / 255)
        light_orange = (249 / 255, 157 / 255, 65 / 255)
        light_pink = (255 / 255, 102 / 255, 204 / 255)
        light_purple = (141 / 255, 109 / 255, 176 / 255)
        light_blue = (68 / 255, 168 / 255, 223 / 255)
        light_teal = (0 / 255, 184 / 255, 196 / 255)
        light_green = (106 / 255, 191 / 255, 105 / 255)
        dark_orange = (195 / 255, 74 / 255, 33 / 255)
        dark_blue = (0 / 255, 84 / 255, 150 / 255)
        dark_green = (0 / 255, 91 / 255, 66 / 255)
        dark_grey = (49 / 255, 53 / 255, 65 / 255)


    palette = Palette()
    clean_style = {
        'font.size': 8.0,
        'font.sans-serif': ['Roboto', 'DejaVu Sans'],
        'xtick.major.size': 2,
        'xtick.major.width': 0.75,
        'xtick.labelsize': 8.0,
        'xtick.direction': 'out',
        'ytick.major.size': 2,
        'ytick.major.width': 0.75,
        'ytick.labelsize': 8.0,
        'ytick.direction': 'out',
        'xtick.major.pad': 2,
        'xtick.minor.pad': 2,
        'ytick.major.pad': 2,
        'ytick.minor.pad': 2,
        'savefig.dpi': 600,
        'axes.linewidth': 0.75,
        'mathtext.default': 'regular',
    }


    def save_and_close(fig, png_file_name, dpi, barcode):
        thumb = f'{op.splitext(png_file_name)[0]}_thumb.png'
        fig.axes[0].set_title(barcode, fontsize=16)
        fig.savefig(png_file_name, bbox_inches='tight', dpi=dpi)
        fig.savefig(thumb, bbox_inches='tight', dpi=20)
        plt.close(fig)
        return op.basename(png_file_name), op.basename(thumb)


    def plot_concatemer_histogram(concat_factors, file_prefix, barcode):
        arraysize = int(np.max(concat_factors))
        heights = np.bincount(concat_factors - 1, minlength=arraysize).tolist()
        total = len(concat_factors)
        percents = [round(count / total * 100, 2) for count in heights]
        with plt.style.context(clean_style):
            fig = plt.figure(figsize=(9, 6))
            axes = fig.add_axes([0.15, 0.125, 0.8, 0.8])
            for index, height in enumerate(heights):
                left = index + 0.5
                axes.add_patch(Rect((left, 0), 1, height, linewidth=0, facecolor=palette.blue))
                axes.text(
                    index + 1,
                    height,
                    f'{percents[index]}%',
                    horizontalalignment='center',
                    verticalalignment='bottom',
                    fontsize='small',
                )
            axes.set_xlim(0.25, arraysize + 0.75)
            axes.set_ylim(0, max(heights) * 1.1)
            axes.set_xlabel('Array Length')
            axes.set_ylabel('Number of Reads')
            axes.set_xticks(range(1, arraysize + 1))
            axes.set_xticklabels(range(1, arraysize + 1))
            return save_and_close(fig, f'{file_prefix}.concatemer_histogram.png', 600, barcode)


    def plot_ligations_heatmap(counts, adapter_names, file_prefix, barcode):
        with plt.style.context(clean_style):
            fig = plt.figure(figsize=(9, 6))
            axes = fig.add_axes([0.1, 0.1, 4.25 / 6, 4.25 / 5])
            scale = fig.add_axes([0.125 + (4.25 / 6), 0.1, 0.06, 4.25 / 5])
            edge_len = len(counts)
            count_max = int(np.max(counts))
            if count_max == 0:
                normalized = np.zeros_like(counts, dtype=int)
            else:
                normalized = ((counts / count_max) * 100).astype(int)
            gradient = [
                palette.white,
                palette.light_purple,
                palette.teal,
                palette.light_orange,
                palette.blue,
            ]
            red = []
            green = []
            blue = []
            steps = 100 // (len(gradient) - 1)
            difference = 101 - (steps * (len(gradient) - 1))
            for index in range(len(gradient) - 1):
                segment_steps = steps + difference if index == len(gradient) - 2 else steps
                red.extend(np.linspace(gradient[index][0], gradient[index + 1][0], segment_steps))
                green.extend(np.linspace(gradient[index][1], gradient[index + 1][1], segment_steps))
                blue.extend(np.linspace(gradient[index][2], gradient[index + 1][2], segment_steps))
            digit_count = len(str(count_max))
            font_size = 5 if digit_count <= 5 else 4 if digit_count == 6 else 3
            for row in range(edge_len):
                for column in range(edge_len):
                    color_index = normalized[column][row]
                    axes.add_patch(
                        Rect(
                            (row, column),
                            1,
                            1,
                            linewidth=0,
                            zorder=9,
                            facecolor=(red[color_index], green[color_index], blue[color_index]),
                        )
                    )
                    axes.text(
                        row + 0.5,
                        column + 0.5,
                        str(counts[column][row]),
                        horizontalalignment='center',
                        verticalalignment='center',
                        fontsize=font_size,
                        zorder=10,
                    )
            for index in range(len(red)):
                scale.add_patch(Rect((0, index), 1, 1, linewidth=0, facecolor=(red[index], green[index], blue[index])))
            axes.set_xlim(0, edge_len)
            axes.set_ylim(0, edge_len)
            axes.set_xticks(np.arange(0.5, edge_len + 0.5))
            axes.set_yticks(np.arange(0.5, edge_len + 0.5))
            axes.set_xticklabels(adapter_names)
            axes.set_yticklabels(adapter_names)
            axes.set_xlabel('Adapter 1')
            axes.set_ylabel('Adapter 2')
            axes.tick_params(
                axis='both',
                which='both',
                bottom=True,
                labelbottom=True,
                left=True,
                labelleft=True,
                right=False,
                labelright=False,
                top=False,
                labeltop=False,
            )
            scale.set_xlim(0, 1)
            scale.set_ylim(0, 100)
            scale.set_yticks(range(0, 101, 20))
            scale.set_yticklabels([int(value) for value in np.linspace(0, count_max, 6)], fontsize=6)
            scale.tick_params(
                axis='both',
                which='both',
                bottom=False,
                labelbottom=False,
                left=False,
                labelleft=False,
                right=True,
                labelright=True,
                top=False,
                labeltop=False,
            )
            return save_and_close(fig, f'{file_prefix}.ligations_heatmap.png', 600, barcode)


    def plot_hifi_read_length(read_lengths, concat_factors, file_prefix, barcode, xmax=25000):
        binsize = 250
        bin_count = xmax // binsize
        arraysize = int(np.max(concat_factors))
        colors = [
            palette.light_orange,
            palette.light_pink,
            palette.light_purple,
            palette.light_blue,
            palette.light_teal,
            palette.light_green,
            palette.orange,
            palette.pink,
            palette.purple,
            palette.blue,
            palette.teal,
            palette.green,
            palette.dark_orange,
            palette.dark_blue,
            palette.dark_green,
            palette.dark_grey,
        ][:arraysize]
        if arraysize > len(colors):
            fail(f'concatenation factor {arraysize} exceeds the supported plot palette')
        bins = np.clip(read_lengths, 0, xmax - 1).astype(np.int64) // binsize
        per_bin_concat = np.bincount(
            bins * arraysize + (concat_factors.astype(np.int64) - 1),
            minlength=bin_count * arraysize,
        ).reshape(bin_count, arraysize)
        heights = per_bin_concat.sum(axis=1).tolist()
        with plt.style.context(clean_style):
            fig = plt.figure(figsize=(9, 6))
            axes = fig.add_axes([0.15, 0.125, 0.65, 0.8])
            legend = fig.add_axes([0.8125, 0.125, 0.1, 0.8], frameon=False)
            for index in range(bin_count):
                left = index * binsize
                if not heights[index]:
                    continue
                last_top = 0
                for concat_index, concat_count in enumerate(per_bin_concat[index]):
                    axes.add_patch(
                        Rect(
                            (left, last_top),
                            binsize,
                            concat_count,
                            linewidth=0.1,
                            edgecolor='black',
                            facecolor=colors[concat_index],
                        )
                    )
                    last_top += concat_count
            axes.set_xlim(0, xmax + binsize)
            axes.set_ylim(0, max(heights) * 1.1)
            axes.set_xlabel('Read length, bp')
            axes.set_ylabel('Number of Reads')
            axes.set_xticks(range(0, xmax + binsize, binsize * 10))
            legend_top = 0
            for index, color in enumerate(colors):
                bottom = index + 0.2 + (0.1 * index)
                legend.add_patch(Rect((0.2, bottom), 1, 1, linewidth=0, color=color))
                legend.text(
                    1.4,
                    bottom + 0.5,
                    f'{index + 1}x',
                    horizontalalignment='left',
                    verticalalignment='center',
                )
                legend_top = bottom + 1.2
            legend.set_xlim(0, 2)
            legend.set_ylim(0, legend_top)
            legend.tick_params(
                axis='both',
                which='both',
                bottom=False,
                labelbottom=False,
                left=False,
                labelleft=False,
                right=False,
                labelright=False,
                top=False,
                labeltop=False,
            )
            return save_and_close(fig, f'{file_prefix}.hifi_read_length.png', 300, barcode)


    def bounded_read_length_max(read_lengths, maximum=75000):
        rank = int(0.9999 * len(read_lengths))
        upper = min(
            int(1.05 * np.max(read_lengths)),
            int(1.05 * np.partition(read_lengths, rank)[rank]),
        )
        return int(min(maximum, upper)) if upper > 0 else maximum


    def plot_sread_length(read_lengths, file_prefix, barcode):
        read_length_max = bounded_read_length_max(read_lengths)
        fig = plt.figure(figsize=(8, 6))
        axes = fig.add_subplot(111)
        axes.patch.set_facecolor('#ffffff')
        axes.grid(color='#e0e0e0', linewidth=0.5, linestyle='-')
        axes.tick_params(axis='both', which='minor', labelsize=12)
        axes.tick_params(axis='both', which='major', labelsize=12)
        axes.set_axisbelow(True)
        axes.hist(
            [read_lengths],
            bins=50,
            color=['#008eb9'],
            edgecolor=(1, 1, 1, 0),
            label=['HiFi reads'],
            stacked=True,
            range=(0, read_length_max),
            rwidth=0.8,
            linewidth=0.5,
        )
        axes.set_xlabel('S-read length, bp', fontsize=16)
        axes.set_ylabel('Number of S-reads', fontsize=16)
        axes.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
        axes.xaxis.set_major_formatter(ticker.FuncFormatter(lambda value, _position: format(int(value), ',')))
        axes.yaxis.set_major_formatter(ticker.FuncFormatter(lambda value, _position: format(int(value), ',')))
        return save_and_close(fig, f'{file_prefix}.sread_length.png', 350, barcode)


    expected_counts = [int(value) for value in expected_counts_text]
    manifests = [
        read_manifest(hifi_barcodes_manifest, expected_counts[0], 'hifi_barcodes'),
        read_manifest(skera_summary_csvs_manifest, expected_counts[1], 'skera_summary_csvs'),
        read_manifest(skera_summary_jsons_manifest, expected_counts[2], 'skera_summary_jsons'),
        read_manifest(read_lengths_csvs_manifest, expected_counts[3], 'read_lengths_csvs'),
        read_manifest(ligations_csvs_manifest, expected_counts[4], 'ligations_csvs'),
    ]
    manifest_names = (
        'hifi_barcodes',
        'skera_summary_csvs',
        'skera_summary_jsons',
        'read_lengths_csvs',
        'ligations_csvs',
    )
    lengths = [len(values) for values in manifests]
    if any(length == 0 for length in lengths):
        fail(f'aligned report input arrays must be non-empty; observed {dict(zip(manifest_names, lengths))!r}')
    if len(set(lengths)) != 1:
        fail(f'aligned report input arrays have unequal lengths: {dict(zip(manifest_names, lengths))!r}')

    hifi_barcodes, summary_csv_paths, summary_json_paths, read_length_paths, ligation_paths = manifests
    adapter_names = read_adapter_names(segmentation_adapters)
    input_plot_data = []
    input_metrics = []
    dataset_uuids = set()
    report_tags = set()
    uuid_inputs = []

    for input_index, artifacts in enumerate(
        zip(
            hifi_barcodes,
            summary_csv_paths,
            summary_json_paths,
            read_length_paths,
            ligation_paths,
        ),
        start=1,
    ):
        barcode, summary_csv, summary_json, read_lengths_csv, ligations_csv = artifacts
        if not barcode:
            fail(f'input {input_index}: hifi_barcode must be a non-empty string')
        input_label = f'input {input_index} ({barcode})'
        metrics = load_summary_csv(summary_csv, input_label)
        report_uuid, source_dataset_uuids, source_tags, json_s_reads = load_skera_json(summary_json, input_label)
        if json_s_reads != metrics['s_reads']:
            fail(
                f'{input_label}: Skera summary JSON s_reads {json_s_reads} does not match '
                f'Skera summary CSV Segmented Reads (S-Reads) {metrics["s_reads"]}'
            )
        sread_lengths, parent_read_lengths, parent_concat_factors = load_read_lengths(
            read_lengths_csv,
            metrics['s_reads'],
            input_label,
        )
        ligations = load_ligations(ligations_csv, len(adapter_names), input_label)
        input_metrics.append({'barcode': barcode, **metrics})
        input_plot_data.append(
            {
                'barcode': barcode,
                'sread_lengths': sread_lengths,
                'parent_read_lengths': parent_read_lengths,
                'parent_concat_factors': parent_concat_factors,
                'ligations': ligations,
            }
        )
        dataset_uuids.update(source_dataset_uuids)
        report_tags.update(source_tags)
        uuid_inputs.append({'hifi_barcode': barcode, 'report_uuid': report_uuid})

    grouped = OrderedDict()
    for metrics in input_metrics:
        group = grouped.setdefault(metrics['barcode'], [])
        group.append(metrics)
    grouped_plot_data = OrderedDict()
    for plot_data in input_plot_data:
        barcode_plot_data = grouped_plot_data.setdefault(
            plot_data['barcode'],
            {
                'sread_length_arrays': [],
                'parent_read_length_arrays': [],
                'parent_concat_factor_arrays': [],
                'ligations': np.zeros((len(adapter_names), len(adapter_names)), dtype=np.int64),
            },
        )
        barcode_plot_data['sread_length_arrays'].append(plot_data['sread_lengths'])
        barcode_plot_data['parent_read_length_arrays'].append(plot_data['parent_read_lengths'])
        barcode_plot_data['parent_concat_factor_arrays'].append(plot_data['parent_concat_factors'])
        barcode_plot_data['ligations'] += plot_data['ligations']


    grouped_rows = []
    for barcode, group in grouped.items():
        if len(group) == 1:
            csv_values = group[0]['raw_values']
        else:
            reads = sum(item['reads'] for item in group)
            s_reads = sum(item['s_reads'] for item in group)
            csv_values = [
                str(reads),
                str(s_reads),
                stable_calculation(
                    sum(item['mean_len_s_reads'] * item['s_reads'] for item in group),
                    s_reads,
                ),
                stable_calculation(
                    sum(item['percent_full_array'] * item['reads'] for item in group),
                    reads,
                ),
                stable_calculation(
                    sum(item['mean_array_size'] * item['reads'] for item in group),
                    reads,
                ),
            ]
        grouped_rows.append(
            {
                'barcode': barcode,
                'csv_values': csv_values,
                'numeric_values': [numeric_csv_value(value) for value in csv_values],
            }
        )

    with open(output_csv, 'w', encoding='utf-8', newline='') as output_csv_fh:
        writer = csv.writer(output_csv_fh, lineterminator='\n')
        writer.writerow(csv_header)
        for row in grouped_rows:
            writer.writerow((row['barcode'], *row['csv_values']))

    total_reads = sum(item['reads'] for item in input_metrics)
    total_s_reads = sum(item['s_reads'] for item in input_metrics)
    overall_value_strings = (
        str(total_reads),
        str(total_s_reads),
        stable_calculation(
            sum(item['mean_len_s_reads'] * item['s_reads'] for item in input_metrics),
            total_s_reads,
        ),
        stable_calculation(
            sum(item['percent_full_array'] * item['reads'] for item in input_metrics),
            total_reads,
        ),
        stable_calculation(
            sum(item['mean_array_size'] * item['reads'] for item in input_metrics),
            total_reads,
        ),
    )
    overall_values = [numeric_csv_value(value) for value in overall_value_strings]

    emitted_plots = OrderedDict()
    for barcode_index, (barcode, plot_data) in enumerate(grouped_plot_data.items(), start=1):
        sread_lengths = np.concatenate(plot_data['sread_length_arrays'])
        if len(sread_lengths) == 0:
            continue
        parent_read_lengths = np.concatenate(plot_data['parent_read_length_arrays'])
        parent_concat_factors = np.concatenate(plot_data['parent_concat_factor_arrays'])
        file_prefix = f'barcode_{barcode_index}'
        emitted_plots[barcode] = {
            'barcode_index': barcode_index,
            'plots': {
                'concatemer_histogram': plot_concatemer_histogram(parent_concat_factors, file_prefix, barcode),
                'ligations_heatmap': plot_ligations_heatmap(
                    plot_data['ligations'],
                    adapter_names,
                    file_prefix,
                    barcode,
                ),
                'hifi_read_length': plot_hifi_read_length(
                    parent_read_lengths,
                    parent_concat_factors,
                    file_prefix,
                    barcode,
                ),
                'sread_length': plot_sread_length(sread_lengths, file_prefix, barcode),
            },
        }

    attributes = [
        {
            'id': f'{report_id}.{metric_id}',
            'name': name,
            'value': value,
        }
        for (name, metric_id), value in zip(
            (
                ('Reads', 'reads'),
                ('Segmented reads (S-reads)', 's_reads'),
                ('Mean length of S-reads', 'mean_len_s_reads'),
                ('Percent of reads with full arrays', 'percent_full_array'),
                ('Mean array size (concatenation factor)', 'mean_array_size'),
            ),
            overall_values,
        )
    ]
    table_column_specs = (
        ('hifi_barcode', 'HiFi barcode'),
        ('reads', 'Input Reads'),
        ('s_reads', 'Segmented Reads (S-Reads)'),
        ('mean_len_s_reads', 'Mean Length of S-Reads'),
        ('percent_full_array', 'Percentage of Reads with Full Array'),
        ('mean_array_size', 'Mean Array Size (Concatenation Factor)'),
    )
    table_columns = []
    for column_index, (column_id, header) in enumerate(table_column_specs):
        if column_index == 0:
            values = [row['barcode'] for row in grouped_rows]
        else:
            values = [row['numeric_values'][column_index - 1] for row in grouped_rows]
            if column_id in ('reads', 's_reads', 'mean_len_s_reads'):
                values = [f'{value:,}' for value in values]
            elif column_id == 'percent_full_array':
                values = [f'{value:.2f}'.rstrip('0').rstrip('.') + '%' for value in values]
            elif column_id == 'mean_array_size':
                values = [f'{value:.2f}'.rstrip('0').rstrip('.') for value in values]
        table_columns.append(
            {
                'id': f'{report_id}.segmentation_summary.{column_id}',
                'header': header,
                'values': values,
            }
        )

    plot_groups = []
    if emitted_plots:
        for plot_group_id, group_title, group_plots in plot_specs:
            plots = []
            for barcode, emitted_barcode_plots in emitted_plots.items():
                barcode_index = emitted_barcode_plots['barcode_index']
                barcode_plots = emitted_barcode_plots['plots']
                for plot_id, title, caption in group_plots:
                    image, thumbnail = barcode_plots[plot_id]
                    plots.append(
                        {
                            'id': f'{report_id}.{plot_group_id}.{plot_id}_barcode_{barcode_index}',
                            'title': f'{title} ({barcode})',
                            'caption': f'{caption} for HiFi barcode {barcode}',
                            'plotType': 'image',
                            'image': image,
                            'thumbnail': thumbnail,
                        }
                    )
            plot_groups.append(
                {
                    'id': f'{report_id}.{plot_group_id}',
                    'title': group_title,
                    'legend': None,
                    'thumbnail': None,
                    'plots': plots,
                }
            )

    uuid_seed = json.dumps(uuid_inputs, sort_keys=True, separators=(',', ':'))
    report_uuid = str(
        uuid.uuid5(
            uuid.NAMESPACE_URL,
            f'kinnex-isoseq-wdl/read_segmentation/{uuid_seed}',
        )
    )
    report = {
        '_comment': f'Combined from {len(input_metrics)} Skera summary reports',
        'id': report_id,
        'version': report_version,
        'uuid': report_uuid,
        'title': report_title,
        'attributes': attributes,
        'tables': [
            {
                'id': f'{report_id}.segmentation_summary',
                'title': 'Segmentation Summary',
                'columns': table_columns,
            }
        ],
        'plotGroups': plot_groups,
        'dataset_uuids': sorted(dataset_uuids),
        'tags': sorted(report_tags),
    }
    with open(output_json, 'w', encoding='utf-8') as output_json_fh:
        json.dump(report, output_json_fh, indent=2, sort_keys=True)
        output_json_fh.write('\n')
    PY
  >>>

  output {
    File segmentation_summary_csv = "segmentation.summary.csv"
    File segmentation_summary_report = "segmentation.summary.report.json"
    Array[File] segmentation_summary_report_resources = glob("*.png")
  }

  runtime {
    cpu: threads
    memory: total_mem_mb + " MB"
    docker: runtime_attributes.container_registry + "/pb_wdl_base@sha256:03cb3c01937eccc907f8ad71c87b258581504572205fe3f31a657e318f3564ae"  # 1.1.0
    maxRetries: runtime_attributes.max_retries
  }
}
