version 1.0

import "../../rna_structs.wdl"

task pigeon_gene_saturation {
  meta {
    description: "Subsample pooled supporting reads from the filtered joint classification with Pigeon's native sampling rules and singleton exclusion."
    outputs: {
      pigeon_saturation: {
        description: "Native filtered joint gene-saturation trajectory"
      }
    }
  }

  parameter_meta {
    pigeon_filtered_classification: {
      description: "Filtered Pigeon classification table with sample FL support"
    }
    supporting_reads: {
      description: "Pooled supporting reads in the filtered classification, used to size memory"
    }
    output_prefix: {
      description: "Output prefix"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    File pigeon_filtered_classification
    Int supporting_reads
    String output_prefix = "joint"
    RuntimeAttributes runtime_attributes
  }

  Int threads = 1
  # Estimate 1 GB per 4M reads for Pigeon (96 bytes/read + allocation peak and string heap).
  Int mem_gb = 4 + ceil(supporting_reads / 4000000.0)
  Int total_mem_mb = (mem_gb * 1024) + runtime_attributes.add_memory_mb

  command <<<
    set -euo pipefail

    pigeon report \
      --log-level INFO \
      --log-file "~{output_prefix}.pigeon.saturation.log" \
      --exclude-singletons \
      "~{pigeon_filtered_classification}" \
      "~{output_prefix}.pigeon.saturation.tsv"
  >>>

  output {
    File pigeon_saturation = "~{output_prefix}.pigeon.saturation.tsv"
  }

  runtime {
    cpu: threads
    memory: total_mem_mb + " MB"
    docker: runtime_attributes.container_registry + "/pigeon@sha256:94cf3cd64f600ae974a7956177be2231c1c9ecb468a1cc4b70b05dcacf59a81c"  # 26.2.0_build2
    maxRetries: runtime_attributes.max_retries
  }
}

task generate_pigeon_report {
  meta {
    description: "Generate the joint Transcript Mapping and Classification report from native Pigeon metrics and classification statistics."
    outputs: {
      pigeon_summary_report: {
        description: "Transcript Mapping and Classification JSON report"
      },
      pigeon_summary_report_resources: {
        description: "Six report PNGs and their six thumbnails"
      }
    }
  }

  parameter_meta {
    pigeon_classify_report: {
      description: "Native Pigeon per-sample classification metrics JSON"
    }
    pigeon_filter_report: {
      description: "Native Pigeon per-sample filtered classification metrics JSON"
    }
    pigeon_classification: {
      description: "Unfiltered joint Pigeon classification table"
    }
    pigeon_filtered_classification: {
      description: "Filtered joint Pigeon classification table"
    }
    pigeon_saturation: {
      description: "Native filtered joint gene-saturation trajectory"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    File pigeon_classify_report
    File pigeon_filter_report
    File pigeon_classification
    File pigeon_filtered_classification
    File pigeon_saturation
    RuntimeAttributes runtime_attributes
  }

  Int threads = 1
  Int mem_gb = 4
  Int total_mem_mb = (mem_gb * 1024) + runtime_attributes.add_memory_mb

  command <<<
    set -euo pipefail

    python3 - \
      "~{pigeon_classify_report}" \
      "~{pigeon_filter_report}" \
      "~{pigeon_classification}" \
      "~{pigeon_filtered_classification}" \
      "~{pigeon_saturation}" <<'PY'
    # Upstream sources: internal/docs/upstream-sources.md#pigeon-transcript-mapping-and-classification-report
    # Copyright (c) 2016, Pacific Biosciences. All rights reserved.
    # The Clear BSD terms are provided in this repository's LICENSE.
    import csv
    import gzip
    import json
    import sys
    import uuid

    REPORT_ID = 'isoseq_mapping'
    # Table order follows pigeon.py; plot order deliberately differs.
    CATEGORIES = (
        ('full-splice_match', 'FSM', '#1383c6', '#44a8df'),
        ('incomplete-splice_match', 'ISM', '#e16a2c', '#f99d41'),
        ('novel_in_catalog', 'NIC', '#009d4e', '#6abf69'),
        ('novel_not_in_catalog', 'NNC', '#df1995', '#ff66cc'),
        ('antisense', 'Antisense', '#009ca2', '#00b8c4'),
        ('fusion', 'Fusion', '#5f249f', '#8d6db0'),
        ('moreJunctions', 'More junctions', '#005b42', '#009d4e'),
        ('genic_intron', 'Genic intron', '#005496', '#1383c6'),
        ('genic', 'Genic genomic', '#676f7f', '#9ba7ba'),
        ('intergenic', 'Intergenic', '#c34a21', '#e16a2c'),
    )
    PLOT_ORDER = (0, 1, 2, 3, 8, 4, 5, 9, 7, 6)
    SUMMARY_COLUMNS = (
        ('sample_name', 'Sample Name', None, 'sample_name'),
        ('total_unique_genes', 'Genes', 'classify', 'total_unique_genes'),
        ('total_unique_genes_filtered', 'Genes, filtered', 'filter', 'total_unique_genes'),
        ('total_unique_genes_known', 'Known Genes, filtered', 'filter', 'total_unique_genes_known'),
        ('total_unique_transcripts', 'Isoforms', 'classify', 'total_unique_transcripts'),
        ('total_unique_transcripts_filtered', 'Isoforms, filtered', 'filter', 'total_unique_transcripts'),
        ('total_unique_transcripts_known', 'Known isoforms, filtered', 'filter', 'total_unique_transcripts_known'),
        (
            'total_unique_transcripts_novel_greater_1tpm',
            'Novel isoforms > 1TPM, filtered',
            'filter',
            'total_unique_transcripts_novel_greater_1tpm',
        ),
    )
    CLASS_COLUMNS = (
        ('category', 'Category'),
        ('count', 'Count'),
        ('cage_detected', 'CAGE Detected'),
        ('cage_detected_pct', 'CAGE Detected (%)'),
        ('polya_detected', 'polyA Motif Detected'),
        ('polya_detected_pct', 'polyA Motif Detected (%)'),
    )
    SUPPORT_DESCRIPTION = (
        'CAGE detection is within_CAGE_peak=TRUE (case-insensitive); polyA motif '
        'detection is polyA_motif other than NA. Missing or disabled support '
        'resources yield NA and no detections, not evidence of assessed absence. '
        'Percentages are fractions of joint isoforms in each category.'
    )


    def fail(message):
        raise ValueError('Pigeon report contract error: ' + message)


    def leaf_id(value, label):
        if not isinstance(value, str) or not value:
            fail(f'{label}: missing metric ID')
        return value.rsplit('.', 1)[-1]


    def metric_map(items, value_key, label):
        if not isinstance(items, list):
            fail(f'{label}: expected an array')
        result = {}
        for item in items:
            if not isinstance(item, dict):
                fail(f'{label}: expected metric objects')
            metric = leaf_id(item.get('id'), label)
            if metric in result:
                fail(f'{label}: duplicate metric {metric!r}')
            if value_key not in item:
                fail(f'{label}: missing {value_key} for {metric!r}')
            result[metric] = item[value_key]
        return result


    def native_metrics(path, label, required):
        with open(path, encoding='utf-8') as handle:
            report = json.load(handle)
        if not isinstance(report, dict):
            fail(f'{label}: report must be an object')
        tables = report.get('tables', [])
        if not isinstance(tables, list):
            fail(f'{label}: tables must be an array')
        sample_tables = [
            table
            for table in tables
            if isinstance(table, dict) and leaf_id(table.get('id'), label) == 'classification_summary_by_sample'
        ]
        if len(sample_tables) != 1:
            fail(f'{label}: expected exactly one classification_summary_by_sample table')
        columns = metric_map(sample_tables[0].get('columns'), 'values', label)
        for metric in ('sample_name', *required):
            if metric not in columns:
                fail(f'{label}: missing required metric {metric!r}')
        samples = columns['sample_name']
        if not isinstance(samples, list) or not samples:
            fail(f'{label}: sample_name must be a nonempty array')
        for metric, values in columns.items():
            if not isinstance(values, list) or len(values) != len(samples):
                fail(f'{label}: {metric!r} column length does not match sample_name')
        rows = {}
        for index, sample in enumerate(samples):
            if not isinstance(sample, str) or not sample:
                fail(f'{label}: sample_name must contain nonempty strings')
            if sample in rows:
                fail(f'{label}: duplicate sample_name {sample!r}')
            row = {}
            for metric in required:
                value = columns[metric][index]
                if isinstance(value, bool) or not isinstance(value, int) or value < 0:
                    fail(f'{label}: {metric!r} for {sample!r} must be a nonnegative integer')
                row[metric] = value
            rows[sample] = row
        return samples, rows


    def summary_table(classify_path, filter_path):
        inputs = {}
        for source, path in (('classify', classify_path), ('filter', filter_path)):
            required = tuple(dict.fromkeys(metric for _, _, origin, metric in SUMMARY_COLUMNS if origin == source))
            inputs[source] = native_metrics(path, source, required)
        samples, classify = inputs['classify']
        _, filtered = inputs['filter']
        if classify.keys() != filtered.keys():
            fail(
                'classify/filter sample_name sets differ; '
                f'missing from filter={sorted(classify.keys() - filtered.keys())!r}, '
                f'extra in filter={sorted(filtered.keys() - classify.keys())!r}'
            )
        columns = []
        for column, header, source, metric in SUMMARY_COLUMNS:
            values = list(samples) if source is None else [inputs[source][1][sample][metric] for sample in samples]
            columns.append(
                {
                    'id': f'{REPORT_ID}.summary_table.{column}',
                    'header': header,
                    'values': values,
                }
            )
        return {
            'id': f'{REPORT_ID}.summary_table',
            'title': 'Classification Summary Metrics (All Samples)',
            'columns': columns,
        }


    def tsv_rows(path, required):
        opener = gzip.open if path.endswith('.gz') else open
        with opener(path, 'rt', encoding='utf-8', newline='') as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            if reader.fieldnames is None:
                fail(f'{path}: missing TSV header')
            # Structural fields are case-insensitive; native sample identities are not.
            fields = [field.lower() if field.lower() in required else field for field in reader.fieldnames]
            if len(fields) != len(set(fields)):
                fail(f'{path}: duplicate TSV columns')
            reader.fieldnames = fields
            missing = set(required) - set(fields)
            if missing:
                fail(f'{path}: missing TSV columns {sorted(missing)!r}')
            for row in reader:
                if None in row or any(row[field] is None for field in required):
                    fail(f'{path}:{reader.line_num}: row does not match TSV header')
                yield reader.line_num, row


    def tsv_integer(value, label):
        try:
            number = int(value)
        except ValueError:
            fail(f'{label}: expected a nonnegative integer, got {value!r}')
        if number < 0:
            fail(f'{label}: expected a nonnegative integer, got {value!r}')
        return number


    def classification_stats(path):
        # Fixed-size aggregates: ten categories and nine 1-kb/overflow bins.
        counts = [[0, 0, 0, 0] for _ in CATEGORIES]
        lengths = [[0] * len(CATEGORIES) for _ in range(9)]
        category_index = {category[0]: index for index, category in enumerate(CATEGORIES)}
        required = ('structural_category', 'coding', 'length', 'within_cage_peak', 'polya_motif')
        for line, row in tsv_rows(path, required):
            category = row['structural_category']
            if category not in category_index:
                fail(f'{path}:{line}: unknown structural_category {category!r}')
            index = category_index[category]
            length = tsv_integer(row['length'], f'{path}:{line}: length')
            counts[index][0] += 1
            counts[index][1] += row['coding'] == 'coding'
            counts[index][2] += row['within_cage_peak'].lower() == 'true'
            counts[index][3] += row['polya_motif'].lower() != 'na'
            lengths[min(length // 1000, 8)][index] += 1
        return {'counts': counts, 'lengths': lengths}


    def classification_table(table_id, title, stats):
        rows = []
        for category, (total, coding, cage, polya) in zip(CATEGORIES, stats['counts']):
            rows.append(
                (
                    category[1],
                    total,
                    cage,
                    cage / total if total else 0.0,
                    polya,
                    polya / total if total else 0.0,
                )
            )
        return {
            'id': f'{REPORT_ID}.{table_id}',
            'title': title + ' (All Samples)',
            'columns': [
                {
                    'id': f'{REPORT_ID}.{table_id}.{column}',
                    'header': header,
                    'values': [row[index] for row in rows],
                }
                for index, (column, header) in enumerate(CLASS_COLUMNS)
            ],
        }


    def saturation_stats(path):
        series = {'reads': [], 'unique_genes': [], 'unique_genes_known': []}
        for line, row in tsv_rows(path, tuple(series)):
            for metric in series:
                series[metric].append(tsv_integer(row[metric], f'{path}:{line}: {metric}'))
        if not series['reads']:
            fail(f'{path}: missing native saturation final snapshot')
        return series


    def plotting():
        import matplotlib

        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib import ticker

        plt.rcParams.update(
            {
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
                'axes.linewidth': 0.75,
                'mathtext.default': 'regular',
            }
        )
        return plt, ticker


    def save_plot(plt, figure, filename):
        thumbnail = filename[:-4] + '_thumb.png'
        figure.savefig(filename, dpi=300, bbox_inches='tight')
        figure.savefig(thumbnail, dpi=60, bbox_inches='tight')
        plt.close(figure)
        return filename, thumbnail


    def plot_object(group, leaf, title, caption, files):
        image, thumbnail = files
        return {
            'id': f'{REPORT_ID}.{group}.{leaf}',
            'title': title,
            'caption': caption,
            'plotType': 'image',
            'image': image,
            'thumbnail': thumbnail,
        }


    def empty_state(axis, message):
        axis.text(
            0.5,
            0.5,
            message,
            transform=axis.transAxes,
            ha='center',
            va='center',
            color='#676f7f',
            bbox={'facecolor': 'white', 'edgecolor': 'none', 'alpha': 0.9},
        )


    def classification_plots(plt, ticker, group, title, stats, filtered):
        total = sum(count[0] for count in stats['counts'])
        suffix = ', filtered' if filtered else ''
        dist_title = f'Isoform Distribution Across Categories{suffix} (All Samples)'
        figure, axis = plt.subplots(figsize=(8, 4.6))
        heights = [stats['counts'][index][0] * 100 / total if total else 0.0 for index in PLOT_ORDER]
        label_offset = max(1, max(heights)) * 0.025
        for position, (index, height) in enumerate(zip(PLOT_ORDER, heights)):
            category = CATEGORIES[index]
            coding = stats['counts'][index][1]
            axis.bar(position, height, width=1, color=category[3], linewidth=0)
            axis.bar(position, coding * 100 / total if total else 0, width=1, color=category[2], linewidth=0)
            if total:
                axis.text(position, height + label_offset, f'{height:.2f}%', ha='center', va='bottom', fontsize=7)
        axis.set_xticks(range(10))
        axis.set_xticklabels(
            [CATEGORIES[index][1] for index in PLOT_ORDER], rotation=40, ha='right', rotation_mode='anchor'
        )
        axis.set_xlim(-0.6, 9.6)
        axis.set_ylim(0, max(1, max(heights)) * 1.2)
        axis.set_ylabel('Percentage of joint isoforms')
        axis.set_title(dist_title, pad=12)
        # The upstream overlay uses the category primary color for coding
        # and the lighter secondary color for noncoding isoforms.
        if not total:
            empty_state(axis, 'No classified isoforms')
        figure.tight_layout()
        dist = plot_object(
            group,
            'isoform_distribution_plot',
            dist_title,
            'Each joint isoform contributes once. Bar heights are percentages '
            'of all joint isoforms; darker category colors overlay coding '
            'isoforms, with lighter colors showing noncoding isoforms.',
            save_plot(plt, figure, group + '.iso_dist.png'),
        )

        length_title = f'Structural Categories by Transcript Length{suffix} (All Samples)'
        figure, axis = plt.subplots(figsize=(8, 4.6))
        bottoms = [0] * 9
        for index in reversed(PLOT_ORDER):
            heights = [row[index] for row in stats['lengths']]
            axis.bar(
                range(9),
                heights,
                bottom=bottoms,
                width=1,
                color=CATEGORIES[index][2],
                linewidth=0,
                label=CATEGORIES[index][1],
            )
            bottoms = [bottom + height for bottom, height in zip(bottoms, heights)]
        axis.set_xlim(-0.6, 8.6)
        axis.set_ylim(0, max(1, max(bottoms)) * 1.15)
        axis.set_xticks(range(9))
        axis.set_xticklabels([str(index) for index in range(8)] + ['≥8'])
        axis.set_xlabel('Isoform length (kb; bin lower bound)')
        axis.set_ylabel('Number of joint isoforms')
        axis.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
        axis.yaxis.set_major_formatter(ticker.StrMethodFormatter('{x:,.0f}'))
        axis.set_title(length_title, pad=12)
        handles, labels = axis.get_legend_handles_labels()
        axis.legend(handles[::-1], labels[::-1], loc='upper left', bbox_to_anchor=(1.01, 1), frameon=False)
        if not total:
            empty_state(axis, 'No classified isoforms')
        figure.tight_layout()
        lengths = plot_object(
            group,
            'structural_categories_plot',
            length_title,
            'Joint isoform counts stacked by structural category in 1-kb bins '
            '[lower bound, lower bound + 1 kb). The ≥8 bin includes every '
            'isoform at least 8 kb long. ' + SUPPORT_DESCRIPTION,
            save_plot(plt, figure, group + '.iso_lens.png'),
        )
        return {
            'id': f'{REPORT_ID}.{group}',
            'title': title + ' (All Samples)',
            'legend': None,
            'thumbnail': None,
            'plots': [dist, lengths],
        }


    def saturation_plots(plt, ticker, series):
        group = 'gene_saturation_plotgroup'
        support_policy = (
            'Singletons excluded: a transcript identity contributes only after more than one sampled supporting read.'
        )
        plots = []
        for index, (metric, label) in enumerate(
            (
                ('unique_genes', 'all genes'),
                ('unique_genes_known', 'known genes only'),
            )
        ):
            title = f'Gene Saturation, {label}, filtered (All Samples)'
            figure, axis = plt.subplots(figsize=(6.4, 4))
            x, y = series['reads'], series[metric]
            if x[-1] == 0:
                empty_state(axis, 'No supporting reads')
            else:
                # Unclipped above the spines so an all-zero series stays visible on the x-axis.
                # A line through one distinct point, including Pigeon's duplicate final
                # snapshot at an increment boundary, has no length and needs a marker.
                axis.plot(
                    x,
                    y,
                    color='#df1995',
                    linewidth=1,
                    marker='o' if len(set(zip(x, y))) == 1 else None,
                    markersize=5,
                    clip_on=False,
                    zorder=3,
                )
            axis.set_xlim(0, max(1, max(x)) * 1.05)
            axis.set_ylim(0, max(1, max(y)) * 1.1)
            axis.set_xlabel('Total sampled supporting reads')
            axis.set_ylabel('Total number of unique genes')
            axis.set_title(title, pad=12)
            for coordinate in (axis.xaxis, axis.yaxis):
                coordinate.set_major_locator(ticker.MaxNLocator(integer=True))
                coordinate.set_major_formatter(ticker.StrMethodFormatter('{x:,.0f}'))
            plt.setp(axis.get_xticklabels(), rotation=45, ha='right', rotation_mode='anchor')
            figure.tight_layout()
            plots.append(
                plot_object(
                    group,
                    f'gene_saturation_plot_{index}',
                    title,
                    'Native Pigeon subsampling of pooled supporting reads from the '
                    'filtered joint classification; all sampled reads remain on '
                    'the x-axis. ' + support_policy,
                    save_plot(plt, figure, f'gene_saturation_plot-{index}.png'),
                )
            )
        return {
            'id': f'{REPORT_ID}.{group}',
            'title': 'Gene Saturation (All Samples)',
            'legend': None,
            'thumbnail': None,
            'plots': plots,
        }


    def main(argv):
        if len(argv) != 5:
            fail('expected classify JSON, filter JSON, unfiltered TSV, filtered TSV, saturation TSV')
        classify_path, filter_path, unfiltered_path, filtered_path, saturation_path = argv
        summary = summary_table(classify_path, filter_path)
        unfiltered = classification_stats(unfiltered_path)
        filtered = classification_stats(filtered_path)
        saturation = saturation_stats(saturation_path)
        tables = [
            summary,
            classification_table('transcript_classification', 'Transcript Classification', unfiltered),
            classification_table('transcript_classification_filtered', 'Transcript Classification, filtered', filtered),
        ]
        # Canonical report values and complete plot aggregates, never file paths,
        # source UUIDs, timestamps, or image bytes, determine the report identity.
        seed = json.dumps(
            {
                'tables': tables,
                'unfiltered': unfiltered,
                'filtered': filtered,
                'saturation': saturation,
                'exclude_singletons': True,
            },
            sort_keys=True,
            separators=(',', ':'),
            allow_nan=False,
        )
        plt, ticker = plotting()
        groups = [
            classification_plots(
                plt, ticker, 'transcript_classification_plotgroup', 'Transcript Classification Plots', unfiltered, False
            ),
            classification_plots(
                plt,
                ticker,
                'transcript_classification_filtered_plotgroup',
                'Transcript Classification Plots, filtered',
                filtered,
                True,
            ),
            saturation_plots(plt, ticker, saturation),
        ]
        report = {
            '_comment': SUPPORT_DESCRIPTION,
            'id': REPORT_ID,
            'title': 'Transcript Mapping and Classification',
            'version': '1.0.1',
            'uuid': str(uuid.uuid5(uuid.NAMESPACE_URL, f'kinnex-isoseq-wdl/{REPORT_ID}/{seed}')),
            'attributes': [],
            'tables': tables,
            'plotGroups': groups,
            'dataset_uuids': [],
            'tags': [],
        }
        with open('isoseq_mapping.report.json', 'w', encoding='utf-8') as handle:
            json.dump(report, handle, indent=2, sort_keys=True, allow_nan=False)
            handle.write('\n')


    if __name__ == '__main__':
        try:
            main(sys.argv[1:])
        except ValueError as error:
            sys.exit(str(error))
    PY
  >>>

  output {
    File pigeon_summary_report = "isoseq_mapping.report.json"
    Array[File] pigeon_summary_report_resources = [
      "transcript_classification_plotgroup.iso_dist.png",
      "transcript_classification_plotgroup.iso_dist_thumb.png",
      "transcript_classification_plotgroup.iso_lens.png",
      "transcript_classification_plotgroup.iso_lens_thumb.png",
      "transcript_classification_filtered_plotgroup.iso_dist.png",
      "transcript_classification_filtered_plotgroup.iso_dist_thumb.png",
      "transcript_classification_filtered_plotgroup.iso_lens.png",
      "transcript_classification_filtered_plotgroup.iso_lens_thumb.png",
      "gene_saturation_plot-0.png",
      "gene_saturation_plot-0_thumb.png",
      "gene_saturation_plot-1.png",
      "gene_saturation_plot-1_thumb.png"
    ]
  }

  runtime {
    cpu: threads
    memory: total_mem_mb + " MB"
    docker: runtime_attributes.container_registry + "/pb_wdl_base@sha256:03cb3c01937eccc907f8ad71c87b258581504572205fe3f31a657e318f3564ae"  # 1.1.0
    maxRetries: runtime_attributes.max_retries
  }
}
