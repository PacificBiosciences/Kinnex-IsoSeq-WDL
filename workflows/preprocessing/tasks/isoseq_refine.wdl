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
    String source_dataset_name
    String hifi_barcode
    String cdna_barcode
    String bio_sample
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
    description: "Gather per-partition isoseq refine reports into one table-oriented PacBio report."
    outputs: {
      refine_summary_report: {
        description: "Combined isoseq refine summary report"
      }
    }
  }

  parameter_meta {
    refine_report_rows: {
      description: "Serialized per-partition isoseq refine report rows"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    Array[File] refine_report_rows
    RuntimeAttributes runtime_attributes
  }

  Int threads = 1
  Int mem_gb = 4
  Int total_mem_mb = (mem_gb * 1024) + runtime_attributes.add_memory_mb

  command <<<
    set -euo pipefail

    printf '%s\n' "~{sep="\" \"" refine_report_rows}" > refine_report_rows.txt

    python3 - refine_report_rows.txt isoseq_refine.report.json <<'PY'
    import json
    import re
    import sys
    import uuid

    report_rows_manifest, output_json = sys.argv[1:]

    report_id = 'isoseq_refine'
    report_title = 'Iso-Seq Refine Report'
    table_id = 'isoseq_refine_samples'
    table_title = 'Summary Metrics'
    row_fields = {
        'report_json',
        'bio_sample',
        'source_dataset',
        'hifi_barcode',
        'cdna_barcode',
    }
    identity_columns = (
        ('bio_sample', 'Bio Sample'),
        ('source_dataset', 'Source Dataset'),
        ('hifi_barcode', 'HiFi Barcode'),
        ('cdna_barcode', 'cDNA Barcode'),
    )
    sample_name_metric_id = 'sample_name'
    identifier_re = re.compile(r'^[a-z0-9_]+$')


    def fail(message):
        sys.exit(message)


    def read_paths(path):
        with open(path, 'r', encoding='utf-8') as fh:
            return [line.rstrip('\n') for line in fh if line.rstrip('\n')]


    def require_nonempty_string(value, label):
        if not isinstance(value, str) or not value:
            fail(f'{label} must be a non-empty string')
        return value


    def require_uuid(value, label):
        value = require_nonempty_string(value, label)
        try:
            uuid.UUID(value)
        except ValueError:
            fail(f'{label} is not a valid UUID: {value!r}')
        return value


    def metric_type(value, label):
        if isinstance(value, bool) or value is None:
            fail(f'{label} has unsupported value {value!r}')
        if isinstance(value, int):
            return 'int'
        if isinstance(value, float):
            return 'float'
        if isinstance(value, str):
            return 'string'
        fail(f'{label} has unsupported value type {type(value).__name__}')


    def require_empty_collection(report, field, row_label):
        value = report.get(field, [])
        if not isinstance(value, list):
            fail(f'{row_label}: report {field} must be an array')
        if value:
            fail(f'{row_label}: report {field} must be empty for table-only gathering')


    row_paths = read_paths(report_rows_manifest)
    if not row_paths:
        fail('refine_report_rows must contain at least one report row')

    rows = []
    seen_keys = set()
    expected_version = None
    expected_source_metric_schema = None
    expected_metric_schema = None
    dataset_uuids = set()
    report_tags = set()
    uuid_inputs = []

    for row_index, row_path in enumerate(row_paths, start=1):
        with open(row_path, 'r', encoding='utf-8') as row_fh:
            row = json.load(row_fh)
        if not isinstance(row, dict):
            fail(f'refine report row {row_index} must be a JSON object')
        if set(row) != row_fields:
            fail(f'refine report row {row_index} has fields {sorted(row)!r}; expected {sorted(row_fields)!r}')

        identity = {
            field: require_nonempty_string(row[field], f'refine report row {row_index} {field}')
            for field, _header in identity_columns
        }
        row_key = (
            identity['source_dataset'],
            identity['hifi_barcode'],
            identity['cdna_barcode'],
        )
        row_label = '/'.join(row_key)
        if row_key in seen_keys:
            fail(f'duplicate refine report row key (Source Dataset, HiFi Barcode, cDNA Barcode)={row_key!r}')
        seen_keys.add(row_key)

        report_json = require_nonempty_string(
            row['report_json'],
            f'{row_label}: report_json',
        )
        try:
            report = json.loads(report_json)
        except json.JSONDecodeError as exc:
            fail(f'{row_label}: report_json is invalid JSON: {exc}')
        if not isinstance(report, dict):
            fail(f'{row_label}: report_json must contain a JSON object')
        if report.get('id') != report_id:
            fail(f'{row_label}: report id must be {report_id!r}; observed {report.get("id")!r}')

        version = require_nonempty_string(report.get('version'), f'{row_label}: report version')
        if expected_version is None:
            expected_version = version
        elif version != expected_version:
            fail(f'{row_label}: report version {version!r} does not match first report version {expected_version!r}')

        report_uuid = require_uuid(report.get('uuid'), f'{row_label}: report uuid')
        uuid_inputs.append(
            {
                'row_key': row_key,
                'report_uuid': report_uuid,
            }
        )

        require_empty_collection(report, 'tables', row_label)
        require_empty_collection(report, 'plotGroups', row_label)

        attributes = report.get('attributes')
        if not isinstance(attributes, list) or not attributes:
            fail(f'{row_label}: report attributes must be a non-empty array')

        source_metric_schema = []
        metric_schema = []
        metric_values = []
        seen_metric_ids = set()
        for attribute_index, attribute in enumerate(attributes, start=1):
            attribute_label = f'{row_label}: attribute {attribute_index}'
            if not isinstance(attribute, dict):
                fail(f'{attribute_label} must be an object')
            raw_metric_id = require_nonempty_string(attribute.get('id'), f'{attribute_label} id')
            metric_id = raw_metric_id.rsplit('.', 1)[-1]
            if not identifier_re.fullmatch(metric_id):
                fail(
                    f'{attribute_label} id leaf {metric_id!r} must contain only '
                    'lowercase alphanumeric characters and underscores'
                )
            if metric_id in seen_metric_ids:
                fail(f'{row_label}: duplicate metric id leaf {metric_id!r}')
            if metric_id in dict(identity_columns):
                fail(f'{row_label}: metric id {metric_id!r} conflicts with an identity column')
            seen_metric_ids.add(metric_id)

            metric_name = require_nonempty_string(attribute.get('name'), f'{attribute_label} name')
            if 'value' not in attribute:
                fail(f'{attribute_label} is missing value')
            value = attribute['value']
            value_type = metric_type(value, f'{attribute_label} value')
            source_metric_schema.append((metric_id, metric_name, value_type))
            if metric_id == sample_name_metric_id:
                if value_type != 'string':
                    fail(f'{attribute_label} must be a string')
                if value != identity['bio_sample']:
                    fail(f'{row_label}: refine Sample Name {value!r} does not match Bio Sample {identity["bio_sample"]!r}')
                continue
            metric_schema.append((metric_id, metric_name, value_type))
            metric_values.append(value)

        if expected_source_metric_schema is None:
            expected_source_metric_schema = source_metric_schema
            expected_metric_schema = metric_schema
        elif source_metric_schema != expected_source_metric_schema:
            fail(
                f'{row_label}: metric schema {source_metric_schema!r} does not match '
                f'first report schema {expected_source_metric_schema!r}'
            )

        row_dataset_uuids = report.get('dataset_uuids', [])
        if not isinstance(row_dataset_uuids, list):
            fail(f'{row_label}: report dataset_uuids must be an array')
        for dataset_index, dataset_uuid in enumerate(row_dataset_uuids, start=1):
            dataset_uuids.add(
                require_uuid(
                    dataset_uuid,
                    f'{row_label}: dataset_uuids[{dataset_index}]',
                )
            )

        tags = report.get('tags', [])
        if not isinstance(tags, list) or any(not isinstance(tag, str) for tag in tags):
            fail(f'{row_label}: report tags must be an array of strings')
        report_tags.update(tags)

        rows.append(
            {
                **identity,
                'metric_values': metric_values,
            }
        )

    columns = [
        {
            'id': f'{report_id}.{table_id}.{column_id}',
            'header': header,
            'values': [row[column_id] for row in rows],
        }
        for column_id, header in identity_columns
    ]
    columns.extend(
        {
            'id': f'{report_id}.{table_id}.{metric_id}',
            'header': metric_name,
            'values': [row['metric_values'][metric_index] for row in rows],
        }
        for metric_index, (metric_id, metric_name, _value_type) in enumerate(expected_metric_schema)
    )

    expected_column_length = len(rows)
    for column in columns:
        if len(column['values']) != expected_column_length:
            fail(
                f'internal error: column {column["id"]!r} has {len(column["values"])} '
                f'values for {expected_column_length} rows'
            )

    uuid_seed = json.dumps(
        uuid_inputs,
        sort_keys=True,
        separators=(',', ':'),
    )
    combined_uuid = str(
        uuid.uuid5(
            uuid.NAMESPACE_URL,
            f'kinnex-isoseq-wdl/{report_id}/{uuid_seed}',
        )
    )
    combined_report = {
        '_comment': f'Combined from {len(rows)} per-partition isoseq refine reports',
        'id': report_id,
        'version': expected_version,
        'uuid': combined_uuid,
        'title': report_title,
        'attributes': [],
        'tables': [
            {
                'id': f'{report_id}.{table_id}',
                'title': table_title,
                'columns': columns,
            }
        ],
        'plotGroups': [],
        'dataset_uuids': sorted(dataset_uuids),
        'tags': sorted(report_tags),
    }
    with open(output_json, 'w', encoding='utf-8') as out_fh:
        json.dump(combined_report, out_fh, indent=2, sort_keys=True)
        out_fh.write('\n')
    PY
  >>>

  output {
    File refine_summary_report = "isoseq_refine.report.json"
  }

  runtime {
    cpu: threads
    memory: total_mem_mb + " MB"
    docker: runtime_attributes.container_registry + "/pb_wdl_base@sha256:03cb3c01937eccc907f8ad71c87b258581504572205fe3f31a657e318f3564ae"
    maxRetries: runtime_attributes.max_retries
  }
}
