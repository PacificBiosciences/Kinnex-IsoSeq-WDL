version 1.0

import "../../rna_structs.wdl"

task group_flnc_bams_by_sm {
  meta {
    description: "Group FLNC BAMs by the single original SM tag present in each BAM header and derive sanitized filename prefixes."
    outputs: {
      sample_names: {
        description: "Sample names"
      },
      sample_prefixes: {
        description: "Sample prefixes"
      },
      group_sizes: {
        description: "FLNC BAM counts per sample group"
      },
      group_bam_indices: {
        description: "Input FLNC BAM indices per sample group"
      },
      flnc_grouping_report: {
        description: "FLNC grouping report TSV"
      }
    }
  }

  parameter_meta {
    flnc_bams: {
      description: "FLNC BAMs"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    Array[File] flnc_bams
    RuntimeAttributes runtime_attributes
  }

  Int threads = 1
  Int mem_gb = 4
  Int total_mem_mb = (mem_gb * 1024) + runtime_attributes.add_memory_mb

  command <<<
    set -euo pipefail

    printf '%s\n' "~{sep="\" \"" flnc_bams}" > flnc_bams.txt

    python3 - flnc_bams.txt <<'PY'
    import collections
    import json
    import re
    import sys

    import pysam

    bam_list_path = sys.argv[1]
    max_sample_prefix_length = 40
    pysam.set_verbosity(0)


    def fail(message):
        sys.exit(message)


    def sanitize_sample_prefix(sample_name):
        sample_prefix = re.sub(r'[^A-Za-z0-9]+', '_', sample_name).strip('_')
        sample_prefix = sample_prefix[:max_sample_prefix_length].strip('_')
        if not sample_prefix:
            fail(
                f'@RG SM value {sample_name!r} produces an empty sanitized sample prefix; '
                'sample names must contain at least one alphanumeric character'
            )
        return sample_prefix


    def read_single_sm_value(bam):
        try:
            with pysam.AlignmentFile(bam, 'rb', check_sq=False) as bam_fh:
                read_groups = bam_fh.header.to_dict().get('RG', [])
        except (OSError, ValueError) as exc:
            fail(f'failed to read BAM header for {bam}: {exc}')

        if not read_groups:
            fail(f'{bam} has no @RG records')

        observed_samples = set()
        for index, read_group in enumerate(read_groups, start=1):
            read_group_id = read_group.get('ID', f'record {index}')
            if 'SM' not in read_group or read_group['SM'] == '':
                fail(f'{bam} @RG {read_group_id!r} is missing SM tag')
            sample_name = read_group['SM']
            if not sample_name.strip():
                fail(f'{bam} @RG {read_group_id!r} has an empty or whitespace-only SM value')
            if any(char in sample_name for char in '\t\r\n'):
                fail(f'{bam} has an SM value with tab, carriage return, or newline characters: {sample_name!r}')
            observed_samples.add(sample_name)

        if len(observed_samples) != 1:
            fail(f'{bam} has multiple @RG SM values: ' + ', '.join(repr(sample) for sample in sorted(observed_samples)))

        sample_name = next(iter(observed_samples))
        return sample_name


    with open(bam_list_path, 'r', encoding='utf-8') as fh:
        flnc_bams = [line.rstrip('\n') for line in fh if line.rstrip('\n')]

    if not flnc_bams:
        fail('flnc_bams must contain at least one FLNC BAM')

    counts = collections.Counter(flnc_bams)
    duplicate_bams = sorted(bam for bam, count in counts.items() if count > 1)
    if duplicate_bams:
        fail('duplicate FLNC BAM inputs are not allowed: ' + ', '.join(duplicate_bams))

    sample_to_bams = collections.OrderedDict()
    sample_to_bam_indices = {}
    sample_to_prefix = {}
    prefix_to_sample = {}

    for bam_index, bam in enumerate(flnc_bams):
        sample_name = read_single_sm_value(bam)

        if sample_name not in sample_to_bams:
            sample_prefix = sanitize_sample_prefix(sample_name)
            if sample_prefix in prefix_to_sample:
                fail(
                    'different @RG SM values sanitize to the same sample prefix '
                    f'{sample_prefix!r}: {prefix_to_sample[sample_prefix]!r}, {sample_name!r}'
                )
            prefix_to_sample[sample_prefix] = sample_name
            sample_to_bams[sample_name] = []
            sample_to_bam_indices[sample_name] = []
            sample_to_prefix[sample_name] = sample_prefix
        sample_to_bams[sample_name].append(bam)
        sample_to_bam_indices[sample_name].append(bam_index)

    group_sizes = []
    group_bam_indices = []
    with (
        open('sample_names.txt', 'w', encoding='utf-8') as sample_fh,
        open('sample_prefixes.txt', 'w', encoding='utf-8') as prefix_fh,
        open('flnc_grouping_report.tsv', 'w', encoding='utf-8') as report_fh,
    ):
        report_fh.write('sample_name\tsample_prefix\tgroup_size\taction\tinput_bam\n')
        for sample_name, bams in sample_to_bams.items():
            bam_indices = sample_to_bam_indices[sample_name]
            sample_prefix = sample_to_prefix[sample_name]
            action = 'merge' if len(bams) > 1 else 'passthrough'
            sample_fh.write(sample_name + '\n')
            prefix_fh.write(sample_prefix + '\n')
            group_sizes.append(len(bams))
            group_bam_indices.append(bam_indices)
            for bam in bams:
                report_fh.write(f'{sample_name}\t{sample_prefix}\t{len(bams)}\t{action}\t{bam}\n')

    with (
        open('group_sizes.json', 'w', encoding='utf-8') as size_fh,
        open('group_bam_indices.json', 'w', encoding='utf-8') as group_index_fh,
    ):
        json.dump(group_sizes, size_fh)
        json.dump(group_bam_indices, group_index_fh)
    PY
  >>>

  output {
    Array[String] sample_names = read_lines("sample_names.txt")
    Array[String] sample_prefixes = read_lines("sample_prefixes.txt")
    Array[Int] group_sizes = read_json("group_sizes.json")
    Array[Array[Int]] group_bam_indices = read_json("group_bam_indices.json")
    File flnc_grouping_report = "flnc_grouping_report.tsv"
  }

  runtime {
    cpu: threads
    memory: total_mem_mb + " MB"
    docker: runtime_attributes.container_registry + "/pb_wdl_base@sha256:03cb3c01937eccc907f8ad71c87b258581504572205fe3f31a657e318f3564ae"
    maxRetries: runtime_attributes.max_retries
  }
}

task samtools_merge {
  meta {
    description: "Merge one FLNC sample group into a single grouped BAM with samtools."
    outputs: {
      grouped_flnc_bam: {
        description: "Grouped FLNC BAM"
      }
    }
  }

  parameter_meta {
    sample_name: {
      description: "Sample name"
    }
    sample_prefix: {
      description: "Sample prefix"
    }
    flnc_bams: {
      description: "FLNC BAMs"
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
    String sample_name
    String sample_prefix
    Array[File] flnc_bams
    Int threads = 8
    Int mem_gb = 32
    RuntimeAttributes runtime_attributes
  }

  Int total_mem_mb = (mem_gb * 1024) + runtime_attributes.add_memory_mb
  Int effective_threads = if (threads > runtime_attributes.nproc)
    then runtime_attributes.nproc
    else threads

  File sample_name_file = write_lines([
    sample_name
  ])
  String output_prefix = sample_prefix

  command <<<
    set -euo pipefail

    sample_name="$(cat "~{sample_name_file}")"
    printf '%s\n' "~{sep="\" \"" flnc_bams}" > selected_flnc_bams.txt

    group_size="$(grep -cve '^[[:space:]]*$' selected_flnc_bams.txt || true)"
    output_bam="~{output_prefix}.flnc.bam"

    if [ "${group_size}" -lt 1 ]; then
      echo "samtools_merge requires at least one FLNC BAM for sample ${sample_name}; observed ${group_size}" >&2
      exit 1
    fi

    samtools_threads_args=()
    if [ "~{effective_threads}" -gt 1 ]; then
      samtools_threads_args=(--threads "$((~{effective_threads} - 1))")
    fi

    samtools merge \
      "${samtools_threads_args[@]}" \
      -c \
      -p \
      -o "${output_bam}" \
      -b selected_flnc_bams.txt

    action="merge"

    {
      printf 'sample_name\tsample_prefix\tgroup_size\taction\toutput_bam\n'
      printf '%s\t%s\t%s\t%s\t%s\n' \
        "${sample_name}" \
        "~{sample_prefix}" \
        "${group_size}" \
        "${action}" \
        "${output_bam}"
    } > "~{output_prefix}.samtools_merge.log"
  >>>

  output {
    File grouped_flnc_bam = "~{output_prefix}.flnc.bam"
  }

  runtime {
    cpu: effective_threads
    memory: total_mem_mb + " MB"
    docker: runtime_attributes.container_registry + "/pb_wdl_base@sha256:03cb3c01937eccc907f8ad71c87b258581504572205fe3f31a657e318f3564ae"
    maxRetries: runtime_attributes.max_retries
  }
}

task pbmm2_align_flnc {
  meta {
    description: "Align one grouped FLNC sample BAM with pbmm2 using the ISOSEQ preset."
    outputs: {
      aligned_bam: {
        description: "Aligned FLNC BAM"
      },
      aligned_bam_index: {
        description: "Aligned FLNC BAM index"
      }
    }
  }

  parameter_meta {
    sample_prefix: {
      description: "Sample prefix"
    }
    flnc_bam: {
      description: "FLNC BAM"
    }
    genome_fasta: {
      description: "Reference FASTA"
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
    String sample_prefix
    File flnc_bam
    File genome_fasta
    Int threads = 16
    Int mem_gb = 48
    RuntimeAttributes runtime_attributes
  }

  Int total_mem_mb = (mem_gb * 1024) + runtime_attributes.add_memory_mb
  Int effective_threads = if (threads > runtime_attributes.nproc)
    then runtime_attributes.nproc
    else threads

  String output_prefix = sample_prefix

  command <<<
    set -euo pipefail

    pbmm2 align \
      --preset ISOSEQ \
      --sort \
      --bam-index BAI \
      --unmapped \
      --num-threads ~{effective_threads} \
      --log-level INFO \
      --log-file "~{output_prefix}.pbmm2.log" \
      "~{genome_fasta}" \
      "~{flnc_bam}" \
      "~{output_prefix}.aligned.bam"
  >>>

  output {
    File aligned_bam = "~{output_prefix}.aligned.bam"
    File aligned_bam_index = "~{output_prefix}.aligned.bam.bai"
  }

  runtime {
    cpu: effective_threads
    memory: total_mem_mb + " MB"
    docker: runtime_attributes.container_registry + "/pbmm2@sha256:0c21f29f1ee429dbafe5c332d4abffcfed5efbf5448b4ef06dd189b5323a7051"  # 26.2.0_build1
    maxRetries: runtime_attributes.max_retries
  }
}
