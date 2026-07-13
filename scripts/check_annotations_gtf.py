#!/usr/bin/env python3
"""Check or canonicalize GTF files for Pigeon classification."""

from __future__ import annotations

import argparse
import gzip
from collections import OrderedDict
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path
from typing import TextIO

GENE_FEATURES = {'gene', 'protein_coding_gene'}
TRANSCRIPT_FEATURES = {'transcript', 'mRNA'}
SKIP_PREFIXES = ('#', 'browser', 'track')


@dataclass
class GtfRecord:
    line: str
    fields: list[str]
    attributes: dict[str, str]
    line_number: int
    order: int

    @property
    def ref_name(self) -> str:
        return self.fields[0]

    @property
    def feature(self) -> str:
        return self.fields[2]

    @property
    def start(self) -> int:
        return int(self.fields[3])

    @property
    def end(self) -> int:
        return int(self.fields[4])

    @property
    def gene_id(self) -> str | None:
        return self.attributes.get('gene_id')

    @property
    def transcript_id(self) -> str | None:
        return self.attributes.get('transcript_id')


@dataclass
class TranscriptBucket:
    transcript_records: list[GtfRecord] = field(default_factory=list)
    child_records: list[GtfRecord] = field(default_factory=list)
    order: int = 0


@dataclass
class GeneBucket:
    gene_records: list[GtfRecord] = field(default_factory=list)
    transcripts: OrderedDict[str, TranscriptBucket] = field(default_factory=OrderedDict)
    orphan_records: list[GtfRecord] = field(default_factory=list)


@dataclass
class CanonicalizeStats:
    input_records: int = 0
    output_lines: int = 0
    skipped_orphans: int = 0
    skipped_transcript_children: int = 0


def open_text_for_read(path: Path) -> TextIO:
    if path.suffix == '.gz':
        return gzip.open(path, 'rt', encoding='utf-8')
    return path.open('rt', encoding='utf-8')


def open_text_for_write(path: Path) -> TextIO:
    if path.suffix == '.gz':
        return gzip.open(path, 'wt', encoding='utf-8')
    return path.open('wt', encoding='utf-8')


def parse_attributes(attributes: str) -> dict[str, str]:
    parsed = {}
    for raw_item in attributes.split(';'):
        item = raw_item.strip()
        if not item:
            continue
        if '=' in item and (' ' not in item or item.index('=') < item.index(' ')):
            key, value = item.split('=', 1)
        else:
            parts = item.split(None, 1)
            if len(parts) != 2:
                continue
            key, value = parts
        parsed[key] = value.strip().strip('"')
    return parsed


def parse_gtf_record(line: str, line_number: int, order: int) -> GtfRecord:
    fields = line.split('\t')
    if len(fields) != 9:
        raise ValueError(f'line {line_number}: expected 9 columns, found {len(fields)}')

    try:
        int(fields[3])
        int(fields[4])
    except ValueError as exc:
        raise ValueError(f'line {line_number}: start and end must be integer coordinates') from exc

    return GtfRecord(
        line=line,
        fields=fields,
        attributes=parse_attributes(fields[8]),
        line_number=line_number,
        order=order,
    )


def iter_gtf_records(path: Path) -> Iterator[GtfRecord]:
    order = 0
    with open_text_for_read(path) as source:
        for line_number, line in enumerate(source, 1):
            line = line.rstrip('\n')
            if not line or line.startswith(SKIP_PREFIXES):
                continue
            yield parse_gtf_record(line, line_number, order)
            order += 1


def check_pigeon_gtf(path: Path, max_errors: int) -> Iterator[str]:
    current_ref = None
    current_gene = None
    current_gene_record: GtfRecord | None = None
    current_transcript = None
    current_transcript_record: GtfRecord | None = None
    seen_refs = set()
    first_record = True
    errors = 0

    try:
        records = iter_gtf_records(path)
        for record in records:
            if first_record:
                first_record = False
                if record.feature not in GENE_FEATURES | TRANSCRIPT_FEATURES:
                    errors += 1
                    yield (
                        f'line {record.line_number}: first data record must be a gene '
                        f'or transcript record, found {record.feature}'
                    )
                    if errors >= max_errors:
                        return

            if record.ref_name != current_ref:
                if record.ref_name in seen_refs:
                    errors += 1
                    yield (
                        f'line {record.line_number}: reference {record.ref_name} '
                        'appears in multiple non-contiguous blocks'
                    )
                    if errors >= max_errors:
                        return
                seen_refs.add(record.ref_name)
                current_ref = record.ref_name
                current_gene = None
                current_gene_record = None
                current_transcript = None
                current_transcript_record = None

            if record.feature in GENE_FEATURES:
                current_ref = record.ref_name
                current_gene = record.gene_id
                current_gene_record = record
                current_transcript = None
                current_transcript_record = None
                continue

            if record.feature in TRANSCRIPT_FEATURES:
                current_ref = record.ref_name
                current_gene = record.gene_id
                current_transcript = record.transcript_id
                current_transcript_record = record

                if not current_transcript:
                    errors += 1
                    yield f'line {record.line_number}: transcript record is missing transcript_id'
                if current_gene_record and record.ref_name == current_gene_record.ref_name:
                    if record.gene_id != current_gene_record.gene_id:
                        errors += 1
                        yield (
                            f'line {record.line_number}: transcript gene_id '
                            f'{record.gene_id or "<missing gene_id>"} does not match '
                            f'current gene block {current_gene_record.gene_id or "<missing gene_id>"}'
                        )
                    if record.start < current_gene_record.start:
                        errors += 1
                        yield (
                            f'line {record.line_number}: transcript starts before current '
                            f'gene block ({record.start} < {current_gene_record.start})'
                        )
                    if record.end > current_gene_record.end:
                        errors += 1
                        yield (
                            f'line {record.line_number}: transcript ends after current '
                            f'gene block ({record.end} > {current_gene_record.end})'
                        )
                if errors >= max_errors:
                    return
                continue

            if not record.transcript_id:
                continue

            if (
                record.ref_name != current_ref
                or record.gene_id != current_gene
                or record.transcript_id != current_transcript
            ):
                errors += 1
                yield (
                    f'line {record.line_number}: {record.feature} belongs to '
                    f'{record.transcript_id or "<missing transcript_id>"}, but current '
                    f'transcript block is {current_transcript or "<none>"}'
                )
                if errors >= max_errors:
                    return
                continue

            if not current_transcript_record:
                errors += 1
                yield (f'line {record.line_number}: {record.feature} appears before a parent transcript record')
                if errors >= max_errors:
                    return
                continue

            if record.start < current_transcript_record.start:
                errors += 1
                yield (
                    f'line {record.line_number}: {record.feature} starts before current '
                    f'transcript block ({record.start} < {current_transcript_record.start})'
                )
            if record.end > current_transcript_record.end:
                errors += 1
                yield (
                    f'line {record.line_number}: {record.feature} ends after current '
                    f'transcript block ({record.end} > {current_transcript_record.end})'
                )
            if errors >= max_errors:
                return
        if first_record:
            yield 'no GTF data records found'
    except ValueError as exc:
        yield str(exc)


def read_headers(path: Path) -> list[str]:
    headers = []
    with open_text_for_read(path) as source:
        for line in source:
            line = line.rstrip('\n')
            if line.startswith('##'):
                headers.append(line)
                continue
            if not line or line.startswith('#'):
                continue
            break
    return headers


def record_sort_key(record: GtfRecord) -> tuple[int, int, int]:
    return (record.start, record.end, record.order)


def transcript_sort_key(item: tuple[str, TranscriptBucket]) -> tuple[int, int, int]:
    _, transcript = item
    if transcript.transcript_records:
        return record_sort_key(transcript.transcript_records[0])
    return (10**18, 10**18, transcript.order)


def gene_bucket(
    records_by_ref: OrderedDict[str, OrderedDict[str, GeneBucket]],
    ref_name: str,
    gene_id: str,
) -> GeneBucket:
    ref_bucket = records_by_ref.setdefault(ref_name, OrderedDict())
    if gene_id not in ref_bucket:
        ref_bucket[gene_id] = GeneBucket()
    return ref_bucket[gene_id]


def transcript_bucket(gene: GeneBucket, transcript_id: str) -> TranscriptBucket:
    if transcript_id not in gene.transcripts:
        gene.transcripts[transcript_id] = TranscriptBucket(order=len(gene.transcripts))
    return gene.transcripts[transcript_id]


def canonicalize_pigeon_gtf(input_path: Path, output_path: Path) -> CanonicalizeStats:
    records_by_ref: OrderedDict[str, OrderedDict[str, GeneBucket]] = OrderedDict()
    stats = CanonicalizeStats()

    for record in iter_gtf_records(input_path):
        stats.input_records += 1
        gene_id = record.gene_id or f'__missing_gene_id_{record.ref_name}_{record.order}'
        gene = gene_bucket(records_by_ref, record.ref_name, gene_id)

        if record.feature in GENE_FEATURES:
            gene.gene_records.append(record)
        elif record.feature in TRANSCRIPT_FEATURES and record.transcript_id:
            transcript_bucket(gene, record.transcript_id).transcript_records.append(record)
        elif record.transcript_id:
            transcript_bucket(gene, record.transcript_id).child_records.append(record)
        else:
            gene.orphan_records.append(record)

    headers = read_headers(input_path)
    with open_text_for_write(output_path) as output:
        for header in headers:
            output.write(header + '\n')
            stats.output_lines += 1

        for genes in records_by_ref.values():
            for gene in genes.values():
                gene_records = sorted(gene.gene_records, key=record_sort_key)
                if gene_records:
                    output.write(gene_records[0].line + '\n')
                    stats.output_lines += 1
                    stats.skipped_orphans += len(gene.orphan_records)
                else:
                    stats.skipped_orphans += len(gene.orphan_records)

                for _, transcript in sorted(gene.transcripts.items(), key=transcript_sort_key):
                    transcript_records = sorted(transcript.transcript_records, key=record_sort_key)
                    if not transcript_records:
                        stats.skipped_transcript_children += len(transcript.child_records)
                        continue

                    output.write(transcript_records[0].line + '\n')
                    stats.output_lines += 1
                    for child in sorted(transcript.child_records, key=record_sort_key):
                        output.write(child.line + '\n')
                        stats.output_lines += 1

    return stats


def require_existing_file(path: Path) -> int | None:
    if path.exists():
        return None
    print(f'FAILED: file not found: {path}')
    return 2


def command_check(args: argparse.Namespace) -> int:
    if error_code := require_existing_file(args.gtf):
        return error_code

    errors = list(check_pigeon_gtf(args.gtf, args.max_errors))
    if errors:
        for error in errors:
            print(error)
        print(f'FAILED: found at least {len(errors)} block-order problem(s)')
        return 1

    print('OK: GTF appears block-ordered for Pigeon')
    return 0


def command_canonicalize(args: argparse.Namespace) -> int:
    if error_code := require_existing_file(args.input_gtf):
        return error_code
    if args.output_gtf.exists() and not args.force:
        print(f'FAILED: output exists: {args.output_gtf} (use --force to overwrite)')
        return 2

    stats = canonicalize_pigeon_gtf(args.input_gtf, args.output_gtf)
    errors = list(check_pigeon_gtf(args.output_gtf, args.max_errors))
    if errors:
        for error in errors:
            print(error)
        print(f'FAILED: canonicalized output is still not block-ordered: {args.output_gtf}')
        return 1

    print(f'Canonicalized {stats.input_records} annotation records into {stats.output_lines} output lines')
    print(
        'Skipped '
        f'{stats.skipped_orphans} gene-level orphan records and '
        f'{stats.skipped_transcript_children} child records without transcript records'
    )
    print(f'OK: wrote Pigeon block-ordered GTF to {args.output_gtf}')
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description='Check or canonicalize GTF files for Pigeon classification.')
    subparsers = parser.add_subparsers(dest='command', required=True)

    check_parser = subparsers.add_parser(
        'check',
        help='Check whether a .gtf or .gtf.gz is block-ordered for Pigeon',
    )
    check_parser.add_argument('gtf', type=Path, help='Input .gtf or .gtf.gz file')
    check_parser.add_argument(
        '--max-errors',
        type=int,
        default=20,
        help='Maximum number of diagnostics to print before stopping',
    )
    check_parser.set_defaults(func=command_check)

    canonicalize_parser = subparsers.add_parser(
        'canonicalize',
        aliases=['fix'],
        help='Rewrite a .gtf or .gtf.gz into Pigeon block order',
    )
    canonicalize_parser.add_argument(
        'input_gtf',
        type=Path,
        help='Input .gtf or .gtf.gz file',
    )
    canonicalize_parser.add_argument(
        'output_gtf',
        type=Path,
        help='Output .gtf or .gtf.gz file',
    )
    canonicalize_parser.add_argument(
        '--force',
        action='store_true',
        help='Overwrite output_gtf if it already exists',
    )
    canonicalize_parser.add_argument(
        '--max-errors',
        type=int,
        default=20,
        help='Maximum number of validation diagnostics to print after writing',
    )
    canonicalize_parser.set_defaults(func=command_canonicalize)

    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    return args.func(args)


if __name__ == '__main__':
    raise SystemExit(main())
