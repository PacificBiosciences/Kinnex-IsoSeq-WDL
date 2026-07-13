# Reference map file specification

The `ref_map_file` input is a two-column TSV map with no header. The first
column is the key, and the second column is the value used by the workflow.

```text
name\tGRCh38
genome_fasta\t/path/to/human_GRCh38_no_alt_analysis_set.fasta
genome_fasta_index\t/path/to/human_GRCh38_no_alt_analysis_set.fasta.fai
```

File paths in the map must be visible to the backend jobs. For public HPC runs,
start from the [HPC reference-map template](../backends/hpc/GRCh38.ref_map.v0p1p0.hpc.tsv)
and replace `<local_path_prefix>` with the parent path of the extracted
`kinnex-isoseq-wdl-resources-v0.1.0` bundle on your filesystem.

The resource bundle containing the GRCh38 reference and other files used in
this workflow can be downloaded from Zenodo: [Zenodo record](https://zenodo.org/records/10839617).
See the [resource bundle layout](./resource_bundle.md) for the package directory
structure and for which bundled files are passed through workflow inputs
instead of `ref_map_file`.

## Stage requirements

| Workflow | Required keys |
| -------- | ------------- |
| `isocall` | `name`, `annotation_gtf_gz`, `genome_fasta`, `genome_fasta_index` |
| `isoform_classification` | `name`, `annotation_gtf_gz`, `genome_fasta`, `genome_fasta_index`, `pigeon_poly_a`, `pigeon_cage_peak_bed`, `pigeon_junction_coverage` |
| `secondary_analysis` | union of FLNC alignment task, `isocall`, and `isoform_classification` keys |
| `kinnex_isoseq` | same `ref_map_file` keys as `secondary_analysis`; preprocessing assay resources are separate workflow inputs |

The `pigeon_use_polya`, `pigeon_use_cage_peak`, and `pigeon_use_junction` inputs
control whether the corresponding `pigeon classify` flags are passed. The
classification support resources are still expected to be present in
`ref_map_file`.

## Keys

| Type | Key | Description | Used by |
| ---- | --- | ----------- | ------- |
| String | `name` | Short reference name emitted as `reference_name` in workflow outputs. | All workflows that consume `ref_map_file` |
| File | `genome_fasta` | Reference genome FASTA. | FLNC alignment task, isocall, isoform classification |
| File | `genome_fasta_index` | Reference FASTA index. | Isocall, isoform classification |
| File | `annotation_gtf_gz` | Gzip-compressed genome annotation GTF used by `isocall prep-isoforms` and decompressed internally for `pigeon classify`. The decompressed records must already be in Pigeon input order. | Isocall, isoform classification |
| File | `pigeon_poly_a` | PolyA motif list passed to `pigeon classify --poly-a` when enabled. | Isoform classification |
| File | `pigeon_cage_peak_bed` | CAGE peak BED prepared internally and passed to `pigeon classify --cage-peak` when enabled. | Isoform classification |
| File | `pigeon_junction_coverage` | Junction coverage TSV prepared internally and passed to `pigeon classify --coverage` when enabled. | Isoform classification |

## Annotation GTF structure

The `annotation_gtf_gz` file must be organized as nested gene and transcript
blocks. Records must be ordered like this:

```text
gene
transcript for gene
child records for that transcript, such as exon and CDS
transcript for gene
child records for that transcript, such as exon and CDS
gene
...
```

The workflow treats `gene` and `protein_coding_gene` records as gene starts,
and `transcript` and `mRNA` records as transcript starts. Child records with a
`transcript_id` must appear after their parent transcript record and before the
next transcript or gene record. They must also carry the same `transcript_id` as
the current transcript block and fall within the transcript coordinates. Records
for each reference sequence must be contiguous. The workflow decompresses
`annotation_gtf_gz` and validates that the resulting GTF uses this block
structure before indexing and classification.

For a quick local check, run:

```bash
scripts/check_annotations_gtf.py check annotation.gtf
scripts/check_annotations_gtf.py check annotation.gtf.gz
```

To rewrite a coordinate-sorted GTF into the expected block order:

```bash
scripts/check_annotations_gtf.py canonicalize annotation.gtf.gz annotation.pigeon.gtf.gz
```
