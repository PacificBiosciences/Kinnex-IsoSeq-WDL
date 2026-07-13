# Resource bundle layout

The Kinnex Iso-Seq resource bundle contains reusable static files for workflow
runs. It includes both reference resources and preprocessing resources.
These files are consumed through different workflow input paths: reference
resources are selected through the reference map, while preprocessing
resources are passed as explicit workflow inputs.

## Downloading the bundle

Download and extract the Zenodo archive in the parent directory that you will
use as `<local_path_prefix>` in the workflow input templates:

```bash
mkdir -p /path/to/resources
cd /path/to/resources
curl -L -o kinnex-isoseq-wdl-resources-v0.1.0.tar.gz \
  "https://zenodo.org/records/10839617/files/kinnex-isoseq-wdl-resources-v0.1.0.tar.gz?download=1"
echo "15c3fca8bfdcb78d68342277e69d2e1e kinnex-isoseq-wdl-resources-v0.1.0.tar.gz" | md5sum --check -
tar -xzf kinnex-isoseq-wdl-resources-v0.1.0.tar.gz
```

```text
kinnex-isoseq-wdl-resources-v0.1.0/
├── GRCh38/
│   ├── annotation/
│   │   └── gencode.v49.annotation.gtf.gz
│   ├── classification/
│   │   ├── intropolis.v1.hg19_with_liftover_to_hg38.tsv.min_count_10.modified2.sorted.tsv
│   │   ├── polyA.list.txt
│   │   └── refTSS_v3.3_human_coordinate.hg38.sorted.bed
│   ├── human_GRCh38_no_alt_analysis_set.fasta
│   └── human_GRCh38_no_alt_analysis_set.fasta.fai
├── kinnex/
│   ├── isoseq_v2_barcoded_primers/
│   │   └── IsoSeq_v2_primers_12.fasta
│   ├── kinnex_hifi_barcodes/
│   │   └── kinnex_hifi_barcodes.fasta
│   └── kinnex_primers/
│       ├── kinnex_12fold_primers.fasta
│       ├── kinnex_16fold_primers.fasta
│       └── kinnex_8fold_primers.fasta
└── GRCh38.ref_map.v0p1p0.template.tsv
```

## How the bundle maps to workflow inputs

| Bundle content | Workflow input | Notes |
| --- | --- | --- |
| Genome FASTA and FASTA index | `ref_map_file` keys | Used by FLNC alignment, isocall, and isoform classification. |
| Selected annotation GTF | `ref_map_file` key `annotation_gtf_gz` | The bundle can contain multiple annotations, but one reference map selects the active annotation for a run. |
| Classification polyA, CAGE, and junction resources | `ref_map_file` keys | Used by isoform classification when the matching `pigeon_use_*` inputs are enabled. |
| `kinnex/kinnex_hifi_barcodes/kinnex_hifi_barcodes.fasta` | `hifi_demux_barcodes` | Used by preprocessing for upstream HiFi demultiplexing. |
| `kinnex/kinnex_primers/kinnex_8fold_primers.fasta` | `skera_adapters` | Used by preprocessing for `skera split`. Choose another bundled Kinnex fold primer file when appropriate for the assay. |
| `kinnex/isoseq_v2_barcoded_primers/IsoSeq_v2_primers_12.fasta` | `barcoded_primers` | Used by preprocessing for cDNA `lima` and `isoseq refine`. |

## Preparing the reference map

Start from `GRCh38.ref_map.v0p1p0.template.tsv` and make a run-local copy, for
example `GRCh38.ref_map.v0p1p0.hpc.tsv`. Replace the path placeholder with a
filesystem path visible to the workflow engine and backend jobs. In the bundled
template, `<prefix>` is the parent directory that contains
`kinnex-isoseq-wdl-resources-v0.1.0/`.

The public HPC templates use `<local_path_prefix>` for the parent directory that
contains `kinnex-isoseq-wdl-resources-v0.1.0/`. With that convention, the
reference map path is:

```text
<local_path_prefix>/kinnex-isoseq-wdl-resources-v0.1.0/GRCh38.ref_map.v0p1p0.hpc.tsv
```

## Typical preprocessing resource inputs

For `kinnex_isoseq` and `preprocessing` runs that use the bundled Kinnex
resources, set:

```json
{
  "kinnex_isoseq.hifi_demux_barcodes": "<local_path_prefix>/kinnex-isoseq-wdl-resources-v0.1.0/kinnex/kinnex_hifi_barcodes/kinnex_hifi_barcodes.fasta",
  "kinnex_isoseq.skera_adapters": "<local_path_prefix>/kinnex-isoseq-wdl-resources-v0.1.0/kinnex/kinnex_primers/kinnex_8fold_primers.fasta",
  "kinnex_isoseq.barcoded_primers": "<local_path_prefix>/kinnex-isoseq-wdl-resources-v0.1.0/kinnex/isoseq_v2_barcoded_primers/IsoSeq_v2_primers_12.fasta",
  "kinnex_isoseq.ref_map_file": "<local_path_prefix>/kinnex-isoseq-wdl-resources-v0.1.0/GRCh38.ref_map.v0p1p0.hpc.tsv"
}
```

For the standalone `preprocessing` entrypoint, use the same file paths with the
`preprocessing.` input namespace.
