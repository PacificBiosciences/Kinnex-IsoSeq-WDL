# Biosample CSV sample sheet specification

The `biosample_csv` input is the preprocessing sample sheet used by
`preprocessing.wdl` and the end-to-end `kinnex_isoseq.wdl`. It maps each
outer SMRTbell barcode and cDNA Barcode combination to the biological sample
name that should be written to the final FLNC BAM `SM` tag.

The biosample CSV shall be in CSV format with column headers `HiFi Barcode`,
`cDNA Barcode`, and `Bio Sample`, written as the exact three-column header:

```csv
HiFi Barcode,cDNA Barcode,Bio Sample
bcM0001--bcM0001,IsoSeqX_bc01_5p--IsoSeqX_3p,sample_1
bcM0001--bcM0001,IsoSeqX_bc02_5p--IsoSeqX_3p,sample_2
bcM0002--bcM0002,IsoSeqX_bc01_5p--IsoSeqX_3p,sample_3
```

## Columns

| Column | Description |
| ------ | ----------- |
| `HiFi Barcode` | Outer SMRTbell barcode pair for the source dataset, written as an exact lima barcode-pair name such as `bcM0001--bcM0001`. |
| `cDNA Barcode` | cDNA Barcode pair used by downstream cDNA `lima`, written as an exact barcode-pair name such as `IsoSeqX_bc01_5p--IsoSeqX_3p`. |
| `Bio Sample` | Biological sample name. This value becomes the `SM` tag in each final FLNC BAM. |

HiFi Barcode names are exact lima barcode-pair names. The workflow does not
infer or normalize names such as `bcM0001` into `bcM0001--bcM0001`.

## Format expectations

- The header must be exactly `HiFi Barcode,cDNA Barcode,Bio Sample`.
- Empty rows are ignored.
- Each non-empty row must have exactly three non-empty fields.
- Leading and trailing whitespace is stripped from each field.
- `HiFi Barcode` and `cDNA Barcode` values must each be exact two-part barcode pairs
  separated by `--`.
- `Bio Sample` names must be 40 characters or fewer and may contain only
  alphanumeric characters (`A-Z`, `a-z`, `0-9`), underscores (`_`), and hyphens
  (`-`).
- Within a single `HiFi Barcode` group, each `cDNA Barcode` must be unique.
- Different `HiFi Barcode` groups may reuse the same `cDNA Barcode` and/or `Bio Sample`.

Choose stable `Bio Sample` names that are safe for downstream grouping and output
prefixes. FLNC alignment groups BAMs by the exact `SM` value. Preprocessing keeps
this strict `Bio Sample` naming contract, while the public `secondary_analysis`
entrypoint is more permissive for external FLNC BAM `SM` values and sanitizes
those values into filename prefixes.

## Run modes

In HiFi demux mode, the workflow runs upstream HiFi demux first. Each component
of every `HiFi Barcode` pair must appear in `hifi_demux_barcodes`, and the full
pair must match the outer-barcode name produced by `lima --split-named`.
Because the workflow uses the symmetric-adapter HiFi demux preset, each
`HiFi Barcode` pair in this mode must use the same barcode on both sides, such
as `bcM0001--bcM0001`. `hifi_sources` must contain exactly one raw HiFi BAM and
must omit `hifi_barcode`.

In cDNA demux-only mode, each input BAM is treated as already demuxed on
instrument. Every `hifi_sources` entry must contain `hifi_barcode`. Each value
must match a `HiFi Barcode` value in `biosample_csv` and each component must
appear in `hifi_demux_barcodes`. Multiple source BAMs may use the same barcode
pair, for example when the same barcode set was reused across movies or runs.

For both modes, each component of every `cDNA Barcode` pair is checked against
the FASTA record names in `barcoded_primers`.

## Derived cDNA biosample CSVs

For each normalized outer-barcode dataset, the workflow filters
`biosample_csv` to that `HiFi Barcode` group and derives a lima-compatible
two-column CSV:

```csv
Barcodes,Bio Sample
IsoSeqX_bc01_5p--IsoSeqX_3p,sample_1
IsoSeqX_bc02_5p--IsoSeqX_3p,sample_2
```

This derived file is passed to downstream cDNA `lima`.

## Examples

### Splitting one raw HiFi BAM

For a HiFi BAM in HiFi demux mode, the three-column sheet can contain multiple
outer barcode groups:

```csv
HiFi Barcode,cDNA Barcode,Bio Sample
bcM0001--bcM0001,bcU0001_5p--bcU0001_3p,sample_1
bcM0001--bcM0001,bcU0002_5p--bcU0002_3p,sample_2
bcM0001--bcM0001,bcU0003_5p--bcU0003_3p,sample_3
bcM0002--bcM0002,bcU0001_5p--bcU0001_3p,sample_4
bcM0002--bcM0002,bcU0002_5p--bcU0002_3p,sample_5
bcM0003--bcM0003,bcU0001_5p--bcU0001_3p,sample_6
bcM0004--bcM0004,bcU0001_5p--bcU0001_3p,sample_1
```

The workflow groups rows by `HiFi Barcode` and derives one cDNA biosample CSV per
outer barcode:

For demultiplexing outputs, the workflow preserves first-seen `HiFi Barcode`
order from this CSV and cDNA barcode row order within each derived cDNA
biosample CSV.

`bcM0001--bcM0001.cdna_biosample.csv`:

```csv
Barcodes,Bio Sample
bcU0001_5p--bcU0001_3p,sample_1
bcU0002_5p--bcU0002_3p,sample_2
bcU0003_5p--bcU0003_3p,sample_3
```

`bcM0002--bcM0002.cdna_biosample.csv`:

```csv
Barcodes,Bio Sample
bcU0001_5p--bcU0001_3p,sample_4
bcU0002_5p--bcU0002_3p,sample_5
```

`bcM0003--bcM0003.cdna_biosample.csv`:

```csv
Barcodes,Bio Sample
bcU0001_5p--bcU0001_3p,sample_6
```

`bcM0004--bcM0004.cdna_biosample.csv`:

```csv
Barcodes,Bio Sample
bcU0001_5p--bcU0001_3p,sample_1
```

The repeated `sample_1` assignment is allowed because it occurs under different
outer barcodes. Those FLNC BAMs will carry the same `SM` value and can be merged
later during FLNC alignment.

### Repeated samples across outer barcodes

When the same biological samples appear under multiple outer barcodes, use the
same `Bio Sample` value in each group:

```csv
HiFi Barcode,cDNA Barcode,Bio Sample
bcM0001--bcM0001,IsoSeqX_bc01_5p--IsoSeqX_3p,sample_1
bcM0001--bcM0001,IsoSeqX_bc02_5p--IsoSeqX_3p,sample_2
bcM0001--bcM0001,IsoSeqX_bc03_5p--IsoSeqX_3p,sample_3
bcM0002--bcM0002,IsoSeqX_bc01_5p--IsoSeqX_3p,sample_1
bcM0002--bcM0002,IsoSeqX_bc02_5p--IsoSeqX_3p,sample_2
bcM0002--bcM0002,IsoSeqX_bc03_5p--IsoSeqX_3p,sample_3
```

Preprocessing produces one FLNC BAM per demultiplexed cDNA Barcode, but the BAM
headers carry the sample names from the `Bio Sample` column:

```text
m21003...bcM0001--bcM0001.IsoSeqX_bc01_5p--IsoSeqX_3p.flnc.bam  SM:sample_1
m21003...bcM0002--bcM0002.IsoSeqX_bc01_5p--IsoSeqX_3p.flnc.bam  SM:sample_1
m21003...bcM0001--bcM0001.IsoSeqX_bc02_5p--IsoSeqX_3p.flnc.bam  SM:sample_2
m21003...bcM0002--bcM0002.IsoSeqX_bc02_5p--IsoSeqX_3p.flnc.bam  SM:sample_2
```

The FLNC-level `secondary_analysis` groups BAMs by the exact `SM` tag before
alignment. A grouping report for the example above
would include rows like:

```tsv
sample_name sample_prefix group_size action input_bam
sample_1 sample_1 2 merge m21003...bcM0001--bcM0001.IsoSeqX_bc01_5p--IsoSeqX_3p.flnc.bam
sample_1 sample_1 2 merge m21003...bcM0002--bcM0002.IsoSeqX_bc01_5p--IsoSeqX_3p.flnc.bam
sample_2 sample_2 2 merge m21003...bcM0001--bcM0001.IsoSeqX_bc02_5p--IsoSeqX_3p.flnc.bam
sample_2 sample_2 2 merge m21003...bcM0002--bcM0002.IsoSeqX_bc02_5p--IsoSeqX_3p.flnc.bam
```

Every sample group is materialized with `samtools merge` into a single
`<sample_prefix>.flnc.bam`, including one-BAM groups. Each grouped sample is
then aligned once with `pbmm2`, producing outputs such as
`sample_1.aligned.bam`, `sample_2.aligned.bam`, and `sample_3.aligned.bam`.
