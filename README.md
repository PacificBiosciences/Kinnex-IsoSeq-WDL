<h1 align="center">PacBio Kinnex Iso-Seq Pipeline</h1>

> [!WARNING]
> **BETA:** This workflow is under active development. Interfaces, defaults,
> runtime behavior, and output bundles may change in wide-sweeping ways; use
> tagged releases for reproducible testing.

Workflow for analyzing human PacBio Kinnex Iso-Seq sequencing data using
[Workflow Description Language (WDL)](https://openwdl.org/).

- Docker images used by this workflow are defined in [the wdl-dockerfiles repo](https://github.com/PacificBiosciences/wdl-dockerfiles). Images are hosted in PacBio's [quay.io repository](https://quay.io/organization/pacbio).
- Task container images are selected in WDL runtime blocks and documented in the [tools and containers](docs/tools_containers.md) documentation.
- User-facing workflow entrypoints live in the [workflows directory](workflows).

## Workflow

This repository exposes three user-facing workflows:

The `kinnex_isoseq` workflow is designed for end-to-end analysis of Kinnex
Isoseq HiFi data. It takes as input HiFi bams from a single acquisition and
produces classified, filtered isoforms.

The `preprocessing` workflow is design to generate transcript BAMs, also known
as FLNC BAMs, from Kinnex HiFi data. It takes as input HiFi bams from a single
acquisition and generates sample-annotated transcript BAMs which contain poly-A
trimmed, 5' -> 3' transcript sequences.

The `secondary_analysis` workflow takes a set of sample-annotated transcript
BAMs, from any number of acquisitions. It produces classified, filtered isoform
for each sample, as annotated in the BAM header. If multiple BAMs are provided
for a sample, they are merged.

| <img src="docs/figures/workflow_diagram.png" width="800" alt="PacBio Kinnex Iso-Seq Pipeline"> |
| :---: |

**Workflow entrypoints**:

- [End-to-end HiFi BAM workflow](workflows/kinnex_isoseq.wdl)
- [Preprocessing workflow](workflows/preprocessing.wdl)
- [FLNC-level secondary analysis workflow](workflows/secondary_analysis.wdl)

## Quick Start Guide

In this section, we describe a quick way to get started with:

- `miniwdl` as the workflow engine,
- `slurm HPC` as the backend,
- `apptainer` as the container run time,
- the PacBio provided human reference data bundle,
- an instrument acquisition which was not demultiplexed on-instrument.

### Setup

- Clone this repository

```bash
git clone https://github.com/PacificBiosciences/Kinnex-IsoSeq-WDL.git
cd Kinnex-IsoSeq-WDL
```

- Install [miniwdl](docs/backend-hpc.md#miniwdl)
- Install [apptainer](https://apptainer.org/docs/admin/main/installation.html#install-from-pre-built-packages)
- Download and decompress the human reference data bundle from [Zenodo](https://zenodo.org/records/10839617)

### Fill out your inputs.json

Copy the `no_instrument_demux` kinnex_isoseq template:
[inputs.json](backends/hpc/kinnex_isoseq.no_instrument_demux.hpc.inputs.json)

Replace every `<local_path_prefix>` placeholder with paths visible to the
workflow jobs. Edit the referenced `ref_map_file` to the reference data bundle
you downloaded and decompressed.

### Run the workflow

```bash
miniwdl run --cfg ~/.config/miniwdl.cfg \
  -i kinnex_isoseq.no_instrument_demux.hpc.inputs.json \
  workflows/kinnex_isoseq.wdl
```

## Beyond Quick Start

1. [Select a backend environment](#selecting-a-backend)
2. [Configure a workflow execution engine and container runtime](#configuring-a-workflow-engine-and-container-runtime)
3. [Set up your reference](#setting-up-your-reference-datasets)
4. [Fill out the inputs JSON file for your run](#filling-out-the-inputs-json)
5. [Run the workflow](#running-the-workflow)

### Selecting a backend

The currently supported public backend target is Slurm HPC. We provide
maintained example configurations for running the workflow with miniwdl through
miniwdl-slurm, and with Sprocket with its Slurm/Apptainer backend. Cromwell
can use the same WDL entrypoints and HPC input templates, but this repository
does not provide a Cromwell backend configuration.

| Backend | Status | Documentation |
| :- | :- | :- |
| HPC | Supported on Slurm with miniwdl, Sprocket, or Cromwell | [HPC backend guide](docs/backend-hpc.md) |
| AWS HealthOmics | No maintained template | Not currently provided |
| Azure | No maintained template | Not currently provided |
| GCP | No maintained template | Not currently provided |

### Configuring a workflow engine and container runtime

An execution engine is required to run WDL workflows. For the supported public
HPC path, install [`miniwdl`](https://miniwdl.readthedocs.io/en/latest/) with
[`miniwdl-slurm`](https://github.com/miniwdl-ext/miniwdl-slurm),
[`sprocket`](https://github.com/stjude-rust-labs/sprocket), or [`Cromwell`](https://cromwell.readthedocs.io/en/stable/backends/Backends/).
Make sure Slurm commands and a Singularity-compatible container runtime are available
on the submit host and cluster nodes. See the [HPC backend guide](docs/backend-hpc.md)
for setup details.

### Setting up your reference datasets

The [resource bundle](docs/resource_bundle.md) contains both reference
and preprocessing resources and can be downloaded from Zenodo: [Zenodo record](https://zenodo.org/records/10839617).
Some resources are selected through the [reference-map TSV](docs/ref_map.md), including the
reference genome, annotation, and classification support files.
Other resources are provided as explicit workflow inputs, including
the preprocessing barcode, adapter, and primer FASTA files.

Start from the [HPC reference-map template](backends/hpc/GRCh38.ref_map.v0p1p0.hpc.tsv) when
running on a Slurm cluster. See the [resource bundle layout](docs/resource_bundle.md)
for the expected directory structure and path-placeholder conventions.

### Filling out the inputs JSON

The input to a workflow run is defined in JSON format. Public workflow input
templates live next to the workflow entrypoints:

- [workflows/kinnex_isoseq.inputs.json](workflows/kinnex_isoseq.inputs.json)
- [workflows/preprocessing.inputs.json](workflows/preprocessing.inputs.json)
- [workflows/secondary_analysis.inputs.json](workflows/secondary_analysis.inputs.json)

HPC input templates are available in the [HPC backend directory](backends/hpc):

- [backends/hpc/kinnex_isoseq.no_instrument_demux.hpc.inputs.json](backends/hpc/kinnex_isoseq.no_instrument_demux.hpc.inputs.json)
- [backends/hpc/kinnex_isoseq.from_instrument_demux.hpc.inputs.json](backends/hpc/kinnex_isoseq.from_instrument_demux.hpc.inputs.json)
- [backends/hpc/preprocessing.hpc.inputs.json](backends/hpc/preprocessing.hpc.inputs.json)
- [backends/hpc/secondary_analysis.hpc.inputs.json](backends/hpc/secondary_analysis.hpc.inputs.json)
- [backends/hpc/GRCh38.ref_map.v0p1p0.hpc.tsv](backends/hpc/GRCh38.ref_map.v0p1p0.hpc.tsv)

For the end-to-end workflow, use the `no_instrument_demux` template for one raw
HiFi BAM, or the `from_instrument_demux` template for BAMs that have been
demultiplexed using the HiFi barcode.

Copy the template for the entrypoint you want to run, replace every
`<local_path_prefix>` placeholder with paths visible to the workflow jobs, and
copy or edit the referenced [`ref_map_file`](docs/ref_map.md) so its reference, annotation, and
classification support-resource paths point to files on your filesystem. See
the [resource bundle layout](docs/resource_bundle.md) for the split between
files referenced by [`ref_map_file`](docs/ref_map.md) and files passed directly as preprocessing
inputs.

### Running the workflow

Run with `miniwdl` from the repository root using the input JSON you copied and
edited. See the [HPC backend guide](docs/backend-hpc.md) for additional
engine-specific execution notes.

End-to-end HiFi BAM entrypoint, using the `no_instrument_demux` HPC template:

```bash
miniwdl run --cfg ~/.config/miniwdl.cfg \
  -i backends/hpc/kinnex_isoseq.no_instrument_demux.hpc.inputs.json \
  workflows/kinnex_isoseq.wdl
```

Preprocessing entrypoint:

```bash
miniwdl run --cfg ~/.config/miniwdl.cfg \
  -i backends/hpc/preprocessing.hpc.inputs.json \
  workflows/preprocessing.wdl
```

FLNC-level secondary analysis entrypoint:

```bash
miniwdl run --cfg ~/.config/miniwdl.cfg \
  -i backends/hpc/secondary_analysis.hpc.inputs.json \
  workflows/secondary_analysis.wdl
```

## Workflow inputs

At a high level, the workflows use three input-file categories:

- *maps* are TSV files that describe shared reference and annotation inputs used by one or more stages. The primary map is [`ref_map_file`](docs/ref_map.md).
- *sample sheets* are CSV files that describe run-specific sample and barcode assignments. `kinnex_isoseq` and `preprocessing` use [`biosample_csv`](docs/biosample_csv.md).
- *inputs.json* files define the datasets, sample files, run mode, tuning parameters, and backend settings for a workflow run.

The end-to-end `kinnex_isoseq` entrypoint starts from `hifi_sources`.
Omit `hifi_barcode` for HiFi demux mode, or provide one `hifi_barcode` per
source BAM for cDNA demux-only mode.

The standalone `preprocessing` entrypoint uses the same `hifi_sources`,
preprocessing resource, and `biosample_csv` inputs as the end-to-end workflow,
and stops after producing FLNC BAMs.

The FLNC-level `secondary_analysis` entrypoint starts from `flnc_bams`. Each
FLNC BAM must contain `@RG` records with exactly one distinct non-empty `SM`
value. BAMs are grouped by the exact original `SM`; sanitized sample prefixes
are used for grouped FLNC and aligned BAM filenames.

Detailed workflow input and output interfaces are described in:

- [End-to-end HiFi BAM workflow interface](docs/kinnex_isoseq.md)
- [Preprocessing workflow interface](docs/preprocessing.md)
- [FLNC-level secondary analysis workflow interface](docs/secondary_analysis.md)
- [Resource bundle layout](docs/resource_bundle.md)
- [Reference-map TSV specification](docs/ref_map.md)
- [Biosample CSV specification](docs/biosample_csv.md)

## Tool versions and Docker images

Task runtime blocks define the container images used by each workflow task.
Registry-relative PacBio image references use the workflow `container_registry`
input when provided, and otherwise default to `"quay.io/pacbio"`.

Tool and container details are documented in [tools and containers](docs/tools_containers.md).

## Resource requirements

Resource requirements depend on input size, sample count, and the entrypoint
being run. Task-level CPU, base memory, retry, and image settings are defined in
WDL runtime blocks.

## Version information

Current version: **0.2.0**.

For a complete changelog, see the [changelog](CHANGELOG.md) or the git history.

## DISCLAIMER

TO THE GREATEST EXTENT PERMITTED BY APPLICABLE LAW, THIS WEBSITE AND ITS
CONTENT, INCLUDING ALL SOFTWARE, SOFTWARE CODE, SITE-RELATED SERVICES, AND DATA,
ARE PROVIDED "AS IS," WITH ALL FAULTS, WITH NO REPRESENTATIONS OR WARRANTIES OF
ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, ANY
WARRANTIES OF MERCHANTABILITY, SATISFACTORY QUALITY, NON-INFRINGEMENT OR FITNESS
FOR A PARTICULAR PURPOSE. ALL WARRANTIES ARE REJECTED AND DISCLAIMED. YOU ASSUME
TOTAL RESPONSIBILITY AND RISK FOR YOUR USE OF THE FOREGOING. PACBIO IS NOT
OBLIGATED TO PROVIDE ANY SUPPORT FOR ANY OF THE FOREGOING, AND ANY SUPPORT
PACBIO DOES PROVIDE IS SIMILARLY PROVIDED WITHOUT REPRESENTATION OR WARRANTY OF
ANY KIND. NO ORAL OR WRITTEN INFORMATION OR ADVICE SHALL CREATE A REPRESENTATION
OR WARRANTY OF ANY KIND. ANY REFERENCES TO SPECIFIC PRODUCTS OR SERVICES ON THE
WEBSITES DO NOT CONSTITUTE OR IMPLY A RECOMMENDATION OR ENDORSEMENT BY PACBIO.
