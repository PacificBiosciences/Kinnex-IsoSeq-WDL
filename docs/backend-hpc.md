# HPC backend execution

The supported HPC target is a Slurm cluster. The repository includes maintained
example configurations for `miniwdl` with `miniwdl-slurm` and for Sprocket with
its Slurm/Apptainer backend.

Select the workflow engine outside WDL by choosing its configuration. The
public miniwdl configuration is `backends/hpc/miniwdl.cfg`, and the Sprocket
configuration is `backends/hpc/sprocket.toml`. Cromwell users must provide
their own configuration.

Execution is container based, and workflow inputs set the shared backend
options. Every engine needs access to Slurm submission commands and a
Singularity-compatible container runtime on the cluster nodes.

> **Slurm launch host:** Start the workflow engine from a login node or submit
> node with access to Slurm submission commands and the shared filesystem. The
> engine itself submits and monitors workflow tasks through Slurm, so do not
> start it from inside an interactive `srun` session.

## Install a workflow engine

### miniwdl

Install `miniwdl` and the Slurm plugin on the submit host where you start
workflow runs. Use `uv` when it is available on your cluster or can be
installed. It creates an isolated environment and installs the Slurm plugin in
that environment.

```bash
uv tool install --python 3.12 'miniwdl>=1.13.1,<2' --with miniwdl-slurm
miniwdl --version
```

If `miniwdl` is not on `PATH` after installation, run `uv tool update-shell`
and start a new shell session, or add the directory shown by
`uv tool dir --bin` to your `PATH` using your site's preferred shell setup.

If your cluster does not support `uv`, use an isolated Python or Mamba/Conda
environment. For example:

```bash
python3 -m venv .venv-miniwdl
source .venv-miniwdl/bin/activate
python -m pip install --upgrade pip
python -m pip install 'miniwdl>=1.13.1,<2' miniwdl-slurm
```

### Sprocket

Install Sprocket on the host where you start workflow runs.

To set up Sprocket, either download a release binary from
<https://github.com/stjude-rust-labs/sprocket/releases>, or use the Rust
toolchain to build it yourself. For example, with Cargo:

```bash
cargo install sprocket --version 0.30.1 --locked
sprocket --version
```

Building Sprocket 0.30.1 with Cargo requires Rust 1.95 or newer.

### Cromwell

Install Cromwell and Java on the host or service that starts workflow runs.

```bash
java -jar /path/to/cromwell.jar --version
```

## Configuring an engine for Slurm

### miniwdl

Start from the example configuration:

```bash
mkdir -p ~/.config
cp backends/hpc/miniwdl.cfg ~/.config/miniwdl.cfg
```

Edit `~/.config/miniwdl.cfg` for your cluster before running.
At minimum, review these settings:

- `[scheduler] task_concurrency`: upper bound on concurrently monitored tasks.
- `[singularity] exe`: path to the container runtime on your cluster.
- `[singularity] run_options`: site-specific container isolation and GPU flags.
- `[singularity] image_cache`: shared cache location reachable by submit and
  compute nodes.
- `[slurm] extra_args`: partition, account, QoS, reservation, and other
  cluster-specific submission flags.
- `[call_cache] dir`: call-cache location for repeat or resumed runs.

Keep the following setting enabled:

```ini
[file_io]
allow_any_input = true
```

Public and internal examples may use input files from shared locations outside
the workflow run directory, so the maintained miniwdl configuration keeps
`allow_any_input = true`.

### Sprocket

> **Warning:** Sprocket is under active development. If you use the example
> `backends/hpc/sprocket.toml` configuration from this repository, use it with
> the Sprocket version recommended here.

Start from the public example configuration:

```bash
cp backends/hpc/sprocket.toml sprocket.hpc.toml
```

Edit `sprocket.hpc.toml` for your cluster before running. At
minimum, review these settings:

- `[run] output_dir`: Sprocket run output directory.
- `[server.database] url`: Sprocket SQLite provenance database for run records.
- `[run.task] cache_dir`: Sprocket task call-cache directory.
- `[run.task] digests`: content-digest strategy used for task call caching.
- `[run.http] cache_dir`: Sprocket HTTP download cache directory.
- `[run.backends.default] default_slurm_partition.name`: default Slurm
  partition for workflow tasks.
- `[run.backends.default] apptainer.executable`: path or command name for
  Apptainer or Singularity.
- `[run.backends.default] apptainer.extra_args`: site-specific container
  isolation flags.
- `[run.backends.default] sbatch.args`: account, QoS, reservation, and
  other cluster-specific submission flags.
- `[run.backends.default] apptainer.image_cache_dir`: shared `.sif` image cache
  location reachable by submit and compute nodes.

The output, provenance database, call-cache, HTTP-cache, and image-cache paths
are relative to the directory where you run `sprocket`. Run from the repository
root as shown below, or edit them to absolute shared paths. Cache formats are
engine specific and should not be shared between engines.

The public Sprocket template enables task call caching by default. Review
`[run.task] cache_dir` before running so repeat or resumed runs use a cache
location that is visible from the submit host and compute nodes. It uses
`digests = "strongish"`, which hashes file metadata and the first 10 MiB of
each file to improve cache invalidation without fully hashing large input
files.

The Sprocket Slurm/Apptainer backend currently requires experimental execution
features. Enable them with `experimental_features_enabled = true`.

All filesystem inputs, including typed reference overrides, must be visible
from the submit host and compute nodes through the same shared paths. When a
reference container supplies any defaults, compute nodes must also be able to
authenticate to and pull it.

### Cromwell

This repository does not provide or maintain a Cromwell configuration. Use
your own Cromwell configuration and workflow options to define the Slurm
backend, container runtime, localization behavior, call-cache settings, and
output directories.

You can use the `.hpc.inputs.json` templates with Cromwell after replacing the
file paths. Keep the workflow input `backend` set to `"HPC"`. This
workflow-level runtime attribute selector does not need to match the backend
name in your Cromwell configuration.

## Prepare inputs

Public workflow input templates live next to the workflow entrypoints:

- [workflows/kinnex_isoseq.inputs.json](../workflows/kinnex_isoseq.inputs.json)
- [workflows/preprocessing.inputs.json](../workflows/preprocessing.inputs.json)
- [workflows/secondary_analysis.inputs.json](../workflows/secondary_analysis.inputs.json)

HPC templates live under [backends/hpc](../backends/hpc):

- [backends/hpc/kinnex_isoseq.no_instrument_demux.hpc.inputs.json](../backends/hpc/kinnex_isoseq.no_instrument_demux.hpc.inputs.json)
- [backends/hpc/kinnex_isoseq.from_instrument_demux.hpc.inputs.json](../backends/hpc/kinnex_isoseq.from_instrument_demux.hpc.inputs.json)
- [backends/hpc/preprocessing.hpc.inputs.json](../backends/hpc/preprocessing.hpc.inputs.json)
- [backends/hpc/secondary_analysis.hpc.inputs.json](../backends/hpc/secondary_analysis.hpc.inputs.json)
- [backends/hpc/biosamples.example.csv](../backends/hpc/biosamples.example.csv),
  a 4 × 12 barcode-combination sample-sheet template for preprocessing and
  end-to-end runs

Copy the matching `.hpc.inputs.json` template and replace every
`<local_path_prefix>` placeholder with paths visible from the Slurm jobs.
The HPC templates demonstrate container-only reference resolution with
`ref_name: "GRCh38_gencode49"`. The shared resolver maps that name to the
immutable container published for the workflow release, so compute nodes need
registry access. For hybrid or fully custom references, add the typed
`reference_overrides` object. The default named container is not unpacked when
every effective resource is overridden. See the
[reference-resource contract](./reference_container.md).

For `preprocessing` or `kinnex_isoseq` runs, prepare the
[`biosample_csv`](./biosample_csv.md) sample sheet for your run. Use the
checked-in [`biosamples.example.csv`](../backends/hpc/biosamples.example.csv)
as a starting point. Set the optional `segmentation_adapter_set` input to `8-fold`,
`12-fold`, or `16-fold` to select a packaged Kinnex segmentation-adapter FASTA;
omission selects `8-fold`. Set `isoseq_primers_set` to `IsoSeq-v2` or
`IsoSeq96` to select packaged indexed primers; omission selects `IsoSeq-v2`.
Each selector conflicts with its matching custom override,
`reference_overrides.segmentation_adapters` or
`reference_overrides.indexed_primers`.

Each preprocessing resource comes from its override or the reference
container.

For the end-to-end `kinnex_isoseq` workflow, choose the mode-specific template:
use `kinnex_isoseq.no_instrument_demux.hpc.inputs.json` for one raw HiFi BAM
that still needs upstream HiFi demux, or
`kinnex_isoseq.from_instrument_demux.hpc.inputs.json` for already HiFi-demuxed
BAMs that each provide a `hifi_barcode`.

Each public template includes `backend: "HPC"`. If a run needs more memory, add
the workflow's `add_memory_mb` input to your input JSON to allocate extra memory
per task when submitting jobs to the compute backend. Optionally set
`container_registry` in your input JSON if your site mirrors registry-relative
PacBio images to a different registry; omit it to use `"quay.io/pacbio"`. Do
not add per-task Docker or resource blocks to the public templates.
`container_registry` applies to named reference containers but does not rewrite
an explicit `reference_container` URI or filesystem override paths.

Review the [tools and containers](./tools_containers.md) documentation before
executing pure-container runs.

## Run workflows

Run from the repository root, using the input JSON you copied and edited.

### miniwdl

Omit `--cfg` if the configuration is installed at `~/.config/miniwdl.cfg`.

```bash
miniwdl run --cfg ~/.config/miniwdl.cfg \
  -i backends/hpc/kinnex_isoseq.no_instrument_demux.hpc.inputs.json \
  workflows/kinnex_isoseq.wdl
```

### Sprocket

Pass input JSON files to Sprocket with an `@` prefix. Sprocket accepts the same
HPC templates as miniwdl.

```bash
sprocket run --config sprocket.hpc.toml \
  workflows/kinnex_isoseq.wdl \
  @backends/hpc/kinnex_isoseq.no_instrument_demux.hpc.inputs.json
```

### Cromwell

Cromwell can run the same WDL entrypoints and HPC input templates when supplied
with a Cromwell backend configuration and workflow options.

```bash
java -Dconfig.file=/path/to/cromwell.conf \
  -jar /path/to/cromwell.jar run \
  workflows/kinnex_isoseq.wdl \
  -i backends/hpc/kinnex_isoseq.no_instrument_demux.hpc.inputs.json \
  -o /path/to/cromwell.options.json
```

Use the same command pattern for the other user-facing entrypoints listed
above.
