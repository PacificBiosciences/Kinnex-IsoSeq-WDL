# Resource bundle

The Kinnex Iso-Seq static resource bundle contains the GRCh38 genome,
annotation and Pigeon support files, plus Kinnex barcode, adapter, and primer
FASTAs. Workflow version `0.3.0` can resolve these files from a versioned
[reference container](./reference_container.md), typed per-file overrides, or
both.

The packaged bundle is version `0.1.0` from the
[Kinnex Iso-Seq Zenodo dataset](https://doi.org/10.5281/zenodo.10839617).

## Workflow contract

All three user-facing workflows—`kinnex_isoseq`, `preprocessing`, and
`secondary_analysis`—accept an optional `reference_container` and a
default-empty `ReferenceOverrides` object. A subset of
these reusable resources are used by the selected workflow:

- `hifi_demux_barcodes`;
- `skera_adapters`;
- `barcoded_primers`;
- genome FASTA and its packaged FASTA index;
- annotation GTF;
- Pigeon polyA, CAGE, and junction resources.

When a required resource has an override, that file takes precedence over its
container counterpart. An overridden genome must be an uncompressed FASTA (the
workflow generates its FASTA index), while a container genome uses its packaged
index. A fully custom run omits the container and supplies all required
user-selectable files. `kinnex_isoseq` and `preprocessing` additionally accept
optional `kinnex_primers_set` to select the packaged 8-, 12-, or 16-fold
Kinnex primer/Skera adapter FASTA when `reference_overrides.skera_adapters` is
absent; omission selects `8fold`. Run-specific BAMs, the biosample CSV, optional dataset XMLs,
and tuning parameters remain ordinary workflow inputs.

See the [reference-container contract](./reference_container.md) for the full
manifest, URI requirements, provenance outputs, and current limitations.
