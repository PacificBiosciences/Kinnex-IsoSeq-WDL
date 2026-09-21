# Resource bundle

The Kinnex Iso-Seq static resource bundle contains the GRCh38 genome,
annotation and Pigeon support files, plus Kinnex barcode, adapter, and primer
FASTAs. Workflow version `0.4.0` can resolve these files from a versioned
[reference container](./reference_container.md), typed per-file overrides, or
both.

The currently published packaged bundle is version `0.2.0` from the
[Kinnex Iso-Seq Zenodo dataset](https://zenodo.org/records/22255332).

## Workflow specification

All five reference-consuming entrypoints accept `ref_name`, an optional
explicit `reference_container`, and a default-empty `ReferenceOverrides`
object. `ref_name` defaults to `GRCh38_gencode49`; an explicit container URI
takes precedence. Each workflow uses only the reusable resources it needs:

- `hifi_demux_barcodes`;
- `segmentation_adapters`;
- `indexed_primers`;
- genome FASTA and its packaged FASTA index;
- annotation GTF;
- Pigeon polyA, CAGE, and junction resources.

When a required resource has an override, that file takes precedence over its
container counterpart. An overridden genome must be an uncompressed FASTA (the
workflow generates its FASTA index), while a container genome uses its packaged
index. A fully custom run omits the container and supplies all required
user-selectable files. `kinnex_isoseq` and `preprocessing` also accept optional
`segmentation_adapter_set` to select the packaged 8-, 12-, or 16-fold Kinnex
segmentation-adapter FASTA when `reference_overrides.segmentation_adapters` is
absent; omission selects `8-fold`. Their optional `isoseq_primers_set` selects
the packaged `IsoSeq-v2` or `IsoSeq96` indexed-primer FASTA when
`reference_overrides.indexed_primers` is absent; omission selects `IsoSeq-v2`.
Each selector conflicts with its matching custom override. Run-specific BAMs,
the biosample CSV, and tuning parameters remain ordinary workflow inputs.

See the [reference-container specification](./reference_container.md) for the full
manifest, URI requirements, provenance outputs, and current limitations.
