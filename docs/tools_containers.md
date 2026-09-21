# Tool versions and containers

Containers package the tools and dependencies needed to run the workflow
reproducibly with Docker, Singularity, or Podman.

PacBio registry-relative images, including named reference containers, use the
workflow-level `container_registry` input when provided and default to
`"quay.io/pacbio"`.

Production task images are pinned by digest for reproducibility and better
compatibility with Cromwell and miniwdl call caching. A task's `runtime` block
identifies its Docker image with the `docker` key.

## Analysis-tool containers

This table lists the production containers for tools directly invoked by the
Kinnex Iso-Seq WDL tasks.

| Container | Major tool versions | Dockerfile | Container |
| --------: | ------------------- | :---: | :---: |
| `pb_wdl_base` | <ul><li>[python3 3.14.5](https://github.com/python/cpython/tree/v3.14.5)</li><li>[pysam 0.24.0](https://github.com/pysam-developers/pysam/releases/tag/v0.24.0)</li><li>[samtools 1.23.1](https://github.com/samtools/samtools/releases/tag/1.23.1)</li><li>container tag: build4</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/f3817ebb378256a87f95f2388a1a340b8b327c91/docker/pb_wdl_base) | [`pb_wdl_base@sha256:03cb3c01937eccc907f8ad71c87b258581504572205fe3f31a657e318f3564ae`](https://quay.io/repository/pacbio/pb_wdl_base/manifest/sha256:03cb3c01937eccc907f8ad71c87b258581504572205fe3f31a657e318f3564ae) |
| `pbmm2` | <ul><li>[pbmm2 26.2.0](https://github.com/PacificBiosciences/pbmm2/releases/tag/v26.2.0)</li><li>container tag: 26.2.0_build1</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/01afb4104fc5c6fa1d43984ef5a8d557f3ab3ae5/docker/pbmm2) | [`pbmm2@sha256:0c21f29f1ee429dbafe5c332d4abffcfed5efbf5448b4ef06dd189b5323a7051`](https://quay.io/repository/pacbio/pbmm2/manifest/sha256:0c21f29f1ee429dbafe5c332d4abffcfed5efbf5448b4ef06dd189b5323a7051) |
| `pbsamoa` | <ul><li>[pbsamoa 20260811](https://github.com/PacificBiosciences/pbsamoa/releases/tag/20260811)</li><li>container tag: 20260811_build1</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/334f27df4ebae330d4469964a2abe187a396aa6f/docker/pbsamoa) | [`pbsamoa@sha256:cf70c89422d63c3a4e4b2ffd912160c3e7481303a211cd6134236fe6e9f9461f`](https://quay.io/repository/pacbio/pbsamoa/manifest/sha256:cf70c89422d63c3a4e4b2ffd912160c3e7481303a211cd6134236fe6e9f9461f) |
| `isocall` | <ul><li>[isocall 1.3.0](https://github.com/PacificBiosciences/isocall/releases/tag/v1.3.0)</li><li>container tag: 1.3.0_build1</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/e3777cff54bc8e4833a96897092e0f8475b276b2/docker/isocall) | [`isocall@sha256:34cbf179e342515ddfdd65e7f170ab76e4f6fd3993fcd06f236c13712eb8e882`](https://quay.io/repository/pacbio/isocall/manifest/sha256:34cbf179e342515ddfdd65e7f170ab76e4f6fd3993fcd06f236c13712eb8e882) |
| `lima` | <ul><li>[lima 26.2.1](https://github.com/PacificBiosciences/barcoding/releases/tag/v26.2.1)</li><li>container tag: 26.2.1_build3</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/b95da0667aef8122febfbed3905cef9178d69512/docker/lima) | [`lima@sha256:c49a067e447fdd2a13b4ebb41cfd66d1242032233671dbe4a4fdbe3e62ac82ef`](https://quay.io/repository/pacbio/lima/manifest/sha256:c49a067e447fdd2a13b4ebb41cfd66d1242032233671dbe4a4fdbe3e62ac82ef) |
| `skera` | <ul><li>[skera 1.4.0](https://github.com/PacificBiosciences/skera/releases/tag/v1.4.0)</li><li>container tag: 1.4.0_build3</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/3192d1e7823d7f62259db4e359d59206f6a97e01/docker/skera) | [`skera@sha256:f823bd64f2beec351cae82e67c9c6f257b4896ddd59f8b2050568d59f3a165f6`](https://quay.io/repository/pacbio/skera/manifest/sha256:f823bd64f2beec351cae82e67c9c6f257b4896ddd59f8b2050568d59f3a165f6) |
| `isoseq` | <ul><li>[isoseq 26.2.0](https://github.com/PacificBiosciences/IsoSeq/releases/tag/v26.2.0)</li><li>container tag: 26.2.0_build2</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/0cd11b56bc46b9aa42010f6cb3a713ac976325eb/docker/isoseq) | [`isoseq@sha256:d5e4bfd5dd570b26b9769b7b340beeefd618324685662ce2fc898b9a817e465e`](https://quay.io/repository/pacbio/isoseq/manifest/sha256:d5e4bfd5dd570b26b9769b7b340beeefd618324685662ce2fc898b9a817e465e) |
| `pbcore` | <ul><li>[python3 3.14.5](https://github.com/python/cpython/tree/v3.14.5)</li><li>[pbcore 2.6.0](https://github.com/PacificBiosciences/pbcore/tree/2.6.0)</li><li>container tag: 2.6.0_build3</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/cd3ba5f0f8710d72cde07092774a503df92d86ac/docker/pbcore) | [`pbcore@sha256:47fc6f1174605be9a8a932f17c15d1a29fc5e6b89ea6809f32995318c7afe1e7`](https://quay.io/repository/pacbio/pbcore/manifest/sha256:47fc6f1174605be9a8a932f17c15d1a29fc5e6b89ea6809f32995318c7afe1e7) |
| `pigeon` | <ul><li>[pigeon 26.2.0](https://github.com/PacificBiosciences/pigeon/releases/tag/v26.2.0)</li><li>container tag: 26.2.0_build2</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/e56e4ea36077d92ac46240e31a8b07dc4deeed6b/docker/pigeon) | [`pigeon@sha256:94cf3cd64f600ae974a7956177be2231c1c9ecb468a1cc4b70b05dcacf59a81c`](https://quay.io/repository/pacbio/pigeon/manifest/sha256:94cf3cd64f600ae974a7956177be2231c1c9ecb468a1cc4b70b05dcacf59a81c) |

## Reference containers

The [reference container](./reference_container.md) is normally selected with
`ref_name: "GRCh38_gencode49"`. This registry-relative selection follows
`container_registry`. You can instead supply a complete digest-pinned URI,
which `container_registry` does not rewrite. Reference containers package
static workflow data rather than analysis tools. Typed filesystem overrides can
replace any packaged reference resource.

The current GRCh38 container was built from the
[Kinnex Iso-Seq reference data bundle on Zenodo](https://zenodo.org/records/22255332).

| Reference bundle | Packaged resources | Container |
| --- | --- | --- |
| Kinnex Iso-Seq GRCh38 `0.2.0` | GRCh38 genome and annotation, Pigeon support files, and Kinnex barcode, adapter, and primer FASTAs | [`workflow-data-container-kinnex-isoseq-wdl-grch38@sha256:c68ba123e395f41176dd7a8d500bc84a7d1d097b08a850dc8ff57664a67c811a`](https://quay.io/repository/pacbio/workflow-data-container-kinnex-isoseq-wdl-grch38/manifest/sha256:c68ba123e395f41176dd7a8d500bc84a7d1d097b08a850dc8ff57664a67c811a) |
