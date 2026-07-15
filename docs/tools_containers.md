# Tool versions and containers

Containers are used to package tools and their dependencies so that the
workflow can run reproducibly on systems with a container runtime, including
Docker, Singularity, or Podman.

PacBio registry-relative images use the workflow-level `container_registry`
input when provided, and default to `"quay.io/pacbio"`.

Production task images are pinned by digest for reproducibility and to improve
compatibility with Cromwell and miniwdl call caching. The Docker image used by a
particular task can be identified from the task `runtime` block's `docker` key.

This table lists the production containers for tools directly invoked by the
Kinnex Iso-Seq WDL tasks.

| Container | Major tool versions | Dockerfile | Container |
| --------: | ------------------- | :---: | :---: |
| `pb_wdl_base` | <ul><li>python3 3.14.5</li><li>pysam 0.24.0</li><li>samtools 1.23.1</li><li>container tag: build4</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/f3817ebb378256a87f95f2388a1a340b8b327c91/docker/pb_wdl_base) | [`pb_wdl_base@sha256:03cb3c01937eccc907f8ad71c87b258581504572205fe3f31a657e318f3564ae`](https://quay.io/repository/pacbio/pb_wdl_base/manifest/sha256:03cb3c01937eccc907f8ad71c87b258581504572205fe3f31a657e318f3564ae) |
| `pbmm2` | <ul><li>pbmm2 26.2.0</li><li>container tag: 26.2.0_build1</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/01afb4104fc5c6fa1d43984ef5a8d557f3ab3ae5/docker/pbmm2) | [`pbmm2@sha256:0c21f29f1ee429dbafe5c332d4abffcfed5efbf5448b4ef06dd189b5323a7051`](https://quay.io/repository/pacbio/pbmm2/manifest/sha256:0c21f29f1ee429dbafe5c332d4abffcfed5efbf5448b4ef06dd189b5323a7051) |
| `isocall` | <ul><li>isocall 1.1.0</li><li>container tag: 1.1.0_build1</li></ul> | TBD | [`isocall@sha256:1d45a7256f2f5e172b4722473d6feb604d104f694840c5ab7b4d4d5202b00c9b`](https://quay.io/repository/pacbio/isocall/manifest/sha256:1d45a7256f2f5e172b4722473d6feb604d104f694840c5ab7b4d4d5202b00c9b) |
| `lima` | <ul><li>lima 26.2.1</li><li>container tag: 26.2.1_build3</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/b95da0667aef8122febfbed3905cef9178d69512/docker/lima) | [`lima@sha256:c49a067e447fdd2a13b4ebb41cfd66d1242032233671dbe4a4fdbe3e62ac82ef`](https://quay.io/repository/pacbio/lima/manifest/sha256:c49a067e447fdd2a13b4ebb41cfd66d1242032233671dbe4a4fdbe3e62ac82ef) |
| `skera` | <ul><li>skera 1.4.0</li><li>container tag: 1.4.0_build3</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/3192d1e7823d7f62259db4e359d59206f6a97e01/docker/skera) | [`skera@sha256:f823bd64f2beec351cae82e67c9c6f257b4896ddd59f8b2050568d59f3a165f6`](https://quay.io/repository/pacbio/skera/manifest/sha256:f823bd64f2beec351cae82e67c9c6f257b4896ddd59f8b2050568d59f3a165f6) |
| `isoseq` | <ul><li>isoseq 26.2.0</li><li>container tag: 26.2.0_build2</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/0cd11b56bc46b9aa42010f6cb3a713ac976325eb/docker/isoseq) | [`isoseq@sha256:d5e4bfd5dd570b26b9769b7b340beeefd618324685662ce2fc898b9a817e465e`](https://quay.io/repository/pacbio/isoseq/manifest/sha256:d5e4bfd5dd570b26b9769b7b340beeefd618324685662ce2fc898b9a817e465e) |
| `pbcore` | <ul><li>python3 3.14.5</li><li>pbcore 2.6.0</li><li>container tag: 2.6.0_build3</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/cd3ba5f0f8710d72cde07092774a503df92d86ac/docker/pbcore) | [`pbcore@sha256:47fc6f1174605be9a8a932f17c15d1a29fc5e6b89ea6809f32995318c7afe1e7`](https://quay.io/repository/pacbio/pbcore/manifest/sha256:47fc6f1174605be9a8a932f17c15d1a29fc5e6b89ea6809f32995318c7afe1e7) |
| `pigeon` | <ul><li>pigeon 26.2.0</li><li>container tag: 26.2.0_build2</li></ul> | [Dockerfile](https://github.com/PacificBiosciences/wdl-dockerfiles/tree/e56e4ea36077d92ac46240e31a8b07dc4deeed6b/docker/pigeon) | [`pigeon@sha256:94cf3cd64f600ae974a7956177be2231c1c9ecb468a1cc4b70b05dcacf59a81c`](https://quay.io/repository/pacbio/pigeon/manifest/sha256:94cf3cd64f600ae974a7956177be2231c1c9ecb468a1cc4b70b05dcacf59a81c) |
