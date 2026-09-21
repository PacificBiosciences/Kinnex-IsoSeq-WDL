version 1.0

import "rna_structs.wdl"

workflow backend_configuration {
  meta {
    description: "Build the shared RNA runtime attributes object for supported HPC execution."
    outputs: {
      runtime_attributes: {
        description: "Runtime attribute structure"
      }
    }
  }

  parameter_meta {
    backend: {
      description: "Backend where the workflow will be executed; only HPC is supported"
    }
    max_retries: {
      description: "Maximum retries for failed task attempts"
    }
    add_memory_mb: {
      name: "Add Task Memory (MB)",
      description: "Increasing this number allocates extra memory per task when submitting jobs to the compute backend."
    }
    nproc: {
      description: "Maximum CPU threads per task for development and backend throttling",
      hidden: true
    }
    container_registry: {
      description: "Optional PacBio registry for registry-relative task images and named reference containers; if omitted, quay.io/pacbio is used"
    }
  }

  input {
    String backend = "HPC"
    Int max_retries = 2
    Int add_memory_mb = 0
    Int nproc = 32
    String? container_registry
  }

  String default_container_registry = "quay.io/pacbio"

  if (backend == "HPC") {
    RuntimeAttributes hpc_runtime_attributes = object {
      backend: "HPC",
      max_retries: max_retries,
      add_memory_mb: add_memory_mb,
      nproc: nproc,
      container_registry: select_first([
        container_registry,
        default_container_registry
      ])
    }
  }
  if (backend == "AWS-HealthOmics") {
    RuntimeAttributes aws_runtime_attributes = object {
      backend: "AWS-HealthOmics",
      max_retries: 1,
      add_memory_mb: add_memory_mb,
      nproc: nproc,
      container_registry: select_first([
        container_registry,
        default_container_registry
      ])
    }
  }

  output {
    RuntimeAttributes runtime_attributes = select_first([
      hpc_runtime_attributes,
      aws_runtime_attributes
    ])
  }
}
