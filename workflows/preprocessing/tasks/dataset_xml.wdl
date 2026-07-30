version 1.0

import "../../rna_structs.wdl"

task populate_flnc_dataset_xml {
  meta {
    description: "Package final FLNC BAM/PBI files and create a top-level FLNC ConsensusReadSet XML from a ConsensusReadSet XML metadata template."
    outputs: {
      flnc_dataset_xml: {
        description: "FLNC ConsensusReadSet XML"
      },
      child_flnc_dataset_xmls: {
        description: "Child FLNC ConsensusReadSet XMLs"
      },
      packaged_flnc_bams: {
        description: "Packaged FLNC BAMs"
      },
      packaged_flnc_bam_pbis: {
        description: "Packaged FLNC BAM PBIs"
      }
    }
  }

  parameter_meta {
    consensusreadset_xml: {
      description: "ConsensusReadSet XML template"
    }
    collection_contexts: {
      description: "Collection contexts"
    }
    validated_biosample_csv: {
      description: "Validated biosample CSV"
    }
    flnc_bams: {
      description: "FLNC BAMs"
    }
    flnc_bam_pbis: {
      description: "FLNC BAM PBIs"
    }
    runtime_attributes: {
      description: "Runtime attribute structure"
    }
  }

  input {
    File consensusreadset_xml
    File collection_contexts
    File validated_biosample_csv
    Array[File] flnc_bams
    Array[File] flnc_bam_pbis
    RuntimeAttributes runtime_attributes
  }

  Int threads = 1
  Int mem_gb = 4
  Int total_mem_mb = (mem_gb * 1024) + runtime_attributes.add_memory_mb

  command <<<
    set -euo pipefail

    printf '%s\n' "~{sep="\" \"" flnc_bams}" > flnc_bams.txt
    printf '%s\n' "~{sep="\" \"" flnc_bam_pbis}" > flnc_bam_pbis.txt

    python3 - \
      "~{consensusreadset_xml}" \
      "~{collection_contexts}" \
      "~{validated_biosample_csv}" \
      "flnc_bams.txt" \
      "flnc_bam_pbis.txt" <<'PY'
    import copy
    import csv
    import errno
    import os
    import re
    import shutil
    import sys

    from pbcore.io import BamReader, ConsensusReadSet
    from pbcore.io.dataset.DataSetMembers import (
        BioSampleMetadata,
        BioSamplesMetadata,
        DNABarcode,
        DNABarcodes,
        ExternalResource,
    )

    (
        consensusreadset_xml,
        collection_contexts_file,
        validated_biosample_csv,
        flnc_bams_file,
        flnc_bam_pbis_file,
    ) = sys.argv[1:]
    dataset_xml = 'flnc.consensusreadset.xml'

    with open(collection_contexts_file, 'r', encoding='utf-8') as fh:
        collection_contexts = {line.rstrip('\n') for line in fh if line.rstrip('\n')}
    with open(flnc_bams_file, 'r', encoding='utf-8') as fh:
        flnc_bams = [line.rstrip('\n') for line in fh if line.rstrip('\n')]
    with open(flnc_bam_pbis_file, 'r', encoding='utf-8') as fh:
        flnc_bam_pbis = [line.rstrip('\n') for line in fh if line.rstrip('\n')]

    if not collection_contexts:
        sys.exit('validated consensusreadset XML collection context set is empty')
    if not flnc_bams:
        sys.exit('flnc_bams must contain at least one FLNC BAM')
    if len(flnc_bams) != len(flnc_bam_pbis):
        sys.exit(
            f'flnc_bams and flnc_bam_pbis must contain the same number of files: {len(flnc_bams)} != {len(flnc_bam_pbis)}'
        )


    def exactly_one(values, what, bam_path):
        if len(values) != 1:
            observed = ', '.join(repr(value) for value in sorted(values)) or 'none'
            sys.exit(f'{bam_path} must have exactly one distinct {what}; observed {observed}')
        return next(iter(values))


    def read_required_bam_sample_and_platform_unit(bam_path):
        with BamReader(bam_path) as reader:
            read_groups = list(reader.readGroupTable)
        if not read_groups:
            sys.exit(f'{bam_path} has no @RG records')

        sample = exactly_one(
            {str(rg.SampleName).strip() for rg in read_groups if str(rg.SampleName).strip()},
            '@RG SM sample',
            bam_path,
        )
        platform_unit = exactly_one(
            {str(rg.MovieName).strip() for rg in read_groups if str(rg.MovieName).strip()},
            '@RG PU acquisition/platform-unit',
            bam_path,
        )
        if platform_unit not in collection_contexts:
            sys.exit(
                f'{bam_path} @RG PU acquisition/platform-unit {platform_unit!r} does not match '
                'any validated consensusreadset_xmls CollectionMetadata Context; observed contexts: '
                + ', '.join(repr(value) for value in sorted(collection_contexts))
            )
        return sample, platform_unit


    def read_biosample_csv(path):
        expected_header = ('HiFi Barcode', 'cDNA Barcode', 'Bio Sample')
        cdna_barcodes_by_sample = {}
        with open(path, 'r', newline='', encoding='utf-8-sig') as fh:
            reader = csv.reader(fh)
            try:
                header = tuple(cell.strip() for cell in next(reader))
            except StopIteration:
                sys.exit(f'{path} is empty')
            if header != expected_header:
                sys.exit(f'{path} has header {header!r}; expected {expected_header!r}')

            for line_no, raw_row in enumerate(reader, start=2):
                row = [cell.strip() for cell in raw_row]
                if not row or all(cell == '' for cell in row):
                    continue
                if len(row) != 3:
                    sys.exit(f'{path}:{line_no}: expected 3 columns, got {len(row)}')
                _outer_barcode, cdna_barcode, sample = row
                if not cdna_barcode or not sample:
                    sys.exit(f'{path}:{line_no}: empty cDNA barcode or Bio Sample field')
                cdna_barcodes_by_sample.setdefault(sample, [])
                if cdna_barcode not in cdna_barcodes_by_sample[sample]:
                    cdna_barcodes_by_sample[sample].append(cdna_barcode)

        if not cdna_barcodes_by_sample:
            sys.exit(f'{path} has no data rows')
        return cdna_barcodes_by_sample


    def iter_biosamples(dataset):
        for bio_sample in dataset.metadata.bioSamples:
            yield bio_sample
        for collection in dataset.metadata.collections:
            for bio_sample in collection.wellSample.bioSamples:
                yield bio_sample


    def inherited_barcode_uuids(dataset):
        uuids_by_barcode = {}
        for bio_sample in iter_biosamples(dataset):
            for dna_barcode in bio_sample.DNABarcodes:
                barcode_name = (dna_barcode.name or '').strip()
                barcode_uuid = (dna_barcode.uniqueId or '').strip()
                if barcode_name and barcode_uuid:
                    uuids_by_barcode.setdefault(barcode_name, set()).add(barcode_uuid)
        return {
            barcode_name: next(iter(barcode_uuids))
            for barcode_name, barcode_uuids in uuids_by_barcode.items()
            if len(barcode_uuids) == 1
        }


    def required_consensus_read_set_uuid(dataset, xml_path):
        consensus_read_set_uuids = set()
        for index, collection in enumerate(dataset.metadata.collections, start=1):
            consensus_read_set_uuid = (collection.consensusReadSetRef.uuid or '').strip()
            if not consensus_read_set_uuid:
                sys.exit(f'{xml_path} CollectionMetadata[{index}] is missing ConsensusReadSetRef UniqueId')
            consensus_read_set_uuids.add(consensus_read_set_uuid)
        if len(consensus_read_set_uuids) != 1:
            sys.exit(
                f'{xml_path} must contain exactly one ConsensusReadSetRef UniqueId; observed '
                + ', '.join(repr(value) for value in sorted(consensus_read_set_uuids))
            )
        return next(iter(consensus_read_set_uuids))


    def link_or_copy(src, dest):
        try:
            os.link(src, dest)
        except OSError as exc:
            if exc.errno in (errno.EXDEV, errno.EEXIST, errno.EPERM):
                shutil.copy2(src, dest)
            else:
                raise


    def candidate_bam_name(path, index, attempt):
        basename = os.path.basename(path)
        if basename.endswith('.bam'):
            stem = basename.removesuffix('.bam')
            suffix = '.bam'
        else:
            stem, suffix = os.path.splitext(basename)
        marker = f'{index:03d}' if attempt == 1 else f'{index:03d}-{attempt}'
        return f'{stem}.{marker}{suffix}'


    def make_packaged_bam_name(path, index, colliding_basenames, used_names):
        basename = os.path.basename(path)
        if basename not in colliding_basenames and basename not in used_names:
            return basename
        attempt = 1
        while True:
            candidate = candidate_bam_name(path, index, attempt)
            if candidate not in used_names:
                return candidate
            attempt += 1


    def make_biosamples(samples, cdna_barcodes_by_sample, uuid_by_barcode):
        bio_samples = BioSamplesMetadata()
        for sample in samples:
            bio_sample = BioSampleMetadata()
            bio_sample.name = sample
            dna_barcodes = DNABarcodes()
            for cdna_barcode in cdna_barcodes_by_sample[sample]:
                dna_barcode = DNABarcode()
                dna_barcode.name = cdna_barcode
                barcode_uuid = uuid_by_barcode.get(cdna_barcode)
                if barcode_uuid:
                    dna_barcode.uniqueId = barcode_uuid
                dna_barcodes.append(dna_barcode)
            bio_sample.append(dna_barcodes)
            bio_samples.append(bio_sample)
        return bio_samples


    def replace_biosample_metadata(
        dataset,
        samples,
        cdna_barcodes_by_sample,
        uuid_by_barcode,
        consensus_read_set_uuid,
    ):
        # The consensus_read_set_uuid is the instrument ConsensusReadSetRef UUID and
        # must stay distinct from the generated dataset.uuid
        dataset.metadata.removeChildren('BioSamples')
        dataset.metadata.append(make_biosamples(samples, cdna_barcodes_by_sample, uuid_by_barcode))
        for collection in dataset.metadata.collections:
            collection.consensusReadSetRef.uuid = consensus_read_set_uuid
            well_sample = collection.wellSample
            well_sample.removeChildren('BioSamples')
            well_sample.append(make_biosamples(samples, cdna_barcodes_by_sample, uuid_by_barcode))


    def sanitize_xml_name(value):
        cleaned = re.sub(r'[^A-Za-z0-9._-]+', '_', value.strip())
        cleaned = cleaned.strip('._-')
        return cleaned or 'biosample'


    def make_child_xml_name(sample, sample_index):
        # sample_index is unique per produced biosample, so this name is unique.
        return f'{sample_index:03d}.{sanitize_xml_name(sample)}.flnc.consensusreadset.xml'


    def add_child_dataset_resource(parent_dataset, parent_bam, child_dataset, child_xml):
        parent_bam_keys = {
            parent_bam,
            os.path.abspath(parent_bam),
        }
        parent_resource = None
        for ext_resource in parent_dataset.externalResources:
            resource_keys = {
                ext_resource.resourceId,
                os.path.abspath(ext_resource.resourceId),
                ext_resource.attrib.get('ResourceId', ''),
            }
            if parent_bam_keys & resource_keys:
                parent_resource = ext_resource
                break
        if parent_resource is None:
            sys.exit(f'could not find parent FLNC BAM ExternalResource {parent_bam!r} while linking child dataset XML')

        child_resource = ExternalResource()
        child_resource.resourceId = child_xml
        child_resource.uniqueId = child_dataset.uuid
        child_resource.metaType = child_dataset.datasetType
        parent_resource._setSubResByMetaType(child_dataset.datasetType, child_resource)


    basename_counts = {}
    for bam in flnc_bams:
        basename_counts[os.path.basename(bam)] = basename_counts.get(os.path.basename(bam), 0) + 1
    colliding_basenames = {basename for basename, count in basename_counts.items() if count > 1}

    packaged_bams = []
    sample_order = []
    bams_by_sample = {}
    used_packaged_names = set()
    for index, (bam, pbi) in enumerate(zip(flnc_bams, flnc_bam_pbis, strict=True), start=1):
        packaged_bam_name = make_packaged_bam_name(
            bam,
            index,
            colliding_basenames,
            used_packaged_names,
        )
        used_packaged_names.add(packaged_bam_name)

        packaged_bam = packaged_bam_name
        packaged_pbi = packaged_bam + '.pbi'
        link_or_copy(bam, packaged_bam)
        link_or_copy(pbi, packaged_pbi)

        sample, _platform_unit = read_required_bam_sample_and_platform_unit(packaged_bam)

        packaged_bams.append(packaged_bam)
        if sample not in bams_by_sample:
            sample_order.append(sample)
            bams_by_sample[sample] = []
        bams_by_sample[sample].append(packaged_bam)

    ds_template = ConsensusReadSet(consensusreadset_xml, skipCounts=True, skipMissing=True)
    cdna_barcodes_by_sample = read_biosample_csv(validated_biosample_csv)
    missing_samples = [sample for sample in sample_order if sample not in cdna_barcodes_by_sample]
    if missing_samples:
        sys.exit(
            'final FLNC BAM @RG SM sample(s) are absent from the preprocessing biosample CSV: '
            + ', '.join(repr(sample) for sample in missing_samples)
        )
    uuid_by_barcode = inherited_barcode_uuids(ds_template)
    instrument_consensus_read_set_uuid = required_consensus_read_set_uuid(
        ds_template,
        consensusreadset_xml,
    )

    template_tags = [tag.strip() for tag in str(ds_template.tags or '').split(',') if tag.strip()]
    if 'flnc' not in {tag.lower() for tag in template_tags}:
        template_tags.append('flnc')
    flnc_tags = ','.join(template_tags)


    def build_flnc_dataset(dataset, samples, name):
        dataset.updateCounts()
        dataset.metadata.collections = copy.deepcopy(ds_template.metadata.collections)
        dataset.name = name
        dataset.tags = flnc_tags
        dataset.newUuid(random=True)
        replace_biosample_metadata(
            dataset,
            samples,
            cdna_barcodes_by_sample,
            uuid_by_barcode,
            consensus_read_set_uuid=instrument_consensus_read_set_uuid,
        )


    ds_out = ConsensusReadSet(*packaged_bams)
    build_flnc_dataset(ds_out, sample_order, ds_template.name)
    top_level_uuid = ds_out.uuid

    child_xmls = []
    for sample_index, sample in enumerate(sample_order, start=1):
        child_ds = ConsensusReadSet(*bams_by_sample[sample])
        build_flnc_dataset(child_ds, [sample], f'{ds_template.name} {sample} FLNC')
        child_ds.metadata.addParentDataSet(
            top_level_uuid,
            ds_out.datasetType,
            ds_out.timeStampedName or ds_out.name or dataset_xml,
        )
        child_xml = make_child_xml_name(sample, sample_index)
        child_xmls.append(child_xml)
        child_ds.write(child_xml, relPaths=True)
        add_child_dataset_resource(ds_out, bams_by_sample[sample][0], child_ds, child_xml)

    ds_out.write(dataset_xml, relPaths=True)


    def write_manifest(path, values):
        with open(path, 'w', encoding='utf-8') as out_fh:
            out_fh.writelines(value + '\n' for value in values)


    write_manifest('child_flnc_dataset_xmls.txt', child_xmls)
    write_manifest('packaged_flnc_bams.txt', packaged_bams)
    write_manifest('packaged_flnc_bam_pbis.txt', [bam + '.pbi' for bam in packaged_bams])
    PY
  >>>

  output {
    File flnc_dataset_xml = "flnc.consensusreadset.xml"
    Array[File] child_flnc_dataset_xmls = read_lines("child_flnc_dataset_xmls.txt")
    Array[File] packaged_flnc_bams = read_lines("packaged_flnc_bams.txt")
    Array[File] packaged_flnc_bam_pbis = read_lines("packaged_flnc_bam_pbis.txt")
  }

  runtime {
    cpu: threads
    memory: total_mem_mb + " MB"
    docker: runtime_attributes.container_registry + "/pbcore@sha256:47fc6f1174605be9a8a932f17c15d1a29fc5e6b89ea6809f32995318c7afe1e7"  # 2.6.0_build3
    maxRetries: runtime_attributes.max_retries
  }
}
