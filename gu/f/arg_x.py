"""Non-PAR X input handling and explicit biological haplotype identities."""
from __future__ import annotations



def nonpar(position, build):
    # One-based inclusive non-PAR interval between PAR1 and PAR2.
    start, end = {'37': (2699521, 154931043), '38': (2781480, 155701382)}[str(build)]
    return start <= position <= end


def read_sexes(path):
    with open(path) as handle:
        header = handle.readline().lstrip('#').split()
        if 'IID' not in header or 'SEX' not in header:
            raise ValueError(f'X requires IID and SEX in PSAM: {path}')
        result = {}
        for line in handle:
            if not line.strip():
                continue
            row = dict(zip(header, line.split()))
            if row['IID'] in result or row['SEX'] not in ('1', '2'):
                raise ValueError(f'Duplicate IID or unknown X sex: {row["IID"]}')
            result[row['IID']] = int(row['SEX'])
    return result


def haploid_zarr(source, destination):
    """Flatten only real GT copies, never infer the -2 padding as a sample.

    Uniform haploid storage also works with tsinfer's partitioned matcher.
    Scan every GT chunk to enforce constant per-individual ploidy.
    """
    import numpy as np
    import numcodecs
    import zarr

    src = zarr.open(str(source), mode='r')
    gt = src['call_genotype']
    names = [str(x) for x in src['sample_id'][:]]
    present = np.asarray(gt[0]) != -2
    ploidies = present.sum(axis=1)
    if np.any((ploidies < 1) | (ploidies > 2)) or np.any(~present[:, 0]):
        raise ValueError('Invalid per-sample X ploidy in VCF-Zarr')
    identities = [(name, h + 1) for name, p in zip(names, ploidies) for h in range(int(p))]
    dest = zarr.open_group(str(destination), mode='w')
    for key in ('variant_position', 'variant_allele', 'variant_AA', 'variant_contig', 'contig_id', 'contig_length'):
        if key in src:
            zarr.copy(src[key], dest, name=key)
    ids = np.asarray([f'arg_hap_{i}' for i in range(len(identities))], dtype=object)
    dest.create_dataset('sample_id', data=ids, object_codec=numcodecs.VLenUTF8())
    width = len(identities)
    chunk = min(gt.chunks[0], max(1, 32_000_000 // (gt.shape[1] * gt.shape[2])))
    out = dest.create_dataset('call_genotype', shape=(gt.shape[0], width, 1),
                              chunks=(chunk, min(width, 256), 1), dtype=gt.dtype)
    for start in range(0, gt.shape[0], chunk):
        values = np.asarray(gt[start:start + chunk])
        if np.any((values != -2) != present[None, :, :]):
            raise ValueError(f'X ploidy changes within chromosome at variant chunk {start}')
        if np.any(values < -2) or np.any(values > 1):
            raise ValueError('Non-biallelic GT in prepared X VCF')
        out[start:start + len(values), :, 0] = values[:, present]
    zarr.consolidate_metadata(dest.store)
    return identities


def restore_individuals(ts, identities):
    """Replace storage-only haploid identities with original biological samples."""
    import numpy as np
    import tskit

    if ts.num_samples != len(identities):
        raise ValueError('X ARG sample count differs from real input haplotypes')
    tables = ts.dump_tables()
    individuals = {}
    tables.individuals.clear()
    tables.individuals.metadata_schema = tskit.MetadataSchema.permissive_json()
    node_individual = np.full(ts.num_nodes, tskit.NULL, dtype=np.int32)
    for node, (name, hap) in zip(ts.samples(), identities):
        if name not in individuals:
            individuals[name] = tables.individuals.add_row(metadata={'sample': name})
        node_individual[node] = individuals[name]
    tables.nodes.individual = node_individual
    return tables.tree_sequence()
