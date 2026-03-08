"""BLAST feature generation (currently optional in scoring pipeline)."""

import os
import pickle

import pandas as pd
from Bio.Blast.Applications import NcbiblastnCommandline

from .Variables import DB_HOST_PREFIX, HASH_TABLE


db_host_prefix = os.path.expanduser(DB_HOST_PREFIX)
dict_genome = pickle.load(open(HASH_TABLE, 'rb'))


def cal_perc(x):
    indicator = [0] * x[4]
    for i in range(len(x[2])):
        indicator[x[2][i]:x[3][i]] = [1] * (x[3][i] - x[2][i] + 1)
    return sum(indicator)


def blastSingle(item, query_virus_dir, output_dir, seqid, numThreads):
    query_name = item.split('.')[0]
    query_file = os.path.join(query_virus_dir, item)
    output_file = os.path.join(output_dir, query_name) + '.blast'

    blast_kwargs = dict(
        query=query_file,
        db=db_host_prefix,
        out=output_file,
        outfmt="6 qacc sacc qstart qend qlen",
        evalue=0.01,
        gapopen=10,
        penalty=-1,
        reward=1,
        gapextend=2,
        word_size=11,
        perc_identity=90,
        num_threads=numThreads,
    )
    if seqid is not None:
        blast_kwargs['seqidlist'] = seqid

    blast_call = NcbiblastnCommandline(**blast_kwargs)
    blast_call()

    if os.stat(output_file).st_size == 0:
        return False, None

    with open(query_file) as f:
        query_len = len(f.read())

    query_res = pd.read_table(output_file, header=None)
    query_res[1] = [dict_genome[k] for k in list(query_res[1])]
    df_blast_positions = query_res.groupby([0, 1]).agg(
        {
            2: lambda x: tuple(x - 1),
            3: lambda x: tuple(x - 1),
            4: min,
        }
    )
    df_blast_positions.index = df_blast_positions.index.droplevel()
    df_blast_perc = df_blast_positions.apply(lambda x: cal_perc(x), axis=1) / query_len
    sr_blast = df_blast_perc.groupby(level=0, sort=False).apply(sum)
    return True, pd.DataFrame({query_name: sr_blast}).T


def blast_calculator(query_virus_dir, virus_index, host_index, genomes, output_dir, numThreads=1):
    blast_output_dir = os.path.join(output_dir, 'BLAST/')
    os.makedirs(blast_output_dir, exist_ok=True)

    query_cont = []
    query_list = [
        f for f in os.listdir(query_virus_dir)
        if f.lower().endswith(('.fasta', '.fa', '.fna'))
    ]

    if genomes is not None:
        genome_set = set(genomes)
        seqid_list = [i for i in dict_genome if dict_genome[i] in genome_set]
        seqid = os.path.join(output_dir, 'seqid_list.txt')
        pd.Series(seqid_list).to_csv(seqid, index=None)
    else:
        seqid = None

    for item in query_list:
        print('----Calculating blast feature values for ', item, ' ----')
        ind, df = blastSingle(item, query_virus_dir, blast_output_dir, seqid, numThreads)
        if ind:
            query_cont.append(df)

    df_pseudo = pd.DataFrame(index=virus_index, columns=host_index).fillna(0)
    print('----BLAST intermediate files are stored in ', blast_output_dir, ' ----')
    if query_cont == []:
        return df_pseudo

    merged = pd.concat(query_cont, axis=1, sort=False)
    tmp = merged.T.groupby(level=0, sort=False).sum().T
    return pd.concat([tmp, df_pseudo], sort=False).groupby(level=0, sort=False).sum().fillna(0).loc[virus_index][host_index]
