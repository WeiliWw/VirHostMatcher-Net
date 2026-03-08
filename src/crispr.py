"""CRISPR feature generation utilities."""

import math
import os
from concurrent.futures import ThreadPoolExecutor

import pandas as pd
from Bio.Application import ApplicationError
from Bio.Blast.Applications import NcbiblastnCommandline
from numpy import isnan

from .Variables import DB_HOST_CRISPR_PREFIX, TAXA_INFO


db_host_crispr_prefix = os.path.expanduser(DB_HOST_CRISPR_PREFIX)
taxa_info = pd.read_pickle(TAXA_INFO)


def crisprSingle(item, query_virus_dir, output_dir, numThreads):
    """Run CRISPR BLAST for one query and return a sparse feature table."""
    query_name = item.split('.')[0]
    query_file = os.path.join(query_virus_dir, item)
    output_file = os.path.join(output_dir, query_name) + '.crispr'

    crispr_call = NcbiblastnCommandline(
        query=query_file,
        db=db_host_crispr_prefix,
        out=output_file,
        outfmt="6 qacc sacc evalue",
        evalue=1,
        gapopen=10,
        penalty=-1,
        gapextend=2,
        word_size=7,
        dust='no',
        task='blastn-short',
        perc_identity=90,
        num_threads=numThreads,
    )
    try:
        crispr_call()
    except ApplicationError as e:
        print('Warning: BLAST failed for {}: {}'.format(query_name, e))
        return False, None

    if os.stat(output_file).st_size == 0:
        return False, None

    query_res = pd.read_table(output_file, header=None)
    query_res = query_res[query_res[1].apply(lambda x: x.count("|")) == 2]
    if query_res.shape[0] == 0:
        return False, None

    query_res[0] = query_name
    query_res[1] = query_res[1].apply(lambda x: x.split('|')[-2])
    query_res[2] = -query_res[2].apply(math.log)
    df_crispr = query_res.groupby([0, 1]).max().unstack(fill_value=0)
    return True, df_crispr.set_index([[query_name]])


def crispr_calculator(query_virus_dir, output_dir, numThreads):
    """Calculate CRISPR features for all query fasta files."""
    query_cont = []
    query_list = [
        f for f in os.listdir(query_virus_dir)
        if f.lower().endswith(('.fasta', '.fa', '.fna'))
    ]

    crispr_output_dir = os.path.join(output_dir, 'CRISPR/')
    os.makedirs(crispr_output_dir, exist_ok=True)

    if query_list:
        total_threads = max(1, int(numThreads))
        workers = min(len(query_list), total_threads)
        threads_per_blast = max(1, total_threads // workers)

        if workers > 1:
            print(
                '----CRISPR parallel mode: {} workers x {} blast threads (budget {})----'.format(
                    workers, threads_per_blast, total_threads
                )
            )
            with ThreadPoolExecutor(max_workers=workers) as executor:
                futures = {
                    item: executor.submit(
                        crisprSingle,
                        item,
                        query_virus_dir,
                        crispr_output_dir,
                        threads_per_blast,
                    )
                    for item in query_list
                }
                for item in query_list:
                    print('----Calculating crispr feature values for ', item, ' ----')
                    ind, df = futures[item].result()
                    if ind:
                        query_cont.append(df)
        else:
            for item in query_list:
                print('----Calculating crispr feature values for ', item, ' ----')
                ind, df = crisprSingle(item, query_virus_dir, crispr_output_dir, threads_per_blast)
                if ind:
                    query_cont.append(df)

    print('----CRISPR intermediate files are stored in ', crispr_output_dir, ' ----')
    if query_cont == []:
        return pd.DataFrame()

    merged = pd.concat(query_cont, axis=1, sort=False)
    return merged.T.groupby(level=1, sort=False).sum().T.fillna(0)


def uniGenus(df_input, virus_index, host_index):
    """Expand CRISPR hits to genus level and align to full host index."""
    df_pseudo = pd.DataFrame(index=virus_index, columns=host_index).fillna(0)
    df_full = pd.concat([df_pseudo, df_input], sort=False).groupby(level=0, sort=False).sum().fillna(0)
    df_full = df_full.loc[virus_index][host_index]
    df_full.loc['hostGenus'] = taxa_info.loc[host_index]['hostGenus']

    dict_genera = {}
    for host_name in df_full:
        col = df_full[host_name].rename(None)
        genus = col['hostGenus']
        if type(genus) is not str and isnan(genus):
            continue
        if genus in dict_genera:
            dict_genera[genus] = pd.concat([dict_genera[genus], col], axis=1, sort=False).max(axis=1)
        else:
            dict_genera[genus] = col

    for host_name in df_full:
        genus = df_full[host_name]['hostGenus']
        if type(genus) is not str and isnan(genus):
            continue
        df_full[host_name] = dict_genera[genus]
    return df_full.loc[virus_index, :]
