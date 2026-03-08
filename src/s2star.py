"""s2* feature generation."""

import os
import pickle
import sys

import numpy as np
import pandas as pd

from .tools import kmer_count
from .Variables import TABLE_BENCH, TABLE_HOST, TABLE_INTER


ALPHABET = ['A', 'C', 'G', 'T']
ALPHA_DICT = dict(zip(ALPHABET, range(4)))


def _kmer_count_fallback(seqfile, K, reverse):
    """Fallback k-mer counter when native extension output is unusable."""
    counts = np.zeros(4 ** K, dtype=np.int64)
    mask = (1 << (2 * (K - 1))) - 1

    def _revcomp_index(num):
        rc = 0
        for i in range(K):
            shift = 2 * (K - i - 1)
            rc += (3 - ((num >> shift) & 3)) * (4 ** i)
        return rc

    with open(seqfile) as f:
        num = 0
        j = 0
        for line in f:
            if not line or line[0] == '>':
                num = 0
                j = 0
                continue
            for nuc in line.strip().upper():
                code = ALPHA_DICT.get(nuc)
                if code is None:
                    num = 0
                    j = 0
                    continue
                if j < (K - 1):
                    num = num * 4 + code
                    j += 1
                else:
                    num = ((num & mask) << 2) + code
                    counts[num] += 1
                    if reverse:
                        counts[_revcomp_index(num)] += 1
    return counts


def _is_nonempty_fasta(path):
    with open(path) as f:
        for line in f:
            if line and line[0] != '>' and line.strip():
                return True
    return False


def _native_kmer_supported(probe_file, numThreads, k):
    counts = np.array(kmer_count(probe_file, numThreads, True, k))
    if counts.size == 0 or counts[0] == -1:
        return False
    return np.sum(counts) > 0


def get_transition(count_array):
    transition_array = count_array.reshape(len(count_array) // 4, 4)
    with np.errstate(divide='ignore', invalid='ignore'):
        transition_array = transition_array / np.sum(transition_array, 1)[:, np.newaxis]
        transition_array[np.isnan(transition_array)] = 0
    return transition_array


def get_expect(a_M_count, a_trans, K, M):
    a_expect = a_M_count
    for _ in range(K - M):
        a_trans = np.tile(a_trans, (4, 1))
        a_expect = (a_expect[:, np.newaxis] * a_trans).flatten()
    return a_expect


def get_f(a_K_count, a_expect):
    with np.errstate(divide='ignore', invalid='ignore'):
        f = (a_K_count - a_expect) / np.sqrt(a_expect)
        f[np.isnan(f)] = 0
    return f


def get_all_f(directory, K, order, reverse, numThreads):
    M = order + 1
    sequence_list = [
        f for f in os.listdir(directory)
        if f.lower().endswith(('.fasta', '.fa', '.fna'))
    ]
    f_matrix = np.ones((len(sequence_list), 4 ** K))
    has_error = False

    probe_file = None
    for seq in sequence_list:
        seqfile = os.path.join(directory, seq)
        if os.path.getsize(seqfile) > 0 and _is_nonempty_fasta(seqfile):
            probe_file = seqfile
            break

    use_native_k = probe_file is not None and _native_kmer_supported(probe_file, numThreads, K)
    use_native_m = probe_file is not None and _native_kmer_supported(probe_file, numThreads, M)

    if probe_file is not None and (not use_native_k or not use_native_m):
        print('----kmer_count native path unavailable; using Python fallback in s2*----')

    for i, seq in enumerate(sequence_list):
        seqfile = os.path.join(directory, seq)

        if use_native_k:
            k_count = np.array(kmer_count(seqfile, numThreads, reverse, K))
        else:
            k_count = _kmer_count_fallback(seqfile, K, reverse)
        if np.sum(k_count) == 0 and os.path.getsize(seqfile) > 0:
            k_count = _kmer_count_fallback(seqfile, K, reverse)

        if k_count[0] == -1:
            print('The query file {} contains invalid chars, please make sure it is a valid fasta file.'.format(seq))
            has_error = True
        if np.sum(k_count) == 0:
            print('The query file {} is empty, please double check the file.'.format(seq))
            has_error = True

        if use_native_m:
            m_count = np.array(kmer_count(seqfile, numThreads, reverse, M))
        else:
            m_count = _kmer_count_fallback(seqfile, M, reverse)
        if np.sum(m_count) == 0 and os.path.getsize(seqfile) > 0:
            m_count = _kmer_count_fallback(seqfile, M, reverse)

        trans = get_transition(m_count)
        expect = get_expect(m_count, trans, K, M)
        f_matrix[i] = get_f(k_count, expect)

    if has_error:
        sys.exit('Program terminated. Please check error info above.')

    name_list = [x.rsplit('.', 1)[0] for x in sequence_list]
    return f_matrix, name_list


def cosine_similarity(f1, f2):
    n1 = np.linalg.norm(f1, axis=1, keepdims=True)
    n2 = np.linalg.norm(f2, axis=1, keepdims=True)
    prod = np.dot(f1, f2.T)
    return prod / (np.dot(n1, n2.T))


def s2_query(query_virus_dir, ifShort, numThreads):
    print('----Start calculating s2* part I... ----')
    virus_mat, virus_list = get_all_f(query_virus_dir, 6, 2, True, numThreads)
    print('----Finished calculating s2* part I----')

    print('----Start calculating s2* part II... ----')
    with open(TABLE_BENCH, 'rb') as f:
        bench_mat, bench_list = pickle.load(f)
    s2_virus_bench_mat = cosine_similarity(virus_mat, bench_mat)
    print('----Finished calculating s2* part II----')

    if ifShort:
        return pd.DataFrame(s2_virus_bench_mat, index=virus_list, columns=bench_list), None

    print('----Start calculating s2* part III... ----')
    with open(TABLE_HOST, 'rb') as f:
        host_mat, host_list = pickle.load(f)
    s2_host_virus_mat = cosine_similarity(virus_mat, host_mat)
    print('----Finished calculating s2* part III----')
    return (
        pd.DataFrame(s2_virus_bench_mat, index=virus_list, columns=bench_list),
        pd.DataFrame(s2_host_virus_mat, index=virus_list, columns=host_list),
    )


def s2star_calculator(query_virus_dir, ifShort, numThreads):
    mat_original_interaction = pd.read_csv(TABLE_INTER)
    if ifShort:
        print('----Calculation of s2* is split into two parts----')
    else:
        print('----Calculation of s2* is split into three parts----')

    mat_query_bench, s2star_query_host = s2_query(query_virus_dir, ifShort, numThreads)
    host_index = mat_original_interaction.columns
    query_index = mat_query_bench.index
    if not ifShort:
        s2star_query_host = s2star_query_host.loc[query_index, host_index]
    return s2star_query_host, mat_query_bench, mat_original_interaction


# Backward-compatible alias (legacy typo kept to avoid breakage).
def s2star_caclculator(query_virus_dir, ifShort, numThreads):
    return s2star_calculator(query_virus_dir, ifShort, numThreads)
