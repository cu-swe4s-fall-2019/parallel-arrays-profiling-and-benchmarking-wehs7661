import data_viz
import gzip
import sys
import os
import argparse
import time
sys.path.insert(1, 'hash-tables-wehs7661')
import hash_functions
import hash_tables


def initialize():
    """
    An initializing function

    Parameters
    ----------
    No parameters required

    Return
    ------
    args_parse : namespace
        a namespace containing the argument parser
    """

    parser = argparse.ArgumentParser(
        description='This code plots the gene expression distribution across \
                    either tissue groups (SMTS), or tissue types (SMTSD) for \
                    a target gene')
    parser.add_argument('--output_file',
                        type=str,
                        help='The file name of the output file (figure)',
                        required=True)

    parser.add_argument('--group_type',
                        type=str,
                        help='The type of the group',
                        required=True)

    parser.add_argument('--gene',
                        type=str,
                        help='The name of the target gene',
                        required=True)

    parser.add_argument('--sample_attributes',
                        type=str,
                        help='An input file which contains metadata',
                        required=True)

    parser.add_argument('--gene_reads',
                        type=str,
                        help='The database of genes',
                        required=True)
    parser.add_argument('--data_structure',
                        type=str,
                        help='The type of data structure. Available options are \
                            "parallel" (parallel array) and "hash" (hash table)',
                        required=True)

    args_parse = parser.parse_args()

    return args_parse


def linear_search(key, L):
    """
    This function returns index with matched key using a linear searching
    method

    Parameters
    ----------
    key : str or int or float
        the key for finding the associated value
    L : list
        the dataset in which the key to be searched for

    Returns
    -------
    i : int
        the index of the given key (If no key is found, the function
        returns -1.)
    """
    for i in range(len(L)):
        if key == L[i]:
            return i
    return -1


def linear_search_all_hits(key, L):
    """
    This function returns a list of index with matched key using a linear
    searching method

    Parameters
    ----------
    key : str or int or float
        the key for finding the associated values
    L : list
        the dataset in which the key to be searched for

    Returns
    -------
    hit : list
        a list of index of the given key
    """
    hit = []
    for i in range(len(L)):
        if key == L[i]:
            hit.append(i)
    return hit


def binary_search(key, D):
    """
    This function returns index with matched key using a binary searching
    method

    Parameters
    ----------
    key : str or int or float
        the key for finding the associated value
    D : list
        the dataset in which the key to be searched for

    Returns
    -------
    D[mid][1] : int
        the index of the given key (If no key is found, the function
        returns -1.)
    """
    lo = -1
    hi = len(D)
    while (hi - lo > 1):
        mid = (hi + lo) // 2

        if key == D[mid][0]:
            return D[mid][1]

        if (key < D[mid][0]):
            hi = mid
        else:
            lo = mid

    return -1


def main():
    args = initialize()

    # check if the input files exist
    if (not os.path.exists(args.sample_attributes)):
        print('Metadata file not found')
        sys.exit(1)
    if (not os.path.exists(args.gene_reads)):
        print('Gene data file not found')
        sys.exit(1)

    target_gene_name = args.gene
    metadata_header = None

    if args.data_structure == 'parallel':
        samples, target_group = [], []   # only for parallel array
    elif args.data_structure == 'hash':
        target_group = []
        ht_meta = hash_tables.ChainedHash(100000, hash_functions.h_rolling)
    else:
        print('Please input data structures available.')
        print('Options available include "parallel" and "hash".')
        sys.exit(1)

    for l in open(args.sample_attributes):
        sample_info = l.rstrip().split('\t')

        if metadata_header is None:
            metadata_header = sample_info
            continue

        sample_idx = linear_search('SAMPID', metadata_header)
        target_idx = linear_search(args.group_type, metadata_header)
        if (target_idx == -1):
            break   # no such group

        if args.data_structure == 'parallel':
            samples.append(sample_info[sample_idx])       # ID
            target_group.append(sample_info[target_idx])  # group type
        elif args.data_structure == 'hash':
            key = sample_info[target_idx]                 # group type
            value = sample_info[sample_idx]               # ID 
            search = ht_meta.search(key)
            if search is None:
                ht_meta.add(key, [value]) # map ID and group
                target_group.append(key)
            else:
                search.append(value)

    if len(target_group) == 0:
        print('Group type not found')
        sys.exit(1)

    version, dim, rna_header = None, None, None

    for l in gzip.open(args.gene_reads, 'rt'):

        if version is None:
            version = l
            continue

        if dim is None:
            dim = l
            continue

        if rna_header is None:
            rna_header = l.rstrip().split('\t')
            rna_header_plus_index = []
            for i in range(len(rna_header)):
                rna_header_plus_index.append([rna_header[i], i])
            rna_header_plus_index.sort()
            continue

        rna_counts = l.rstrip().split('\t')
        description_idx = linear_search('Description', rna_header)

        if description_idx == -1:
            print('No genes found in the header')
            sys.exit(1)

        if rna_counts[description_idx] == target_gene_name:
            if args.data_structure == 'parallel':
                attrs = list(set(target_group))
                attrs.sort()
                par_array = []
                # search_start = time.time()
                for attr in attrs:
                    attr_idxs = linear_search_all_hits(attr, target_group)

                    attr_counts = []
                    for attr_idx in attr_idxs:
                        rna_header_idx = linear_search(samples[attr_idx],
                                                    rna_header)
                        # rna_header_idx = binary_search(samples[attr_idx],
                        #                                rna_header_plus_index)
                        if rna_header_idx == -1:
                            continue
                        count = rna_counts[rna_header_idx]
                        attr_counts.append(int(count))
                    par_array.append(attr_counts)
                data_viz.boxplot(par_array, target_group, args.group_type,
                                'Gene read counts', target_gene_name,
                                args.output_file)
                # search_end = time.time()
                # print(search_end - search_start)
                sys.exit(0)
            
            elif args.data_structure == 'hash':
                counts_list = []
                ht_rna = hash_tables.ChainedHash(100000, hash_functions.h_rolling)
                for i in range(description_idx + 1, len(rna_header)):
                    ht_rna.add(rna_header[i], int(rna_counts[i]))  # map ID and counts
                target_group.sort()
                for attr in target_group:
                    attr_counts = []
                    sampID = ht_meta.search(attr)
                    if sampID is None:
                        continue
                    for ID in sampID:
                        count = ht_rna.search(ID)
                        if count is None:
                            continue
                        attr_counts.append(count)
                    counts_list.append(attr_counts)
                data_viz.boxplot(counts_list, target_group, args.group_type,
                                'Gene read counts', target_gene_name,
                                args.output_file)
                sys.exit(0)
    sys.exit(0)

if __name__ == '__main__':
    main()
