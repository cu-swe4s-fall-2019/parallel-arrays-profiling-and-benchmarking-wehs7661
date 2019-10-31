import data_viz
import gzip
import sys
import os
import argparse
import time

def initialize():
    """
    An initializing function
    """

    parser = argparse.ArgumentParser(
        description='This code plots the gene expression distribution across either tissue groups (SMTS), or tissue types (SMTSD) for a target gene')
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

    args_parse = parser.parse_args()

    return args_parse


def linear_search(key, L):
    """ Returns index with matched key using a linear searching method
    """
    for i in range(len(L)):
        if key == L[i]:
            return i
    return -1


def linear_search_all_hits(key, L):
    """ Gives indices not values
    """
    hit = []
    for i in range(len(L)):
        if key == L[i]:
            hit.append(i)
    return hit


def binary_search(key, D):
    """ Returns index with matched key using a binary searching method
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
    samples, target_group = [], []

    for l in open(args.sample_attributes):
        sample_info = l.rstrip().split('\t')

        if metadata_header is None:
            metadata_header = sample_info
            continue

        sample_idx = linear_search('SAMPID', metadata_header)
        target_idx = linear_search(args.group_type, metadata_header)
        if (target_idx == -1):
            return [], []
        samples.append(sample_info[sample_idx])
        target_group.append(sample_info[target_idx])

    if not target_group:
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
            
            # Benchmarking
            # start = time.time()
            rna_header_plus_index.sort()
            # end = time.time()
            # print(end - start)
            continue

        rna_counts = l.rstrip().split('\t')
        description_idx = linear_search('Description', rna_header)

        if description_idx == -1:
            print('No genes found in the header')
            sys.exit(1)

        if rna_counts[description_idx] == target_gene_name:
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
            data_viz.boxplot(par_array, attrs, args.group_type,
                             'Gene read counts', target_gene_name,
                             args.output_file)
            # search_end = time.time()
            # print(search_end - search_start)
            sys.exit(0)

    sys.exit(0)

if __name__ == '__main__':
    main()
