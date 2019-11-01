#!/bin/bash
test -e ssshtest || wget https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest
. ssshtest

run test_style pycodestyle plot_gtex.py
assert_no_stdout
run test_style pycodestyle data_viz.py
assert_no_stdout
run test_style pycodestyle test_data_viz.py
assert_no_stdout
run test_style pycodestyle test_plot_gtex.py
assert_no_stdout

# test for invalid file
run invalidfile python plot_gtex.py --gene_reads AAA.txt --sample_attributes GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt --gene ACTA2 --group_type SMTS --output_file test.png --data_structure parallel
assert_exit_code 1
assert_stdout

# test for invalid attribute
run invalidfile python plot_gtex.py --gene_reads AAA.txt --sample_attributes AAA.txt --gene ACTA2 --group_type SMTS --output_file test.png --data_structure parallel
assert_exit_code 1
assert_stdout

# invalid gene
run invalid_gene python plot_gtex.py --gene_reads GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.acmg_59.gct.gz --sample_attributes GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt --gene AAA --group_type SMTS --output_file test.png --data_structure parallel
assert_exit_code 0
assert_no_stdout

# invalid group type
run invalid_group python plot_gtex.py --gene_reads GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.acmg_59.gct.gz --sample_attributes GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt --gene ACTA2 --group_type AAA --output_file test.png --data_structure parallel
assert_exit_code 1
assert_stdout

# invalid data structure
run box1 python plot_gtex.py --gene_reads GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.acmg_59.gct.gz --sample_attributes GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt --gene ACTA2 --group_type SMTS --output_file test.png --data_structure AAA
assert_exit_code 1
assert_stdout

# normal inputs (hash tables)
run box1 python plot_gtex.py --gene_reads GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.acmg_59.gct.gz --sample_attributes GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt --gene ACTA2 --group_type SMTS --output_file test.png --data_structure hash
assert_exit_code 0
assert_no_stdout
rm test.png
