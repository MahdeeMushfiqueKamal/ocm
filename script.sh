#!/bin/bash
python3 input/breaker.py input/lhg22L20MC5x.fa 5        #input fasta file number of thread
rm ocm
g++ main.cpp -ltbb -pthread -std=c++17 -o ocm
#g++ -std=c++17 -ltbb -lpthread -o ocm main.cpp
#./ocm count [-c = conservative update] [-r = do not canonicalize] -k <kmer-length> [-h <counters-height>] [-w <counters-width>] [-n <num-rounds>] [-t <num-threads>] -o <out-ocm-file> -fa <input-fasta-file>
./ocm test -k 22 -h 7 -w 1048576 -n 4 -t 5 -o output.bin -fa input/lhg22L20MC5x.fa -q input/test_exact_count_lhg22L20MC5x_20000.txt
#./ocm query -f output.bin -q input/test_exact_count_lhg22L20MC5x_20000.txt -o output/query_result.csv
