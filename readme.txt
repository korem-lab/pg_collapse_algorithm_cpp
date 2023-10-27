COMMANDS
conda activate copan
cd Documents/GitHub/pg_collapse_algorithm/src
python setup.py build_ext --inplace
python main.py --contig-list ~/Documents/GitHub/data/sample_5_megahit_20M.pe_ext_all.fasta --graph-name newtest --out-dir ~/Documents/GitHub/data/output --div-thresh 0.02 --window-size 10 -t 12

OUTPUT
Copangraph run configurations:
contig_list: /home/nat/Documents/GitHub/data/sample_5_megahit_20M.pe_ext_all.fasta
graph_name: newtest_possible_fix
out_dir: /home/nat/Documents/GitHub/data/output
div_thresh: 0.02
window_size: 10
threads: 16
asymmetrical: False
sensitive_mode: False
keep_short_extensions: False
kmer_size: 15
max_separation: 75
min_overlap: 100
max_jump: 200
large_gap_cost: 2.0
small_gap_cost: 0.5
high_freq_kmer_filter: 1e-05

aligning contigs for newtest_possible_fix...
Sketching 84461 contigs
window_size: 10
counting kmers...
total kmers: 107857664
total avoided high-frequency kmers 0
Sketching kmers
Transferring 0.0% of sketches
Transferring 10.001065580563811% of sketches
Transferring 20.002131161127622% of sketches
Transferring 30.00319674169143% of sketches
Transferring 40.004262322255244% of sketches
Transferring 50.005327902819054% of sketches
Transferring 60.00639348338286% of sketches
Transferring 70.00745906394667% of sketches
Transferring 80.00852464451049% of sketches
Transferring 90.0095902250743% of sketches
Building kmer map
Indexing 84461 contigs.
Built index.
23638535
num_threads 16
thread 0 analyzed contig 0
thread 4 analyzed contig 5000
thread 1 analyzed contig 10000
thread 9 analyzed contig 15000
thread 1 analyzed contig 20000
thread 4 analyzed contig 25000
thread 2 analyzed contig 30000
thread 14 analyzed contig 35000
thread 7 analyzed contig 40000
thread 10 analyzed contig 45000
thread 15 analyzed contig 50000
thread 0 analyzed contig 55000
thread 3 analyzed contig 60000
thread 11 analyzed contig 65000
thread 15 analyzed contig 70000
thread 11 analyzed contig 75000
thread 2 analyzed contig 80000
Merging per-thread alignments
returning from align
num py_alignments 3465848
generating gluepoints for newtest_possible_fix
construct endpoints from alignments
cluster endpoints
initialize
group breakpoints into gluepoints
build consensus gluepoints
building repeat graph for newtest_possible_fix
writing data...
