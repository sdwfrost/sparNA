# rm -rf output

source ../../venv/bin/activate

# Agnostic HCV multisample
../pipeline.py --fwd-reads-sig _F --rev-reads-sig _R --norm-k-list 25,31 --norm-cov-list 5,10,15,20 --asm-k-list 33,43  --reads-dir data --threads 12 --reference /Users/Bede/Research/Datasets/ref_sequences/hcv/h77.fasta
sleep 1
# ../pipeline.py --fwd-reads /Users/Bede/Research/Analyses/phe_asm/phe_hcv_pipeline/tests/data/hcv_n1000_F.fastq --rev-reads /Users/Bede/Research/Analyses/phe_asm/phe_hcv_pipeline/tests/data/hcv_n1000_R.fastq --fwd-reads-sig _F --rev-reads-sig _R --norm-k-list 25,31 --norm-c-list 5,10,15,20 --asm-k-list 33,43  --threads 12

# --reference --asm-using-ref

tput bel