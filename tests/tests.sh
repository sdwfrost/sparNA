# rm -rf output

source ../../venv/bin/activate

# Fast agnostic HCV multisample
# ../pipeline.py --fwd-reads-sig _F --rev-reads-sig _R --norm-k-list 25,31 --norm-cov-list 5,10,15,20 --asm-k-list 33,43 --reads-dir data --threads 12 --reference /Users/Bede/Research/Datasets/ref_sequences/hcv/h77.fasta

# Thorough agnostic HCV multisample
../pipeline.py --fwd-reads-sig _F --rev-reads-sig _R --norm-k-list 31 --norm-cov-list 1,2,5,10,20 --asm-k-list 21,33,55,77 --reads-dir data --threads 12 --reference /Users/Bede/Research/Datasets/ref_sequences/hcv/h77.fasta

sleep 1
# ../pipeline.py --fwd-reads /Users/Bede/Research/Analyses/phe_asm/phe_hcv_pipeline/tests/data/hcv_n1000_F.fastq --rev-reads /Users/Bede/Research/Analyses/phe_asm/phe_hcv_pipeline/tests/data/hcv_n1000_R.fastq --fwd-reads-sig _F --rev-reads-sig _R --norm-k-list 25,31 --norm-c-list 5,10,15,20 --asm-k-list 33,43  --threads 12

# --reference --asm-using-ref

curl bede.im/notify.php?assembly_finished

tput bel