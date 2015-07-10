# rm -rf output/*

source ../../venv/bin/activate

# Oxford data r1s, thorough 1/2
# ../sparna.py --fwd-reads-sig _F --rev-reads-sig _R --norm-k-list 21,31 --norm-cov-list 1,2,5,15 --asm-k-list 21,33,55,77 --reads-dir data_oxford_r1 --threads 12 --reference /Users/Bede/Research/Datasets/ref_sequences/hcv/h77.fasta

# IVA HIV data
../sparna.py --fwd-reads-sig _1 --rev-reads-sig _2 --norm-k-list 21,31 --norm-cov-list 1,2,5 --asm-k-list 21,33,55,77 --reads-dir data_iva_hiv --threads 6 --reference /Users/bede/Research/Datasets/ref_sequences/hiv/hxb2.fasta

# Fast agnostic HCV multisample
# ../sparna.py --fwd-reads-sig _F --rev-reads-sig _R --norm-k-list 31 --norm-cov-list 1,2,5,10 --asm-k-list 21,33,55,77 --reads-dir data_small --threads 12 --est-ref-len 9678

# Fast agnostic HCV multisample with mapping
# ../sparna.py --fwd-reads-sig _F --rev-reads-sig _R --norm-k-list 31 --norm-cov-list 1,2,5,10 --asm-k-list 21,33,55,77 --reads-dir data_small --threads 12 --reference /Users/Bede/Research/Datasets/ref_sequences/hcv/h77.fasta

# Thorough agnostic HCV multisample
# ../sparna.py --fwd-reads-sig _F --rev-reads-sig _R --norm-k-list 31 --norm-cov-list 1,2,5,10,20 --asm-k-list 21,33,55,77 --reads-dir data --threads 12 --reference /Users/Bede/Research/Datasets/ref_sequences/hcv/h77.fasta

# Not thorough agnostic HCV multisample
# ../sparna.py --fwd-reads-sig _F --rev-reads-sig _R --norm-k-list 31 --norm-cov-list 1,5,20 --asm-k-list 21,33,55,77 --reads-dir data --threads 12 --reference /Users/Bede/Research/Datasets/ref_sequences/hcv/h77.fasta

# Gold standard agnostic HCV multisample
# ../sparna.py --fwd-reads-sig _F --rev-reads-sig _R --norm-k-list 19,25,31 --norm-cov-list 1,2,5,10,20 --asm-k-list 21,33,55,77 --reads-dir data --threads 12 --reference /Users/Bede/Research/Datasets/ref_sequences/hcv/h77.fasta

# Fast agnostic HCV single sample
# ../sparna.py --fwd-reads /Users/Bede/Research/Analyses/phe_asm/phe_hcv_pipeline/tests/data/hcv_n1000_F.fastq --rev-reads /Users/Bede/Research/Analyses/phe_asm/phe_hcv_pipeline/tests/data/hcv_n1000_R.fastq --fwd-reads-sig _F --rev-reads-sig _R --norm-k-list 25,31 --norm-c-list 5,10,15,20 --asm-k-list 33,43  --threads 12

# --reference --asm-using-ref

curl bede.im/notify.php?assembly_finished

tput bel