# rm -rf output/*

source ../../venv/bin/activate

../sparna_after_asm.py --fwd-reads-sig _F --rev-reads-sig _R --norm-k-list 31 --norm-cov-list 1,2,5,10 --asm-k-list 21,33,55,77 --reads-dir data --threads 12 --est-ref-len 9678

tput bel