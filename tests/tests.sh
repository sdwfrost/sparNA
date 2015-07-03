# rm -rf output

source ../../venv/bin/activate

# Agnostic HCV multisample
../pipeline.py --fwd-reads-sig _F --rev-reads-sig _R --norm-k-list 25,31 --norm-c-list 5,10,15,20 --asm-k-list 33,43 --multiple-samples --in-dir data --out-dir output --threads 12

tput bel