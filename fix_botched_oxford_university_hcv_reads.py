#!/usr/bin/env python3
# Fix pairs where /1s and /2s are correctly ordered but mixed between both files
# Potentially destructive but should cause not loss of information in non-SS seq libraries

import os
import sys
from Bio import SeqIO

path_to_reads = sys.argv[1]
if not os.path.exists(path_to_reads + '/fixed'):
    os.mkdir(path_to_reads + '/fixed')

fwd_fastqs = [fn for fn in os.listdir(path_to_reads) if fn.endswith('_F.fastq')]
rev_fastqs = [fn for fn in os.listdir(path_to_reads) if fn.endswith('_R.fastq')]
fastq_pairs = zip(fwd_fastqs, rev_fastqs)

for fastq_pair in fastq_pairs:
	with open(path_to_reads + '/' +  fastq_pair[0], 'rU') as fwd_fastq:
		with open(path_to_reads + '/fixed/' +  fastq_pair[0], 'w') as fixed_fwd_fastq:
			fixed_fwd_records = []
			for fwd_record in SeqIO.parse(fwd_fastq, 'fastq'):
				fwd_record.description = fwd_record.description.replace('/2','/1')
				fixed_fwd_records.append(fwd_record)
			SeqIO.write(fixed_fwd_records, fixed_fwd_fastq, 'fastq')
	# ...