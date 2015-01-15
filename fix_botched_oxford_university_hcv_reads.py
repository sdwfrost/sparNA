#!/usr/bin/env python3
# Fix pairs where /1s and /2s are correctly ordered but mixed between both files
# Potentially destructive but should cause not loss of information in non-SS seq libraries

import os
import sys
from Bio import SeqIO

path_to_reads = sys.argv[1]
if not os.path.exists(path_to_reads + '/fixed'):
    os.mkdir(path_to_reads + '/fixed')

fwd_fastqs = [fn for fn in os.listdir(path_to_reads) if fn.endswith('_F.fastq].fastq')]
rev_fastqs = [fn for fn in os.listdir(path_to_reads) if fn.endswith('_R.fastq].fastq')]
fastq_pairs = zip(fwd_fastqs, rev_fastqs)

for fastq_pair in fastq_pairs:
	with open(path_to_reads + '/' +  fastq_pair[0], 'rU') as fwd_fastq:
		with open(path_to_reads + '/fixed/' +  fastq_pair[0], 'w') as fixed_fwd_fastq:
			fixed_fwd_records = []
			for fwd_record in SeqIO.parse(fwd_fastq, 'fastq'):
				fwd_record.id = fwd_record.id.replace('/2','/1')
				fwd_record.description = ''
				fixed_fwd_records.append(fwd_record)
			SeqIO.write(fixed_fwd_records, fixed_fwd_fastq, 'fastq')
			print('Processed ' + fastq_pair[0])

	with open(path_to_reads + '/' +  fastq_pair[1], 'rU') as rev_fastq:
		with open(path_to_reads + '/fixed/' +  fastq_pair[1], 'w') as fixed_rev_fastq:
			fixed_rev_records = []
			for rev_record in SeqIO.parse(rev_fastq, 'fastq'):
				rev_record.id = rev_record.id.replace('/1','/2')
				rev_record.description = ''
				fixed_rev_records.append(rev_record)
			SeqIO.write(fixed_rev_records, fixed_rev_fastq, 'fastq')
			print('Processed ' + fastq_pair[1])

# # Tests
# with open(path_to_reads + '/fixed/' +  fastq_pair[0], 'rU') as fixed_rev_fastq:
# 	for line in fixed_rev_fastq:
# 		header = line.strip()
# 		if header.startswith('@MIS') and ('/2' in header or len(line) > 60):
# 			print('1> ' + header)

# with open(path_to_reads + '/fixed/' +  fastq_pair[1], 'rU') as fixed_rev_fastq:
# 	for line in fixed_rev_fastq:
# 		header = line.strip()
# 		if header.startswith('@MIS') and ('/1' in header or len(line) > 60):
# 			print('2> ' + header)


