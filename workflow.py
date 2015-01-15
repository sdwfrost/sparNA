#!/usr/bin/env python

# AUTHOR: Bede Constantinides
# SPAdes pipeline for paired-end HepC assembly developed on placement at PHE

# TODO

# | use khmer for deinterleaving (split-paired-reads.py)
# | replace envoy calls with a less broken subprocess wrapper
# | add minimum similarity threshold for reference selection
# | Build command strings using  cleaner printf() style variable susbtituion
# | decent error handling at each stage
# | tests for dependencies
# | command line params
# | mauve/nucmer integration for when a contiguous assembly is unavailable

# DEPENDENCIES
# | python packages:
# |   biopython, envoy, khmer
# | others expected inside $PATH:
# |   blast, samtools, seqtk, spades, quast
# | bundled inside res/ directory:
# |   trimmomatic, fastq_deinterleave
# | bundled but requiring compilation: segemehl

from __future__ import division, print_function
import os
import sys
import time
import envoy
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline

n_threads = 8
pipe_path = os.path.dirname(os.path.realpath(__file__))
job_id = str(int(time.time()))
job_path = pipe_path + '/tmp/' + job_id
os.mkdir(job_path)
os.mkdir(job_path + '/merge/')
os.mkdir(job_path + '/sample/')
os.mkdir(job_path + '/blast/')
os.mkdir(job_path + '/ref/')
os.mkdir(job_path + '/map/')
os.mkdir(job_path + '/trim/')
os.mkdir(job_path + '/norm/')
os.mkdir(job_path + '/asm/')
os.mkdir(job_path + '/remap/')

def merge_reads(fwd_reads_sig, rev_reads_sig):
	print('Merging reads... ')
	fastqs = {'f':[], 'r':[]}
	for fastq in os.listdir(pipe_path + '/input/'):
		if fastq.endswith('.fastq') or fastq.endswith('.fastq'):
			if fwd_reads_sig in fastq:
				fastqs['f'].append(pipe_path + '/input/' + fastq)
			elif rev_reads_sig in fastq:
				fastqs['r'].append(pipe_path + '/input/' + fastq)
	cmd_merge = os.system(''
	+ 'cat ' + ' '.join(fastqs['f']) + ' > ' + job_path + '/merge/raw.r1.fastq && '
	+ 'cat ' + ' '.join(fastqs['r']) + ' > ' + job_path + '/merge/raw.r2.fastq && '
	+ 'interleave-reads.py ' + job_path + '/merge/raw.r1.fastq ' + job_path
	+ '/merge/raw.r2.fastq > ' + job_path + '/merge/raw.r12.fastq')

def count_reads():
	print('Counting reads...')
	cmd_count = envoy.run('wc -l ' + job_path + '/merge/raw.r12.fastq')
	print(cmd_count.std_out)
	n_reads = int(cmd_count.std_out.replace(' ','').split('/')[0])/4
	return n_reads

def sample_reads(n_reads):
	print('Sampling reads...')
	n_reads_sample = n_reads/100
	cmd_sample = os.system(''
	+ 'cat ' + job_path + '/merge/raw.r12.fastq | seqtk sample - ' + str(n_reads_sample)
	+ ' | seqtk seq -a - > ' + job_path + '/sample/sample.fasta')
	return n_reads_sample

def choose_reference(n_reads):
	print('Choosing optimal reference sequence...')
	if not os.path.exists(pipe_path + '/res/hcv_db/db.nhr'):
		cmd_blastn_index = envoy.run(''
		+ 'makeblastdb -dbtype nucl -input_type fasta -in '
		+ pipe_path + '/res/hcv_db/db.fasta -title db')
	cmd_blastn = NcbiblastnCommandline(
	query = job_path + '/sample/sample.fasta',
	db = pipe_path + '/res/hcv_db/db.fasta',
	out = job_path + '/blast/blast.tsv',
	evalue = 1e-10,
	outfmt = 7,
	num_alignments = 1,
	num_threads = n_threads)
	cmd_blastn()

	accession_freqs = {}
	genotype_freqs = {}
	with open(job_path + '/blast/blast.tsv', 'r') as blast_out:
		for line in blast_out:
			if not line.startswith('#'):
			 	accession = line.split('\t')[1]
				genotype = accession.split('_')[1].split('.')[0]
				if accession in accession_freqs.keys():
					accession_freqs[accession] += 1
				else: accession_freqs[accession] = 1
				if genotype in genotype_freqs.keys():
					genotype_freqs[genotype] += 1
				else: genotype_freqs[genotype] = 1
	genotype_composition = {key: genotype_freqs[key]/100 for key in genotype_freqs.keys()}
	best_accession = max(accession_freqs, key=accession_freqs.get)
	genotype_freqs_sorted = reversed(sorted(genotype_freqs.items(),key=lambda x:x[1]))
	print('Strain represention:')
	for items in genotype_freqs_sorted:
		print(items[0] + ': ' + str(round(items[1]/n_reads*100,2)) + '%')

	reference = ''
	print('Fetching chosen reference sequence...')
	with open(pipe_path + '/res/hcv_db/db.fasta', 'r') as references_fa:
		is_inside_best_reference = False
		for line in references_fa:
			if line.startswith('>'):
				if best_accession in line:
					is_inside_best_reference = True
				else: is_inside_best_reference = False
			elif is_inside_best_reference:
				reference += line.strip()
	reference_len = len(reference)
	reference_fa_path = job_path + '/ref/ref.fasta'
	with open(reference_fa_path, 'w') as reference_fa:
		reference_fa.write('>' + best_accession + '\n' + reference)
	return best_accession, reference_fa_path,  reference_len

def map_reads(reference_fa_path):
	print('Aligning with Segemehl... ')
	cmds_align = [
		pipe_path + '/res/segemehl/segemehl.x -d ' + reference_fa_path + ' -x '
		+ reference_fa_path + '.idx',
		pipe_path + '/res/segemehl/segemehl.x -d ' + reference_fa_path + ' -x '
		+ reference_fa_path + '.idx' + ' -q ' + job_path + '/merge/raw.r12.fastq'
		+ ' --threads ' + str(n_threads) + ' -A 80 > ' + job_path + '/map/segemehl_mapped.sam',
		'samtools view -bS ' + job_path + '/map/segemehl_mapped.sam | samtools sort - '
		+ job_path + '/map/segemehl_mapped; ' + 'samtools index ' + job_path
		+ '/map/segemehl_mapped.bam',
		'samtools mpileup -d 1000 -f ' + reference_fa_path + ' ' + job_path
		+ '/map/segemehl_mapped.bam > ' + job_path + '/map/segemehl_mapped.pile',
		'samtools mpileup -ud 1000 -f ' + reference_fa_path + ' ' + job_path + '/map/segemehl_mapped.bam | '
		+ 'bcftools call -c | vcfutils.pl vcf2fq | seqtk seq -a - | fasta_formatter '
		+ '-o ' + job_path + '/map/consensus.fasta']
	for cmd_align in cmds_align:
		os.system(cmd_align)
		# x = envoy.run(cmd_align)
		# print(x.std_out) if x.status_code == 0 else sys.exit(x.std_err)

def assess_coverage(reference_len):
	print('Identifying low coverage regions... ')
	min_depth = 1
	min_coverage = 0.9
	depths = {}
	uncovered_sites = []
	with open(job_path + '/map/segemehl_mapped.pile', 'r') as pileup:
		bases_covered = 0
		for line in pileup:
			site = int(line.split('\t')[1])
			depths[site] = int(line.split('\t')[3])
			if depths[site] < min_depth:
				uncovered_sites.append(site)

	uncovered_region = 0
	uncovered_regions = []
	last_uncovered_site = 0
	largest_uncovered_region = 0
	for uncovered_site in uncovered_sites:
		if uncovered_site == last_uncovered_site + 1:
			uncovered_region += 1
			if uncovered_region > largest_uncovered_region:
				largest_uncovered_region = uncovered_region
		else:
			if uncovered_region > 0:
				uncovered_regions.append(uncovered_region)
			uncovered_region = 1
		last_uncovered_site = uncovered_site

	print('Uncovered sites: ' + str(len(uncovered_sites)))
	print('Uncovered regions: ' + str(len(uncovered_regions)))

	if not uncovered_sites:
		print('All [reference] bases covered!')
	elif len(uncovered_sites) < (1-min_coverage)*reference_len:
		print('Reference coverage safely above threshold')
	else:
		print('Reference coverage below threshold')


def trim():
	print('Trimming with Trimmomatic... ')
	cmd_trim = envoy.run(''
	+ 'java -jar ' + pipe_path + '/res/trimmomatic-0.32.jar PE ' + job_path
	+ '/merge/raw.r1.fastq ' + job_path + '/merge/raw.r2.fastq ' + job_path
	+ '/trim/trim.r1_pe.fastq ' + job_path + '/trim/trim.r1_se.fastq ' + job_path
	+ '/trim/trim.r2_pe.fastq ' + job_path + '/trim/trim.r2_se.fastq ILLUMINACLIP:'
	+ pipe_path + '/res/illumina_adapters.fa:2:30:10 MINLEN:25')
	print(cmd_trim.std_out, cmd_trim.std_err)
	os.system(''
	+ 'cat ' + job_path + '/trim/trim.r1_se.fastq ' + job_path + '/trim/trim.r2_se.fastq > '
	+ job_path + '/trim/trim.se.fastq')
	os.system(''
	+ 'interleave-reads.py ' + job_path + '/trim/trim.r1_pe.fastq ' + job_path + '/trim/trim.r2_pe.fastq > '
	+ job_path + '/trim/trim.r12_pe.fastq')
	cmd_trim_stats = ''.join(cmd_trim.std_err).split('\n')[25]

def normalise():
	print('Normalising... ')
	k_params = [31]
	c_params = [5]
	norm_perms = [{'k':k, 'c':c} for k in k_params for c in c_params]
	for norm_perm in norm_perms:
		k, c = str(norm_perm['k']), str(norm_perm['c'])
		cmd_norm = os.system(''
		+ 'normalize-by-median.py -C ' + c + ' -k ' + k + ' -N 1 -x 1e9 -p ' + job_path + '/trim/trim.r12_pe.fastq '
		+ '-o ' + job_path + '/norm/norm_k' + k + 'c' + c + '.r12_pe.fastq && > /dev/null 2>&1 && '
		+ 'normalize-by-median.py -C ' + c + ' -k ' + k + ' -N 1 -x 1e9 ' + job_path + '/trim/trim.se.fastq '
		+ '-o ' + job_path + '/norm/norm_k' + k + 'c' + c + '.se.fastq && > /dev/null 2>&1 && '
		+ pipe_path + '/res/fastq_deinterleave ' + job_path + '/norm/norm_k' + k + 'c' + c + '.r12_pe.fastq '
		+ job_path + '/norm/norm_k' + k + 'c' + c + '.r1_pe.fastq '
		+ job_path + '/norm/norm_k' + k + 'c' + c + '.r2_pe.fastq && '
		+ 'cat ' + job_path + '/norm/norm_k' + k + 'c' + c + '.r12_pe.fastq '
		+ job_path + '/norm/norm_k' + k + 'c' + c + '.se.fastq > '
		+ job_path + '/norm/norm_k' + k + 'c' + c + '.pe_and_se.fastq')
		print('...Khmer completed (k=' + k + ', c=' + c + ')') if cmd_norm == 0 else sys.exit('FAIL')
	print('DONE')
	return norm_perms

def assemble(norm_perms):
	k_params_str = '21,33,43,55,77,99'
	untrusted_contigs_params = [False]
	asm_perms = [{'k':n['k'], 'c':n['c'], 'uc':uc} for n in norm_perms for uc in untrusted_contigs_params]
	for asm_perm in asm_perms:
		print('Assembling with SPAdes... ')
		k, c, uc = str(asm_perm['k']), str(asm_perm['c']), str(asm_perm['uc'])
		norm_path_prefix = job_path + '/norm/norm_k' + k + 'c' + c		
		asm_path = job_path + '/asm/asm_uc' + uc + '.norm_k' + k + 'c' + c
		cmd_asm = (''
		+ 'spades.py -m 8 -t ' + str(n_threads) + ' -k ' + k_params_str + ' '
		+ '--pe1-1 ' + norm_path_prefix + '.r1_pe.fastq '
		+ '--pe1-2 ' + norm_path_prefix + '.r2_pe.fastq '
		+ '--s1 ' + norm_path_prefix + '.se.fastq '
		+ '-o ' + asm_path + ' --careful')
		if asm_perm['uc']:
			cmd_asm += ' --untrusted-contigs ' + job_path + '/ref/ref.fasta'
		cmd_asm = envoy.run(cmd_asm)

		print('...SPAdes completed ') if cmd_asm.status_code == 0 else sys.exit('...Error (SPAdes) \n' + cmd_asm.std_out +  '\n' + cmd_asm.std_err)

def choose_assembly():
	asm_paths = [job_path + '/asm/' + dir + '/contigs.fasta' for dir in os.listdir(job_path + '/asm') if dir.startswith('asm')]
	print(asm_paths)
	os.system('quast.py ' + ' '.join(asm_paths) + ' -R ' + job_path + '/ref/ref.fasta -o ' + job_path
	+ '/eval' + ' --threads ' + str(n_threads))

def remap_reads():
	pass

def report():
	os.system('open ' + job_path + '/eval/report.html')

def run_pipeline():
	# fwd_reads_sig = '.R1.'
	# rev_reads_sig = '.R2.'
	# fwd_reads_sig = '_F.'
	# rev_reads_sig = '_R.'
	fwd_reads_sig = '_R1.'
	rev_reads_sig = '_R2.'
	
	merge_reads(fwd_reads_sig, rev_reads_sig)
	n_reads = count_reads()
	sample_reads(n_reads)
	reference_accession, reference_path, reference_len = choose_reference(n_reads)
	trim()
	assemble(normalise())
	choose_assembly()
	remap_reads()
	report()

run_pipeline()
