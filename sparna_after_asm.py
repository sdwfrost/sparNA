#!/usr/bin/env python

# Author: Bede Constantinides
# Python (2.7+) pipeline for paired-end HepC assembly developed during a placement at PHE
# Please seek permission prior to use or distribution

# TODO
# | Fix read remapping
# | Interleaved reads (ONE TRUE FORMAT)
# | Python3
# | stop using envoy
# | add minimum similarity threshold for reference selection
# | report on trimming, %remapped
# | increase khmer table size
# | TESTS
# | Which reference to use in QUAST... ref_found?

# DEPENDENCIES
# | python packages:
# |    argh, biopython, envoy, khmer, matplotlib
# | others, expected inside $PATH:
# |    bwa, blast, samtools, vcftools, bcftools, bedtools, seqtk, spades, quast, parallel (GNU)
# | others, bundled inside res/ directory:
# |    trimmomatic
# | others, bundled and requiring compilation: segemehl

# USAGE: ./pipeline.py --threads 12 --fwd-reads-sig _F --rev-reads-sig _R --norm-k-list 31 --norm-c-list 5 --asm-k-list 33 --multiple-samples --in-dir /path/to/fastqs --out-dir /path/to/output
# Input fastq filenames should have an extension and a signature to allow identification of forward and reverse reads

# time ./pipeline.py --hcv --fwd-reads-sig _F --rev-reads-sig _R --norm-k-list 25,31 --norm-c-list 5,10,15,20 --asm-k-list 33,43 --multiple-samples --in-dir /Users/Bede/Research/Analyses/phe_asm/phe_hcv_pipeline/input --out-dir /Users/Bede/Research/Analyses/phe_asm/phe_hcv_pipeline/tmp --threads 12

# min_cov
# segemehl -A -D -E values
# min_depth

from __future__ import division, print_function
import os
import sys
import time
import argh
import envoy
import subprocess
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline

def list_fastqs(fwd_reads_sig, rev_reads_sig, paths):
    print('Identifying input...')
    fastqs = {'f':[],'r':[],'fr':[]}
    if paths['in_dir']:
        for fastq in os.listdir(paths['in_dir']):
            if fastq.endswith('.fastq') or fastq.endswith('.fastq'):
                if fwd_reads_sig in fastq:
                    fastqs['f'].append(fastq)
                elif rev_reads_sig in fastq:
                    fastqs['r'].append(fastq)
        fastqs = zip(fastqs['f'], fastqs['r'])
        fastqs = {os.path.splitext(f[0].replace(fwd_reads_sig,''))[0]: f for f in fastqs}
    elif paths['in_fwd'] and paths['in_rev']:
        fastqs['f'].append(os.path.basename(paths['in_fwd']))
        fastqs['r'].append(os.path.basename(paths['in_rev']))
        fastqs = zip(fastqs['f'], fastqs['r'])
        fastqs = {os.path.splitext(f[0].replace(fwd_reads_sig,''))[0]: f for f in fastqs}
    print('\tDone') if fastqs else sys.exit('ERR_READS')
    print('-' * 40)
    return fastqs

def choose_assembly(est_ref_len, sample_name, paths, threads, i=1):
    print('Choosing best assembly...')
    longest_contigs = {}
    contigs_paths = ([
    paths['o'] + '/asm/' +  dir + '/contigs.fasta' for dir in
    filter(lambda d: d.startswith(str(i)), os.listdir(paths['o'] + '/asm'))])
    longest_contigs = {}
    for contigs_path in contigs_paths:
        asm_name = os.path.split(contigs_path)[0].split('/')[-1]
        with open(contigs_path, 'r') as contigs_file:
            longest_contig_name = None
            longest_contig_len = None
            for record in SeqIO.parse(contigs_file, 'fasta'):
                if len(record.seq) > longest_contig_len:
                    longest_contig_len = len(record.seq) 
                    longest_contig_name = record.id
        longest_contigs[asm_name] = (longest_contig_name, longest_contig_len)
    contig_differences = {s: abs(int(est_ref_len)-int(c[1])) for s, c in longest_contigs.items()}
    best_asm = min(contig_differences, key=lambda k: contig_differences[k])
    best_asm_path = paths['o'] + '/asm/' + best_asm + '/contigs.fasta'
    best_contig = longest_contigs[best_asm]
    
    with open(best_asm_path, 'r') as best_asm_file:
        for record in SeqIO.parse(best_asm_file, 'fasta'):
            if record.id == best_contig[0]:
                with open(paths['o'] + '/remap/' + str(i) + '.contig.fasta', 'w') as asm_ref_file:
                    SeqIO.write(record, asm_ref_file, 'fasta')

    print('\tPutative best assembly: ' + best_asm)
    print('\tPutative best contig name: ' + str(best_contig[0]))
    print('\tPutative best contig length: ' + str(best_contig[1]))
    return best_asm, best_contig[0], best_contig[1]

def map_reads_to_assembly(sample_name, paths, threads, i=1):
    print('Aligning to best assembled contig... (Bowtie2)')
    cmd_vars = {
     'i':str(i),
     'sample_name':sample_name,
     'path_pipe':paths['pipe'],
     'path_o':paths['o'],
     'threads':threads}
    cmds = [
     'bowtie2-build -q {path_o}/remap/{i}.contig.fasta {i}.contig &> /dev/null',
     'bowtie2 -x {i} -S {path_o}/remap/{i}.sam --no-unal --threads 12 --local '
     '-1 {path_o}/merge/{i}.{sample_name}.raw.r1.fastq '
     '-2 {path_o}/merge/{i}.{sample_name}.raw.r2.fastq &> /dev/null',
     'grep -v XS:i: {path_o}/remap/{i}.sam > {path_o}/remap/{i}.uniq.sam',
     'samtools view -bS {path_o}/remap/{i}.uniq.sam | samtools sort - {path_o}/remap/{i}.uniq',
     'samtools index {path_o}/remap/{i}.uniq.bam',
     'samtools mpileup -d 1000 -f {path_o}/remap/{i}.contig.fasta {path_o}/remap/{i}.uniq.bam '
     '2> /dev/null > {path_o}/remap/{i}.uniq.pile',
     'samtools mpileup -ud 1000 -f {path_o}/remap/{i}.contig.fasta {path_o}/remap/{i}.uniq.bam '
     '2> /dev/null | bcftools call -c | vcfutils.pl vcf2fq '
     '| seqtk seq -a - > {path_o}/remap/{i}.denovo.consensus.fasta']
    for j, cmd in enumerate(cmds, start=1):
        cmd_map = os.system(cmd.format(**cmd_vars))
        print('\tDone (' + cmd.split(' ')[0] + ')') if cmd_map == 0 else sys.exit('ERR_REMAP')

def assess_remap_coverage(best_contig_len, paths, i=1):
    print('Identifying low coverage regions...')
    depths = {}
    max_coverage = None
    uncovered_sites = []
    with open(paths['o'] + '/remap/' + str(i) + '.uniq.pile', 'r') as pileup:
        bases_covered = 0
        for line in pileup:
            site = int(line.split('\t')[1])
            depths[site] = int(line.split('\t')[3])
            if depths[site] < 1:
                uncovered_sites.append(site)
            else:
                bases_covered += 1
            if depths[site] > max_coverage:
                max_coverage = depths[site]
    prop_coverage = round(bases_covered/best_contig_len, 3)
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
    
    print('\tBases covered: {}/{} ({})'.format(bases_covered, best_contig_len, prop_coverage))
    print('\tMaximum depth of coverage: ' + str(max_coverage))
    print('\tUncovered regions: ' + str(len(uncovered_regions)))
    print('\tLargest uncovered region: ' + str(largest_uncovered_region) + 'bp')

# First arg previosly ref_found
def evaluate_assemblies(reference, est_ref_len, sample_name, paths, threads, i=1):
    print('Comparing assemblies...')
    asm_dirs = (
     [paths['o'] + '/asm/' + dir + '/contigs.fasta' for dir in
     filter(lambda d: d.startswith(str(i)), os.listdir(paths['o'] + '/asm'))])
    print(asm_dirs)
    cmd_vars = {
     'i':str(i),
     'asm_dirs':' '.join(asm_dirs),
     'ref_len':est_ref_len,
     'path_ref':paths['ref'],
     'sample_name':sample_name,
     'path_o':paths['o'],
     'threads':threads}
    cmd_eval = (
     'quast.py {asm_dirs} -o {path_o}/eval/{i}.{sample_name} '
     '--threads {threads} '
     # '--gene-finding'
     .format(**cmd_vars))
    # if ref_found:
    #     cmd_eval += ' -R {path_o}/ref/{i}.{sample_name}.ref.fasta'.format(**cmd_vars)
    if reference:
        cmd_eval += ' -R {path_ref}'.format(**cmd_vars)
    if est_ref_len:
        cmd_eval += ' --est-ref-size {ref_len}'.format(**cmd_vars)
    cmd_eval += ' --min-contig 50' # Just for testing
    cmd_eval += ' &> /dev/null'
    # print(cmd_eval)
    cmd_eval = os.system(cmd_eval)
    print('\tDone') if cmd_eval == 0 else sys.exit('ERR_EVAL')

def evaluate_all_assemblies(reference, est_ref_len, sample_name, paths, threads, i=1):
    print('Comparing all assemblies...')
    asm_dirs = (
    [paths['o'] + '/asm/' + dir + '/contigs.fasta' for dir in
    filter(lambda d: d[0].isdigit(), os.listdir(paths['o'] + '/asm'))])
    print(asm_dirs)
    cmd_vars = {
     'i':str(i),
     'asm_dirs':' '.join(asm_dirs),
     'ref_len':est_ref_len,
     'path_ref':paths['ref'],
     'sample_name':sample_name,
     'path_o':paths['o'],
     'threads':threads}
    cmd_eval = (
     'quast.py {asm_dirs} -o {path_o}/eval/ '
     '--threads {threads} '
     '--gene-finding'.format(**cmd_vars))
    if reference:
        cmd_eval += ' -R {path_ref}'.format(**cmd_vars)
    if est_ref_len:
        cmd_eval += ' --est-ref-size {ref_len}'.format(**cmd_vars)
    cmd_eval += ' --min-contig 50' # Just for testing
    cmd_eval += ' &> /dev/null'
    cmd_eval = os.system(cmd_eval)
    print('\tDone') if cmd_eval == 0 else sys.exit('ERR_EVAL_SUMMARY')

def report(start_time, end_time, paths):
    elapsed_time = end_time - start_time
    with open(paths['o'] + '/summary.txt', 'w') as report:
        report.write('wall_time\t{}'.format(elapsed_time))
        print('\tWall time: {0:.{1}f}s'.format(elapsed_time, 1))

def main(fwd_reads=None, rev_reads=None, reads_dir=None, out_dir='output', fwd_reads_sig='_F',
    rev_reads_sig='_R', norm_k_list=None, norm_cov_list=None, asm_k_list=None,
    use_segemehl=False, reference=None, reference_guided_asm=False, est_ref_len=None,
    map_before_asm=False, map_after_asm=False,
    hcv=False, threads=1):
   
    multiple_samples = True if reads_dir else False
    if not est_ref_len and reference:
        with open(reference, 'r') as reference_file:
            record = SeqIO.parse(reference_file, 'fasta')
            est_ref_len = len(record[0].seq)

    print('-' * 40)
    print('Run options...')
    print('\tMultiple input samples') if multiple_samples else print('\tSingle sample')
    print('\tPaired read signatures: \'' + fwd_reads_sig + '\', \'' + rev_reads_sig + '\'')
    print('\tVirus agnostic') if not hcv else print('\tUsing HCV-specific features')
    print('\tReference length: ' + str(est_ref_len))
    print('\tUsing Segemehl') if use_segemehl else print('\tUsing BWA')
    print('\t' + str(threads) + ' threads available')

    paths = {
    'in_dir':reads_dir,
    'in_fwd':fwd_reads,
    'in_rev':rev_reads,
    'ref':reference,
    'pipe':os.path.dirname(os.path.realpath(__file__)),
    'o':out_dir + '/run_' + str(1436296420) + '_gold_no_remap',
    'segemehl':os.path.dirname(os.path.realpath(__file__)) + '/res/segemehl/segemehl.x'}

    state = {
    'fastqs':None,
    'n_reads':None,
    'n_reads_sample':None,
    'ref_found':None,
    'top_ref_accession':None,
    'top_hcv_genotype':None,
    'ref_path':None,
    'ref_len':None}

    state['fastqs'] = list_fastqs(fwd_reads_sig, rev_reads_sig, paths)
    i = 0 # counter for file naming
    for sample_name, fastq_names in state['fastqs'].items():
        i += 1
        best_asm, best_contig_id, best_contig_len = choose_assembly(est_ref_len, sample_name, paths, threads, i)
        map_reads_to_assembly(sample_name, paths, threads, i)
        best_contig_len = 9471
        assess_remap_coverage(best_contig_len, paths, i)
        break
        # evaluate_assemblies(reference, est_ref_len, sample_name, paths, threads, i)
    # evaluate_all_assemblies(reference, est_ref_len, sample_name, paths, threads, i)

argh.dispatch_command(main)