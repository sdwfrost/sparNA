#!/usr/bin/env python

# Copyright 2015 Bede Constantinides, University of Manchester (b at bede dot im)
# Python 3.5+ pipeline for assembly of Illumina PE RNA virus amplicons
# Developed in collaboration with Public Health England, Colindale

# TODO
# | Handle no reads mapped by BWA
# | Decide on BWA vs Bowtie2
# | remove Python2 prefix from spades.py call
# | Top100 contigs only? --eval-n-longest-contigs 100
# | GZIP support
# | Fully migrate to bwa for mapping... Base reports on samtools flagsts?
# | best_contig_len undefined
# | Send only best contigs per sample to QUAST for final eval step
# | Write and parse and report Bowtie2 output map and remap
# | Send coverage stats to file
# | Report % read alignment in mapping to ref and contig
# | More pipelining to reduce disk I/O (mainly SAMtools)
# | Consistent use of r12 / fr
# | Fix read remapping
# | Interleaved reads (ONE TRUE FORMAT)
# | add minimum similarity threshold for reference selection
# | report on trimming, %remapped
# | increase khmer table size
# | TESTS
# | Which reference to use in QUAST... ref_found?
# | Usage
# | Determine best asm by mapping reads to all assemblies?
# | Bootstrap/dogfoood assemblies with --trusted-contigs etc

# DEPENDENCIES
# | python packages:
# |    argh, biopython, envoy, khmer, matplotlib
# | others, expected inside $PATH:
# |    bwa, bowtie2, blast, samtools, vcftools, bcftools, bedtools, seqtk, spades, quast, parallel (GNU)
# | others, bundled inside res/ directory:
# |    trimmomatic
# | others, bundled and requiring compilation:

# Input fastq filenames should have an extension and a signature to allow identification of forward and reverse reads
# min_cov
# min_depth

import os
import sys
import time
import argh
import logging
import subprocess
import multiprocessing

from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def run(cmd):
    return subprocess.run(cmd,
                          shell=True,
                          universal_newlines=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT)


def list_fastqs(fwd_reads_sig, rev_reads_sig, paths):
    print('-' * 40)
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
    return fastqs


def import_reads(multiple_samples, sample_name, fastq_names, paths, i=1):
    print('-' * 40)
    print('Importing reads...')
    print('\tSample ID: ' + str(i))
    print('\tSample name: ' + sample_name)
    print('\tFiles:', fastq_names)
    if multiple_samples:
        fastq_path = paths['in_dir']
    else:
        fastq_path = os.path.split(paths['in_fwd'])[0]
    cmd_vars = {
    'i':str(i),
    'sample_name':sample_name,
    'fastq_path_f':fastq_path + '/' + fastq_names[0],
    'fastq_path_r':fastq_path + '/' + fastq_names[1],
    'path_o':paths['o'],
    'path_pipe':paths['pipe']}
    cmd = (
    'cp {fastq_path_f} {path_o}/merge/{i}.{sample_name}.raw.r1.fastq && '
    'cp {fastq_path_r} {path_o}/merge/{i}.{sample_name}.raw.r2.fastq && '
    # 'interleave-reads.py {path_o}/merge/{i}.{sample_name}.raw.r1.fastq '
    '{path_pipe}/res/interleave.py {path_o}/merge/{i}.{sample_name}.raw.r1.fastq '
    '{path_o}/merge/{i}.{sample_name}.raw.r2.fastq > '
    '{path_o}/merge/{i}.{sample_name}.raw.r12.fastq'
    .format(**cmd_vars))
    logger.info(cmd)
    cmd_run = run(cmd)
    logger.info(cmd_run.stdout)
    print('\tDone') if cmd_run.returncode == 0 else sys.exit('ERR_IMPORT')


def count_reads(sample_name, paths, i=1):
    print('Counting reads...')
    cmd_count = ('wc -l {path_o}/merge/{i}.{sample_name}.raw.r12.fastq'
    .format(i=str(i),
            path_o=paths['o'],
            sample_name=sample_name))
    cmd_run = run(cmd_count)
    n_reads = int(cmd_run.stdout.strip().split(' ')[0])/4
    print('\tDone') if cmd_run.returncode == 0 else sys.exit('ERR_COUNT')
    return n_reads


def map_to_reference(ref, sample_name, paths, threads, i=1):
    print('Aligning... (BWA)')
    cmd_vars = {
     'i':str(i),
     'sample_name':sample_name,
     'path_pipe':paths['pipe'],
     'path_o':paths['o'],
     'ref':ref,
     'threads':threads}
    cmds = [
     'cp {ref} {path_o}/premap/ref.fa',
     'bwa index {path_o}/premap/ref.fa',
     'bwa mem -t {threads} {path_o}/premap/ref.fa {path_o}/merge/{i}.{sample_name}.raw.r12.fastq '
     '> {path_o}/premap/{i}.{sample_name}.mapped.sam ',
     'samtools view -bS {path_o}/premap/{i}.{sample_name}.mapped.sam '
     '| samtools sort - {path_o}/premap/{i}.{sample_name}.mapped',
     'samtools index {path_o}/premap/{i}.{sample_name}.mapped.bam',
     'samtools mpileup -d 1000 '
     '-f {path_o}/premap/ref.fa {path_o}/premap/{i}.{sample_name}.mapped.bam ' 
     '> {path_o}/premap/{i}.{sample_name}.mapped.pile ',
     'samtools mpileup -ud 1000 '
     '-f {path_o}/premap/ref.fa {path_o}/premap/{i}.{sample_name}.mapped.bam ' 
     '| bcftools call -c | vcfutils.pl vcf2fq '
     '| seqtk seq -a - > {path_o}/premap/{i}.{sample_name}.consensus.fasta']
    cmds = [cmd.format(**cmd_vars) for cmd in cmds]
    for cmd in cmds:
        logger.info(cmd)
        cmd_run = run(cmd)
        logger.info(cmd_run.stdout)
        print('\tDone') if cmd_run.returncode == 0 else sys.exit('ERR_PREMAP')
    # Return mapping stats


def trim(sample_name, paths, i=1):
    print('Trimming...')
    cmd_vars = {
     'i': str(i), 
     'path_pipe':paths['pipe'],
     'path_o':paths['o'],
     'sample_name':sample_name}
    cmds = [
     'java -jar {path_pipe}/res/trimmomatic-0.32.jar PE '
     '{path_o}/merge/{i}.{sample_name}.raw.r1.fastq '
     '{path_o}/merge/{i}.{sample_name}.raw.r2.fastq '
     '{path_o}/trim/{i}.{sample_name}.trim.r1_pe.fastq '
     '{path_o}/trim/{i}.{sample_name}.trim.r1_se.fastq '
     '{path_o}/trim/{i}.{sample_name}.trim.r2_pe.fastq '
     '{path_o}/trim/{i}.{sample_name}.trim.r2_se.fastq '
     'ILLUMINACLIP:{path_pipe}/res/illumina_adapters.fa:2:30:10 MINLEN:25',
     'cat {path_o}/trim/{i}.{sample_name}.trim.r1_se.fastq '
     '{path_o}/trim/{i}.{sample_name}.trim.r2_se.fastq '
     '> {path_o}/trim/{i}.{sample_name}.trim.se.fastq '
     '&& {path_pipe}/res/interleave.py {path_o}/trim/{i}.{sample_name}.trim.r1_pe.fastq '
     '{path_o}/trim/{i}.{sample_name}.trim.r2_pe.fastq > '
     '{path_o}/trim/{i}.{sample_name}.trim.r12_pe.fastq']
    cmds = [cmd.format(**cmd_vars) for cmd in cmds]
    for cmd in cmds:
        logger.info(cmd)
        cmd_run = run(cmd)
        logger.info(cmd_run.stdout)
        print('\tDone') if cmd_run.returncode == 0 else sys.exit('ERR_TRIM')
    # Return trim stats?


def normalise(norm_k_list, norm_cov_list, sample_name, paths, threads, i=1):
    print('Normalising...')
    ks = norm_k_list.split(',')
    cs = norm_cov_list.split(',')
    norm_perms = [{'k':k, 'c':c} for k in ks for c in cs]
    cmds = []
    for norm_perm in norm_perms:
        cmd_vars = {
         'i':str(i),
         'k':str(norm_perm['k']),
         'c':str(norm_perm['c']),
         'path_pipe':paths['pipe'],
         'path_o':paths['o'],
         'sample_name':sample_name}
        cmd = (
         'normalize-by-median.py -C {c} -k {k} -N 4 -x 1e9 -p '
         '{path_o}/trim/{i}.{sample_name}.trim.r12_pe.fastq '
         '-o {path_o}/norm/{i}.{sample_name}.norm_k{k}c{c}.r12_pe.fastq '
         '&& normalize-by-median.py -C {c} -k {k} -N 1 -x 1e9 '
         '{path_o}/trim/{i}.{sample_name}.trim.se.fastq '
         '-o {path_o}/norm/{i}.{sample_name}.norm_k{k}c{c}.se.fastq '
         '&& split-paired-reads.py '
         '-1 {path_o}/norm/{i}.{sample_name}.norm_k{k}c{c}.r1_pe.fastq '
         '-2 {path_o}/norm/{i}.{sample_name}.norm_k{k}c{c}.r2_pe.fastq '
         '{path_o}/norm/{i}.{sample_name}.norm_k{k}c{c}.r12_pe.fastq '
         '&& cat {path_o}/norm/{i}.{sample_name}.norm_k{k}c{c}.r12_pe.fastq '
         '{path_o}/norm/{i}.{sample_name}.norm_k{k}c{c}.se.fastq > '
         '{path_o}/norm/{i}.{sample_name}.norm_k{k}c{c}.pe_and_se.fastq'
         .format(**cmd_vars))
        cmds.append(cmd)
        print('\tNormalising norm_k={k},norm_c={c}'.format(**cmd_vars))
        logger.info('Normalising norm_k={k},norm_c={c}'.format(**cmd_vars))
    with multiprocessing.Pool(threads) as pool:
        results = pool.map(run, cmds)
    logger.info([result.stdout for result in results])
    print('\tDone')
    # with open(os.devnull, 'w') as devnull:
    #     processes = [subprocess.Popen(cmd, shell=True, stdout=devnull) for cmd in cmds_norm]
    #     for process in processes:
    #         process.wait()
    #         print('\tDone') if process.returncode == 0 else sys.exit('ERR_NORM')
    return norm_perms


def assemble(norm_perms, asm_k_list, reference_guided_asm, reference, sample_name, paths, threads, i=1):
    print('Assembling...')
    if reference and reference_guided_asm:
        asm_perms = [{'k':p['k'],'c':p['c'],'rg':rg} for p in norm_perms for rg in [1, 0]]
    else:
        asm_perms = [{'k':p['k'],'c':p['c'],'rg':rg} for p in norm_perms for rg in [0]]
    cmds_asm = []
    for asm_perm in asm_perms:
        cmd_vars = {
         'i':str(i),
         'k':str(asm_perm['k']),
         'c':str(asm_perm['c']),
         'rg':str(asm_perm['rg']),
         'asm_k_list':asm_k_list,
         'path_o':paths['o'],
         'path_ref':paths['ref'],
         'sample_name':sample_name,
         'threads':threads}
        cmd_asm = (
         'python2 /usr/local/bin/spades.py -m 8 -t {threads} -k {asm_k_list} '
         '--pe1-1 {path_o}/norm/{i}.{sample_name}.norm_k{k}c{c}.r1_pe.fastq '
         '--pe1-2 {path_o}/norm/{i}.{sample_name}.norm_k{k}c{c}.r2_pe.fastq '
         '--s1 {path_o}/norm/{i}.{sample_name}.norm_k{k}c{c}.se.fastq '
         '-o {path_o}/asm/{i}.{sample_name}.norm_k{k}c{c}.asm_k{asm_k_list}.rg{rg} --careful'
         .format(**cmd_vars))
        if asm_perm['rg']:
            cmd_asm += ' --untrusted-contigs {path_ref}'.format(**cmd_vars)
        cmds_asm.append(cmd_asm)
        print('\tAssembling norm_k={k},norm_c={c},asm_k={asm_k_list},rg={rg}'.format(**cmd_vars))
    with open(os.devnull, 'w') as devnull:
        processes = [subprocess.Popen(cmd, shell=True, stdout=devnull) for cmd in cmds_asm]
        for process in processes:
            process.wait()
            print('\tDone') if process.returncode == 0 else sys.exit('ERR_ASM')


def choose_assembly(target_genome_len, sample_name, paths, threads, i=1):
    # Perhaps should map reads to multiple assemblies?
    print('Choosing best assembly...')
    longest_contigs = {}
    contigs_paths = (
    [paths['o'] + '/asm/' +  dir + '/contigs.fasta' for dir in
    filter(lambda d: d.startswith(str(i)), os.listdir(paths['o'] + '/asm'))])
    longest_contigs = {}
    for contigs_path in contigs_paths:
        asm_name = os.path.split(contigs_path)[0].split('/')[-1]
        with open(contigs_path, 'r') as contigs_file:
            longest_contig_name = None
            longest_contig_len = 0
            for record in SeqIO.parse(contigs_file, 'fasta'):
                if len(record.seq) > longest_contig_len:
                    print('truee')
                    longest_contig_len = len(record.seq) 
                    longest_contig_name = record.id
        longest_contigs[asm_name] = (longest_contig_name, longest_contig_len)
    contig_differences = {s: abs(int(target_genome_len)-int(c[1])) for s, c in longest_contigs.items()}
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


def map_to_assembly(sample_name, paths, threads, i=1):
    print('Aligning to best assembled contig... (Bowtie2)')
    cmd_vars = {
     'i':str(i),
     'sample_name':sample_name,
     'path_pipe':paths['pipe'],
     'path_o':paths['o'],
     'threads':threads}
    cmds = [
     'bowtie2-build -q {path_o}/remap/{i}.contig.fasta {path_o}/remap/{i}.contig &> /dev/null',
     'bowtie2 -x {path_o}/remap/{i}.contig -S {path_o}/remap/{i}.sam --no-unal --threads {threads} '
     '--very-sensitive-local '
     '-1 {path_o}/merge/{i}.{sample_name}.raw.r1.fastq '
     '-2 {path_o}/merge/{i}.{sample_name}.raw.r2.fastq '
     '2> {path_o}/remap/{i}.{sample_name}.bt2.stats',
     'grep -v XS:i: {path_o}/remap/{i}.sam > {path_o}/remap/{i}.uniq.sam',
     'samtools view -bS {path_o}/remap/{i}.uniq.sam | samtools sort - {path_o}/remap/{i}.uniq',
     'samtools index {path_o}/remap/{i}.uniq.bam',
     'samtools mpileup -d 1000 -f {path_o}/remap/{i}.contig.fasta {path_o}/remap/{i}.uniq.bam '
     '2> /dev/null > {path_o}/remap/{i}.uniq.pile',
     'samtools mpileup -ud 1000 -f {path_o}/remap/{i}.contig.fasta {path_o}/remap/{i}.uniq.bam '
     '2> /dev/null | bcftools call -c | vcfutils.pl vcf2fq '
     '| seqtk seq -a - > {path_o}/remap/{i}.denovo.consensus.fasta']
    cmds = [cmd.format(**cmd_vars) for cmd in cmds]
    for cmd in cmds:
        logger.info(cmd)
        cmd_run = run(cmd)
        logger.info(cmd_run.stdout)
        cmd_prefix = cmd.split(' ')[0]
        print('\tDone (' + cmd_prefix + ')') if cmd_run.returncode == 0 else sys.exit('ERR_REMAP')
    with open('{path_o}/remap/{i}.{sample_name}.bt2.stats'.format(**cmd_vars), 'r') as bt2_stats:
        bt2_count = float(bt2_stats.read().partition('% overall')[0].split('\n')[-1].strip())/100
    logger.info('Proportion of reads mapped to assembly: {}'.format(bt2_count))
    return bt2_count


def assess_assembly_coverage(best_contig_len, paths, i=1):
    print('Identifying low coverage regions...')
    depths = {}
    max_coverage = 0
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
   
    map_stats = {
     'bases_covered': bases_covered,
     'best_contig_len': best_contig_len,
     'prop_coverage': prop_coverage,
     'max_coverage': max_coverage,
     'uncovered_regions': uncovered_regions,
     'len_uncovered_regions': len(uncovered_regions),
     'largest_uncovered_region': largest_uncovered_region
    }
    return map_stats


def evaluate_assemblies(reference, target_genome_len, sample_name, paths, threads, i=1):
    '''
    Execute QUAST on all generated assemblies
    '''
    print('Comparing assemblies...')
    asm_dirs = (
     [paths['o'] + '/asm/' + dir + '/contigs.fasta' for dir in
     filter(lambda d: d.startswith(str(i)), os.listdir(paths['o'] + '/asm'))])
    cmd_vars = {
     'i':str(i),
     'asm_dirs':' '.join(asm_dirs),
     'ref_len':target_genome_len,
     'path_ref':paths['ref'],
     'sample_name':sample_name,
     'path_o':paths['o'],
     'threads':threads}
    cmd = (
     'python2 /usr/local/bin/quast.py {asm_dirs} -o {path_o}/eval/{i}.{sample_name} '
     '--threads {threads} '
     '--gene-finding'.format(**cmd_vars))
    if reference:
        cmd += ' -R {path_o}/ref/{i}.{sample_name}.ref.fasta'.format(**cmd_vars)
    if reference:
        cmd += ' -R {path_ref}'.format(**cmd_vars)
    if target_genome_len:
        cmd += ' --est-ref-size {ref_len}'.format(**cmd_vars)
    logger.info(cmd)
    cmd_run = run(cmd)
    logger.info(cmd_run.stdout)
    print('\tDone') if cmd_run.returncode == 0 else sys.exit('ERR_EVAL')


def report(start_time, end_time, paths):
    elapsed_time = end_time - start_time
    with open(paths['o'] + '/summary.txt', 'w') as report:
        report.write('wall_time\t{}'.format(elapsed_time))
        print('\tWall time: {0:.{1}f}s'.format(elapsed_time, 1))


def main(
    fwd_reads=None, rev_reads=None, reads_dir=None, out_dir='output',
    fwd_reads_sig='_F', rev_reads_sig='_R',
    norm_k_list=None, norm_cov_list=None,
    asm_k_list=None, reference_guided_asm=False, 
    reference=None, target_genome_len=10000,
    premap=False, no_remap=False,
    restart_from=False,
    threads=1):
   
    multiple_samples = True if reads_dir else False
    if not target_genome_len and reference:
        with open(reference, 'r') as reference_file:
            for record in SeqIO.parse(reference_file, 'fasta'):
                target_genome_len = len(record.seq)

    print('-' * 40)
    print('Run options...')
    print('\tMultiple input samples') if multiple_samples else print('\tSingle sample')
    print('\tPaired read signatures: \'' + fwd_reads_sig + '\', \'' + rev_reads_sig + '\'')
    print('\tTarget genome length: ' + str(target_genome_len))
    print('\t' + str(threads) + ' threads available')
   
    start_time = time.time()
    
    paths = {
    'in_dir':reads_dir,
    'in_fwd':fwd_reads,
    'in_rev':rev_reads,
    'ref':reference,
    'pipe':os.path.dirname(os.path.realpath(__file__)),
    'o':out_dir + '/run_' + str(int(time.time()))}
    
    fastqs = None
    n_reads = None
    n_reads_sample = None

    job_dirs = ['merge', 'sample', 'blast', 'ref', 'premap', 'trim', 'norm', 'asm', 'remap', 'eval']
    for dir in job_dirs:
        os.makedirs(paths['o'] + '/' + dir)

    best_asms = {}
    fastqs = list_fastqs(fwd_reads_sig, rev_reads_sig, paths)
    i = 0 # For sample IDs
    for sample_name, fastq_names in fastqs.items():
        i += 1
        import_reads(multiple_samples, sample_name, fastq_names, paths, i)
        n_reads = count_reads(sample_name, paths, i)
        if reference:
            map_to_reference(reference, sample_name, paths, threads, i)
        trim(sample_name, paths, i)
        assemble(normalise(norm_k_list, norm_cov_list, sample_name, paths, threads, i),
                 asm_k_list, reference_guided_asm, reference, sample_name, paths, threads, i)
        best_asms[sample_name] = choose_assembly(target_genome_len, sample_name, paths, threads, i)
        logger.info('best_asms: {}'.format(best_asms[sample_name]))
        prop_mapped_assembly = map_to_assembly(sample_name, paths, threads, i)
        assembly_map_stats = assess_assembly_coverage(best_asms[sample_name][2], paths, i)
        evaluate_assemblies(reference, target_genome_len, sample_name, paths, threads, i)        
    evaluate_all_assemblies(reference, target_genome_len, sample_name, paths, threads, i)
    report(start_time, time.time(), paths)


argh.dispatch_command(main)