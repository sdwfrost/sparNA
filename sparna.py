#!/usr/bin/env python

# SparNA: A pipeline for assembling deep-sequenced viral amplicon reads
# Copyright 2015 Bede Constantinides, University of Manchester (b|at|bede|dot|im)
# Developed in collaboration with Public Health England, Colindale
# Distributed under the GNU General Public License version 3 (see LICENSE)

# TODO
# | Handle no reads mapped by BWA
# | Decide on BWA vs Bowtie2
# | remove Python2 prefix hack from spades.py and quast.py call
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
# |    bwa, bowtie2, samtools, vcftools, bcftools, bedtools, seqtk, spades, quast
# | others, bundled inside res/ directory:
# |    trimmomatic

import os
import io
import sys
import argh
import time
import pandas
import logging
import requests
import networkx
import matplotlib
import subprocess
import collections
import multiprocessing

# import plotly.plotly as py
# import plotly.graph_objs as go

from Bio import SeqIO


from collections import OrderedDict


logging.basicConfig(level=logging.ERROR)
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


def premap_to_reference(ref, sample_name, paths, threads, i=1):
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
        cmd_prefix = cmd.split(' ')[0]
        print('\tDone (' + cmd_prefix + ')') if cmd_run.returncode == 0 else sys.exit('ERR_PREMAP')
    # Return mapping stats


def trim(sample_name, paths, i=1):
    print('Trimming...')
    cmd_vars = {
     'i': str(i), 
     'path_pipe':paths['pipe'],
     'path_o':paths['o'],
     'sample_name':sample_name}
    cmds = [
     'java -jar {path_pipe}/res/trimmomatic-0.33.jar PE '
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


def assemble(norm_perms, asm_k_list, untrusted_contigs, reference, sample_name, paths, threads, i=1):
    print('Assembling...')
    if reference and untrusted_contigs:
        asm_perms = [{'k':p['k'],'c':p['c'],'uc':uc} for p in norm_perms for uc in [1, 0]]
    else:
        asm_perms = [{'k':p['k'],'c':p['c'],'uc':uc} for p in norm_perms for uc in [0]]
    cmds_asm = []
    for asm_perm in asm_perms:
        cmd_vars = {
         'i':str(i),
         'k':str(asm_perm['k']),
         'c':str(asm_perm['c']),
         'uc':str(asm_perm['uc']),
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
         '-o {path_o}/asm/{i}.{sample_name}.norm_k{k}c{c}.asm_k{asm_k_list}.uc{uc} --careful'
         .format(**cmd_vars))
        if asm_perm['uc']:
            cmd_asm += ' --untrusted-contigs {path_ref}'.format(**cmd_vars)
        cmds_asm.append(cmd_asm)
        print('\tAssembling norm_k={k},norm_c={c},asm_k={asm_k_list},uc={uc}'.format(**cmd_vars))
    with open(os.devnull, 'w') as devnull:
        processes = [subprocess.Popen(cmd, shell=True, stdout=devnull) for cmd in cmds_asm]
        for process in processes:
            process.wait()
            print('\tDone') if process.returncode == 0 else sys.exit('ERR_ASM')


def build_ebi_blast_query(title, sequence):
    '''
    Returns dict of REST params for the EBI BLAST API
    '''
    logger.info('building query')
    return { 'email': 'bede.constantinides@manchester.ac.uk',
             'program': 'blastn',
             'stype': 'dna',
             'database': 'em_rel_vrl',
             'align': 6,
             'match_scores': '1,-3',
             'gapopen': 5, 
             'gapextend': 2,
             'exp': '1e-10',
             'filter': 'T',
             'dropoff': 0,
             'scores': 5,
             'alignments': 5,
             'title': title,
             'sequence': str(sequence) }

def parse_hits(title, raw_hits):
    '''
    Returns list of tuples of BLAST hits
    [(blast, tab, output, fields), (blast, tab, output, fields)]
    '''
    hits = []
    for line in io.StringIO(raw_hits):
        if ':' in line:
            fields = [field.strip() for field in line.split('\t')]
            hit = (title, ) + tuple(fields[1].split(':')) + tuple(fields[2:])
            hits.append(hit)
    return hits

def fetch_annotation(database, accession):
    '''
    Return SeqRecord of annotation for given EMBL accession number
    '''
    query = 'http://www.ebi.ac.uk/Tools/dbfetch/dbfetch/{}/{}'.format(database, accession)
    request = requests.get(query)
    annotation = SeqIO.read(io.StringIO(request.text), 'embl')
    return annotation

def ebi_blast(query):
    '''
    Returns BLAST hits as a tuple containing a list of tuples for each hit
    ('seq', [(blast, tab, output, fields),
             (blast, tab, output, fields)])
    '''
    run_url = 'http://www.ebi.ac.uk/Tools/services/rest/ncbiblast/run/'
    status_url = 'http://www.ebi.ac.uk/Tools/services/rest/ncbiblast/status/'
    results_url = 'http://www.ebi.ac.uk/Tools/services/rest/ncbiblast/result/'

    call = requests.post(run_url, data=query)
    start_time = time.time()
    logger.info('dispatched blast jobid: ' + call.text)
    while True:
        status = requests.get(status_url + call.text)
        if status.text == 'FINISHED':
            hits_r = requests.get(results_url + call.text + '/out')
            hits = parse_hits(query['title'], hits_r.text)
            logger.info(status.text + ' ' + call.text)
            logger.info('Job completed in ' + str(time.time() - start_time))
            break
        elif time.time() - start_time > 120:
            print('blast timeout')
            logger.error('blast timeout')
            break
        elif status.text == 'RUNNING':
            time.sleep(2)
        else:
            print('status: ' + status.text)
            logger.error('status: ' + status.text)
            break
    return (query['title'], hits)

def ebi_annotated_blast(query):
    '''
    Returns BLAST hits as a tuple containing a list of tuples of hit tuples and SeqRecord annotations
    ('seq', [((blast, tab, output, fields), SeqRecord),
             ((blast, tab, output, fields), SeqRecord)])
    '''
    run_url = 'http://www.ebi.ac.uk/Tools/services/rest/ncbiblast/run/'
    status_url = 'http://www.ebi.ac.uk/Tools/services/rest/ncbiblast/status/'
    results_url = 'http://www.ebi.ac.uk/Tools/services/rest/ncbiblast/result/'

    call = requests.post(run_url, data=query)
    start_time = time.time()
    logger.info('dispatched blast jobid: ' + call.text)
    while True:
        status = requests.get(status_url + call.text)
        if status.text == 'FINISHED':
            hits_r = requests.get(results_url + call.text + '/out')
            hits = parse_hits(query['title'], hits_r.text)
            annotations_items = [hit[1] + hit[2] for hit in hits] # all there
            annotations = list(fetch_annotation(hit[1], hit[2]) for hit in hits)
            hits_annotations = list(zip(hits, annotations))
            logger.info(status.text + ' ' + call.text)
            print(time.time() - start_time)
            break
        elif time.time() - start_time > 120:
            print('blast timeout')
            logger.error('blast timeout')
            break
        elif status.text == 'RUNNING':
            time.sleep(2)
        else:
            print('status: ' + status.text)
            logger.error('status: ' + status.text)
            break
    return (query['title'], hits_annotations)

def fasta_blaster(fasta, seq_limit=0):
    '''
    CONTAINS TESTING CODE
    Returns BLAST results as an OrderedDict of ebi_blast() or ebi_annotated_blast() output
    ebi_blast():
    OrderedDict([('seq_1', [(blast, tab, output, fields),
                            (blast, tab, output, fields)],
                 ('seq_2', [(blast, tab, output, fields),
                            (blast, tab, output, fields)])])
    ebi_annotated_blast():
    OrderedDict([('seq_1', [((blast, tab, output, fields), SeqRecord),
                            ((blast, tab, output, fields), SeqRecord)],
                 ('seq_2', [((blast, tab, output, fields), SeqRecord),
                            ((blast, tab, output, fields), SeqRecord)])])
    '''
    records = collections.OrderedDict()
    with open(fasta, 'r') as fasta_file:
        for record in SeqIO.parse(fasta_file, 'fasta'):
            records[record.id] = record.seq
    # JUST FOR TESTING
    queries = [build_ebi_blast_query(title, seq) for title, seq in records.items()][0:6]
    
    with multiprocessing.Pool(30) as pool:
        results_tuple = pool.map(ebi_annotated_blast, queries)
    
    results = collections.OrderedDict(results_tuple)
    return results


def map_to_assemblies(sample_name, paths, threads, i):
    '''
    Map original reads to each assembly with Bowtie2
    Record mapping statistics
    Screen uniquely mapped reads and quantify reads mapped per contig
    Returns dict of tuples containing contig_len and n_reads_mapped for each contig
    '''
    print('Aligning to assemblies... (Bowtie2)')
    asm_names = filter(lambda d: d.startswith(str(i)), os.listdir(paths['o'] + '/asm'))
    asm_names_paths = {a:paths['o'] + '/asm/' + a + '/contigs.fasta' for a in asm_names}
    remap_stats = {}
    for asm_name, asm_path in asm_names_paths.items():
        cmd_vars = {
         'i':str(i),
         'sample_name':sample_name,
         'path_pipe':paths['pipe'],
         'path_o':paths['o'],
         'threads':threads,
         'asm_name':asm_name,
         'path_asm':asm_path}
        cmds = [
         'bowtie2-build -q {path_asm} {path_o}/remap/{asm_name}',
         'bowtie2 -x {path_o}/remap/{asm_name} --no-unal --very-sensitive-local --threads {threads}'
         ' -1 {path_o}/merge/{i}.{sample_name}.raw.r1.fastq'
         ' -2 {path_o}/merge/{i}.{sample_name}.raw.r2.fastq'
         ' -S {path_o}/remap/{asm_name}.sam'
         ' 2> {path_o}/remap/{asm_name}.bt2.stats',
         'grep -v XS:i: {path_o}/remap/{asm_name}.sam > {path_o}/remap/{asm_name}.uniq.sam',
         'samtools view -bS {path_o}/remap/{asm_name}.uniq.sam'
         ' | samtools sort - {path_o}/remap/{asm_name}.uniq',
         'samtools index {path_o}/remap/{asm_name}.uniq.bam',
         'samtools idxstats {path_o}/remap/{asm_name}.uniq.bam'
         ' > {path_o}/remap/{asm_name}.uniq.bam.stats']
        cmds = [cmd.format(**cmd_vars) for cmd in cmds]
        for cmd in cmds:
            logger.info(cmd)
            cmd_run = run(cmd)
            logger.info(cmd_run.stdout)
            cmd_prefix = cmd.split(' ')[0]
            print('\tDone (' +cmd_prefix+ ')') if cmd_run.returncode == 0 else sys.exit('ERR_REMAP')
        
        with open('{path_o}/remap/{asm_name}.bt2.stats'.format(**cmd_vars), 'r') as bt2_stats:
            map_prop = float(bt2_stats.read().partition('% overall')[0].split('\n')[-1].strip())/100
        
        contig_stats = []
        with open('{path_o}/remap/{asm_name}.uniq.bam.stats'.format(**cmd_vars), 'r') as bam_stats:
            for line in bam_stats:
                contig_name, contig_len, reads_mapped = line.strip().split('\t')[0:3]
                contig_stats.append((contig_name, int(contig_len), int(reads_mapped)))
        
        remap_stats[asm_name] = contig_stats
    print(remap_stats)
    return remap_stats


def report(start_time, end_time, paths):
    elapsed_time = end_time - start_time
    with open(paths['o'] + '/summary.txt', 'w') as report:
        report.write('wall_time\t{}'.format(elapsed_time))
        print('\tWall time: {0:.{1}f}s'.format(elapsed_time, 1))


def main(
    fwd_reads=None, rev_reads=None, reads_dir=None, out_dir='output',
    fwd_reads_sig='_F', rev_reads_sig='_R',
    norm_k_list=None, norm_cov_list=None,
    asm_k_list=None, untrusted_contigs=False, 
    reference=None, target_genome_len=10000,
    premap=False, no_remap=False,
    restart_from=False, min_contig_len=500,
    threads=1):
   
    multiple_samples = True if reads_dir else False
    if not target_genome_len and reference:
        with open(reference, 'r') as reference_file:
            for record in SeqIO.parse(reference_file, 'fasta'):
                target_genome_len = len(record.seq)
                print('woohoo')

    print('-' * 40)
    print('Run options...')
    print('\tMultiple input samples') if multiple_samples else print('\tSingle sample')
    print('\tPaired read signatures: \'' + fwd_reads_sig + '\', \'' + rev_reads_sig + '\'')
    print('\tTarget genome length: ' + str(target_genome_len))
    print('\t' + str(threads) + ' threads')
   
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

    job_dirs = ['merge', 'ref', 'premap', 'trim', 'norm', 'asm', 'subgraph', 'remap', 'eval']
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
            premap_to_reference(reference, sample_name, paths, threads, i)
        trim(sample_name, paths, i)
        assemble(normalise(norm_k_list, norm_cov_list, sample_name, paths, threads, i),
                 asm_k_list, untrusted_contigs, reference, sample_name, paths, threads, i)
        # fetch_subgraphs(paths, i)
        best_asms[sample_name] = choose_assembly(target_genome_len, sample_name, paths, threads, i)
        # blast_results = fasta_blaster(paths['o'] + '/asm/' +  best_asms[sample_name][0] + '/contigs.fasta')
        logger.info('best_asms: {}'.format(best_asms[sample_name]))
        remap_stats = map_to_assemblies(sample_name, paths, threads, i)
        # prop_mapped_assembly = map_to_longest_contig(sample_name, paths, threads, i)
        # assembly_map_stats = assess_assembly_coverage(best_asms[sample_name][2], paths, i)
        # evaluate_assemblies(reference, target_genome_len, sample_name, paths, threads, i)        
    report(start_time, time.time(), paths)


argh.dispatch_command(main)