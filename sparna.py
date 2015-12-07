#!/usr/bin/env python3

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
import subprocess
import multiprocessing

from collections import OrderedDict

from Bio import SeqIO


logging.basicConfig(level=logging.ERROR)
logger = logging.getLogger(__name__)


def run(cmd):
    return subprocess.run(cmd,
                          shell=True,
                          universal_newlines=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT)

def name_sample(fwd_fq):
    fwd_fq_prefix = os.path.splitext(os.path.split(fwd_fq)[0])[0]
    return fwd_fq_prefix


def import_reads(fwd_fq, rev_fq, params):
    print('Importing reads...')
    cmd_run = ''
     'ln -s {fwd_fq} {out}/raw/{name}.f.fastq '
     '&& ln -s {rev_fq} {out}/raw/{name}.r.fastq '
     '&& interleave-reads.py {out}/raw/{name}.f.fastq {out}/raw/{name}.r.fastq '
     '> {out}/raw/{name}.fr.fastq'.format(**params)
    logger.info(cmd)
    cmd_run = run(cmd)
    logger.info(cmd_run.stdout)
    print('\tDone') if cmd_run.returncode == 0 else sys.exit('ERR_IMPORT')


def trim(params):
    print('Trimming...')
    cmds = [
     'java -jar {pipe}/res/trimmomatic-0.33.jar PE '
     '{out}/raw/{name}.f.fastq '
     '{out}/raw/{name}.r.fastq '
     '{out}/trim/{name}.trim.f_pe.fastq '
     '{out}/trim/{name}.trim.f_se.fastq '
     '{out}/trim/{name}.trim.r_pe.fastq '
     '{out}/trim/{name}.trim.r_se.fastq '
     'ILLUMINACLIP:{pipe}/res/illumina_adapters.fa:2:30:10 MINLEN:30',
     'cat {out}/trim/{name}.trim.f_se.fastq {out}/trim/{name}.trim.r_se.fastq '
     '> {out}/trim/{name}.trim.se.fastq '
     '&& interleave-reads.py {out}/trim/{name}.trim.f_pe.fastq {out}/trim/{name}.trim.r_pe.fastq '
     '> {out}/trim/{name}.trim.fr_pe.fastq']
    cmds = [cmd.format(**params) for cmd in cmds]
    for cmd in cmds:
        logger.info(cmd)
        cmd_run = run(cmd)
        logger.info(cmd_run.stdout)
        print('\tDone') if cmd_run.returncode == 0 else sys.exit('ERR_TRIM')


def normalise(norm_cov_list, norm_k_list, params, threads):
    print('Normalising...')
    cs = norm_cov_list.split(',')
    ks = norm_k_list.split(',')
    norm_perms = [{'k':k, 'c':c} for k in ks for c in cs]
    cmds = []
    for norm_perm in norm_perms:
        cmd_vars = {
         'k':str(norm_perm['k']),
         'c':str(norm_perm['c']),
         'pipe':params['pipe'],
         'out':params['out'],
         'name':params['name']}
        cmd = (
         'normalize-by-median.py -C {c} -k {k} -N 4 -x 1e9 -p '
         '{out}/trim/{name}.trim.fr_pe.fastq '
         '-o {out}/norm/{name}.norm_k{k}c{c}.fr_pe.fastq '
         '&& normalize-by-median.py -C {c} -k {k} -N 1 -x 1e9 '
         '{out}/trim/{name}.trim.se.fastq '
         '-o {out}/norm/{name}.norm_k{k}c{c}.se.fastq '
         '&& split-paired-reads.py '
         '-1 {out}/norm/{name}.norm_k{k}c{c}.f_pe.fastq '
         '-2 {out}/norm/{name}.norm_k{k}c{c}.r_pe.fastq '
         '{out}/norm/{name}.norm_k{k}c{c}.fr_pe.fastq '
         '&& cat {out}/norm/{name}.norm_k{k}c{c}.fr_pe.fastq '
         '{out}/norm/{name}.norm_k{k}c{c}.se.fastq > '
         '{out}/norm/{name}.norm_k{k}c{c}.pe_and_se.fastq'
         .format(**cmd_vars))
        cmds.append(cmd)
        print('\tNormalising norm_k={k},norm_c={c}'.format(**cmd_vars))
        logger.info('Normalising norm_k={k},norm_c={c}'.format(**cmd_vars))
    with multiprocessing.Pool(threads) as pool:
        results = pool.map(run, cmds)
    logger.info([result.stdout for result in results])
    print('\tDone')
    return norm_perms


def assemble(norm_perms, asm_k_list, params, threads):
    print('Assembling...')
    asm_perms = [{'k':p['k'],'c':p['c']} for p in norm_perms]
    cmds_asm = []
    for asm_perm in asm_perms:
        cmd_vars = {
         'k':str(asm_perm['k']),
         'c':str(asm_perm['c']),
         'asm_k_list':asm_k_list,
         'out':params['out'],
         'name':params['name'],
         'threads':threads}
        cmd_asm = ''
         'python2 /usr/local/bin/spades.py -m 8 -t {threads} -k {asm_k_list} '
         '--pe1-1 {out}/norm/{name}.norm_k{k}c{c}.f_pe.fastq '
         '--pe1-2 {out}/norm/{name}.norm_k{k}c{c}.r_pe.fastq '
         '--s1 {out}/norm/{name}.norm_k{k}c{c}.se.fastq '
         '-o {out}/asm/{name}.norm_k{k}c{c}.asm_k{asm_k_list} --careful'
         .format(**cmd_vars)
        cmds_asm.append(cmd_asm)
        print('\tAssembling norm_k={k},norm_c={c},asm_k={asm_k_list}'.format(**cmd_vars))
    with open(os.devnull, 'w') as devnull:
        processes = [subprocess.Popen(cmd, shell=True, stdout=devnull) for cmd in cmds_asm]
        for process in processes:
            process.wait()
            print('\tDone') if process.returncode == 0 else sys.exit('ERR_ASM')
    return os.listdir(params['out'] + '/asm')


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
        for i, record in SeqIO.parse(fasta_file, 'fasta'):
            records[record.id] = record.seq
            if seq_limit and i >= seq_limit:
                break
    # JUST FOR TESTING
    queries = [build_ebi_blast_query(title, seq) for title, seq in records.items()][0:6]
    
    with multiprocessing.Pool(30) as pool:
        results_tuple = pool.map(ebi_annotated_blast, queries)
    
    results = collections.OrderedDict(results_tuple)
    return results


def map_to_assemblies(params, threads):
    '''
    Map original reads to each assembly with Bowtie2
    Record mapping statistics
    Screen uniquely mapped reads and quantify reads mapped per contig
    Returns dict of tuples containing contig_len and n_reads_mapped for each contig
    '''
    print('Aligning to assemblies... (Bowtie2)')
    asm_names = filter(lambda d: d.startswith(str(i)), os.listdir(params['out'] + '/asm'))
    asm_names_params = {a:params['out'] + '/asm/' + a + '/contigs.fasta' for a in asm_names}
    remap_stats = {}
    for asm_name, asm_path in asm_names_params.items():
        cmd_vars = {
         'i':str(i),
         'name':name,
         'pipe':params['pipe'],
         'path_o':params['out'],
         'threads':threads,
         'asm_name':asm_name,
         'path_asm':asm_path}
        cmds = [
         'bowtie2-build -q {path_asm} {path_o}/remap/{asm_name}',
         'bowtie2 -x {path_o}/remap/{asm_name} --no-unal --very-sensitive-local --threads {threads}'
         ' -1 {path_o}/raw/{name}.f.fastq'
         ' -2 {path_o}/raw/{name}.r.fastq'
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


def report(start_time, end_time, params):
    elapsed_time = end_time - start_time
    with open(params['out'] + '/summary.txt', 'w') as report:
        report.write('wall_time\t{}'.format(elapsed_time))
        print('\tWall time: {0:.{1}f}s'.format(elapsed_time, 1))


def main(
    fwd_fq=None, rev_fq=None, norm_cov_list=None, norm_k_list=None, asm_k_list=None,
    untrusted_contigs=False, out_dir='sparna_out', threads=1):
    
    # Setup
    start_time = int(time.time())
    params = {
     'name': name_sample(),
     'out': out_dir + '/' +  name,
     'pipe': os.path.dirname(os.path.realpath(__file__))}
    for dir in ['in', 'trim', 'norm', 'asm', 'remap', 'eval']:
        os.makedirs(params['out'] + '/' + dir)

    import_reads(fwd_fq, rev_fq, params)
    trim(params)
    norm_perms = normalise(norm_cov_list, norm_k_list, params, threads)
    assemblies = assemble(norm_perms, asm_k_list, params, threads)
    # blast_results = fasta_blaster(params['out'] + '/asm/' +  best_asms[name][0] + '/contigs.fasta')
    remap_stats = map_to_assemblies(name, params, threads, i)
    report(start_time, time.time(), params)


argh.dispatch_command(main)