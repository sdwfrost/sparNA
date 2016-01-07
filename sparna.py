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
# |    argh, biopython, khmer, plotly
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
from Bio import SeqUtils

import plotly.plotly as py
import plotly.graph_objs as go

import pprint


logging.basicConfig(level=logging.ERROR)
logger = logging.getLogger(__name__)


def run(cmd):
    return subprocess.run(cmd,
                          shell=True,
                          universal_newlines=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT)


def name_sample(fwd_fq):
    fwd_fq_prefix = os.path.splitext(os.path.split(fwd_fq)[1])[0]
    return fwd_fq_prefix


def import_reads(fwd_fq, rev_fq, params):
    print('Importing reads...')
    cmd = (
    'ln -s {fwd_fq} {out}/raw/{name}.f.fastq '
    '&& ln -s {rev_fq} {out}/raw/{name}.r.fastq '
    '&& interleave-reads.py {out}/raw/{name}.f.fastq {out}/raw/{name}.r.fastq '
    '> {out}/raw/{name}.fr.fastq'
    .format(**params, fwd_fq=fwd_fq, rev_fq=rev_fq))
    logger.info(cmd)
    cmd_run = run(cmd)
    # logger.info(cmd_run.stdout, cmd_run.stderr)
    print('\tDone') if cmd_run.returncode == 0 else sys.exit('ERR_IMPORT')


def trim(params):
    print('Trimming...')
    cmd = (
    'java -jar {pipe}/res/trimmomatic-0.33.jar PE '
    '{out}/raw/{name}.f.fastq '
    '{out}/raw/{name}.r.fastq '
    '{out}/trim/{name}.trim.f_pe.fastq '
    '{out}/trim/{name}.trim.f_se.fastq '
    '{out}/trim/{name}.trim.r_pe.fastq '
    '{out}/trim/{name}.trim.r_se.fastq '
    'ILLUMINACLIP:{pipe}/res/illumina_adapters.fa:2:30:10 MINLEN:30 '
    '&& cat {out}/trim/{name}.trim.f_se.fastq {out}/trim/{name}.trim.r_se.fastq '
    '> {out}/trim/{name}.trim.se.fastq '
    '&& interleave-reads.py {out}/trim/{name}.trim.f_pe.fastq {out}/trim/{name}.trim.r_pe.fastq '
    '> {out}/trim/{name}.trim.fr_pe.fastq'
    .format(**params))
    logger.info(cmd)
    cmd_run = run(cmd)
    logger.info(cmd_run.stderr)
    print('\tDone') if cmd_run.returncode == 0 else sys.exit('ERR_TRIM')


def normalise(norm_c_list, norm_k_list, params):
    print('Normalising...')
    cs = norm_c_list.split(',')
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
        '{out}/norm/{name}.norm_k{k}c{c}.pe_and_se.fastq'.format(**cmd_vars))
        cmds.append(cmd)
        print('\tNormalising norm_c={c}, norm_k={k}'.format(**cmd_vars))
        logger.info('Normalising norm_c={c}, norm_k={k}'.format(**cmd_vars))
    with multiprocessing.Pool(params['threads']) as pool:
        results = pool.map(run, cmds)
    logger.info([result.stdout + result.stdout for result in results])
    print('\tAll done') if not max([r.returncode for r in results]) else sys.exit('ERR_NORM')
    return norm_perms


def assemble(norm_perms, asm_k_list, params):
    '''
    Performs multiple assemblies and returns and OrderedDict of assembly names and paths
    **PYTHON2**
    '''
    print('Assembling...')
    asm_perms = [{'k':p['k'],'c':p['c']} for p in norm_perms]
    asm_k_list_fmt = 'k' + asm_k_list.replace(',', 'k')
    cmds_asm = []
    for asm_perm in asm_perms:
        cmd_vars = {**params,
                    'k':str(asm_perm['k']),
                    'c':str(asm_perm['c']),
                    'asm_k_list':asm_k_list,
                    'asm_k_list_fmt':asm_k_list_fmt}
        cmd_asm = (
        'python /usr/local/bin/spades.py -m 8 -t {threads} -k {asm_k_list} '
        '--pe1-1 {out}/norm/{name}.norm_k{k}c{c}.f_pe.fastq '
        '--pe1-2 {out}/norm/{name}.norm_k{k}c{c}.r_pe.fastq '
        '--s1 {out}/norm/{name}.norm_k{k}c{c}.se.fastq '
        '-o {out}/asm/{name}.norm_k{k}c{c}.asm_{asm_k_list_fmt} --careful'.format(**cmd_vars))
        cmds_asm.append(cmd_asm)
        print('\tAssembling norm_c={c}, norm_k={k}, asm_k={asm_k_list}'.format(**cmd_vars))
    with multiprocessing.Pool(params['threads']) as pool:
        results = pool.map(run, cmds_asm)
    logger.info([result.stdout for result in results])
    print('\tAll done') if not max([r.returncode for r in results]) else sys.exit('ERR_ASM')
    asms = os.listdir(params['out'] + '/asm')
    asm_paths = [params['out'] + '/asm/' + asm + '/contigs.fasta' for asm in asms]
    return OrderedDict(zip(asms, asm_paths))


def prune_assemblies(asms_paths, min_len, params):
    asms_paths_pruned = OrderedDict()
    for asm, path in asms_paths.items():
        asms_paths_pruned[asm] = path.replace('/asm/', '/asm_prune/')
        records = (r for r in SeqIO.parse(path, 'fasta') if len(r.seq) >= min_len)
        os.makedirs(asms_paths_pruned[asm].replace('/contigs.fasta', ''))
        SeqIO.write(records, asms_paths_pruned[asm], 'fasta')
    return asms_paths_pruned 


def gc_content(asms_paths):
    asms_gc = {}
    for asm, path in asms_paths.items():
        asm_gc = []
        for record in SeqIO.parse(path, 'fasta'):
            asm_gc.append(SeqUtils.GC(record.seq)/100)
        asms_gc[asm] = asm_gc
    return asms_gc


def map_to_assemblies(asms_paths, params):
    '''
    Map original reads to each assembly with Bowtie2
    Record mapping statistics
    Screen uniquely mapped reads and quantify reads mapped per contig
    Returns dict of tuples containing contig_len and n_reads_mapped for each contig
    '''
    print('Aligning to assemblies... (Bowtie2)')
    asms_coverages = {}
    for asm, asm_path in asms_paths.items():
        cmd_vars = {**params,
                    'asm':asm,
                    'asm_path':asm_path}
        cmds = [
        'bowtie2-build -q {asm_path} {out}/remap/{asm}',
        'bowtie2 -x {out}/remap/{asm} --no-unal --very-sensitive-local --threads {threads}'
        ' -1 {out}/raw/{name}.f.fastq'
        ' -2 {out}/raw/{name}.r.fastq'
        ' -S {out}/remap/{asm}.sam'
        ' 2> {out}/remap/{asm}.bt2.stats',
        'grep -v XS:i: {out}/remap/{asm}.sam > {out}/remap/{asm}.uniq.sam',
        'samtools view -bS {out}/remap/{asm}.uniq.sam'
        ' | samtools sort - {out}/remap/{asm}.uniq',
        'samtools index {out}/remap/{asm}.uniq.bam',
        'samtools idxstats {out}/remap/{asm}.uniq.bam'
        ' > {out}/remap/{asm}.uniq.bam.stats']
        cmds = [cmd.format(**cmd_vars) for cmd in cmds]
        for cmd in cmds:
            logger.info(cmd)
            cmd_run = run(cmd)
            logger.info(cmd_run.stdout)
            cmd_prefix = cmd.split(' ')[0]
            print('\tDone (' +cmd_prefix+ ')') if cmd_run.returncode == 0 else sys.exit('ERR_REMAP')
        
        with open('{out}/remap/{asm}.bt2.stats'.format(**cmd_vars), 'r') as bt2_stats:
            map_prop = float(bt2_stats.read().partition('% overall')[0].split('\n')[-1].strip())/100
        
        asm_coverages = []
        with open('{out}/remap/{asm}.uniq.bam.stats'.format(**cmd_vars), 'r') as bam_stats:
            for line in bam_stats:
                reads_mapped = int(line.strip().split('\t')[2])
                asm_coverages.append(int(reads_mapped))
        
        asms_coverages[asm] = asm_coverages
    return asms_coverages


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
            # hits_annotations = list(zip(hits)) TESTING WITHOUT SEQRECORD
            logger.info(status.text + ' ' + call.text)
            print('\t\tQuery ' + query['title'])
            # print(time.time() - start_time)
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

def fasta_blaster(fasta, max_seqs, min_len):
    '''
    NEEDS UPDATING FOR NESTED ORDEREDDICTS
    MIN_LEN NEEDS IMPLEMENTING
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
    records = OrderedDict()
    with open(fasta, 'r') as fasta_file:
        for record in SeqIO.parse(fasta_file, 'fasta'):
            if len(record.seq) >= min_len:
                records[record.id] = record.seq

    queries = [build_ebi_blast_query(title, seq) for title, seq in records.items()]
    with multiprocessing.Pool(30) as pool:
        results = pool.map(ebi_annotated_blast, queries[0:max_seqs+1])
        if len(queries) > max_seqs:
            results += zip([q['title'] for q in queries[max_seqs+1:]], [None]*len(queries[max_seqs+1:]))
    return OrderedDict(results)

def blast_assemblies(asms_paths, max_seqs, min_len):
    '''
    Returns BLAST hit information for a dict of assembly names and corresponding paths 
    '''
    print('BLASTing assemblies...')
    sample_results = OrderedDict()
    for asm_name, asm_path in asms_paths.items():
        print('\tAssembly {}'.format(asm_name))
        sample_results[asm_name] = fasta_blaster(asm_path, max_seqs=5, min_len=300)
    return sample_results

def blast_superkingdoms(blast_results):
    asms_superkingdoms = {}
    for asm, contigs in blast_results.items():
        asm_superkingdoms = []
        for contig, hits in contigs.items():
            if hits: # Hits found
                top_hit_superkingdom = hits[0][1].annotations['taxonomy'][0]
                asm_superkingdoms.append(top_hit_superkingdom)
            elif type(hits) is list: # Zero hits
                asm_superkingdoms.append(False)
            elif hits is None: # Not searched
                asm_superkingdoms.append(None)
            assert hits or type(hits) is list or hits is None
        asms_superkingdoms[asm] = asm_superkingdoms
    return asms_superkingdoms

def blast_summary(blast_results, asms_covs):
    asms_summaries = {}
    for asm, contigs in blast_results.items():
        asm_summaries = []
        for i, (contig, hits) in enumerate(contigs.items()):
            if hits: # Hits found
                description = hits[0][1].description[:40] + (hits[0][1].description[40:] and '…')
                top_hit_summary = (''
                '{2}<br>Coverage: {0} reads<br>{1}<br>Accession: {3}:{4}'
                '<br>Identity: {5}%<br>Alignment length: {6}<br>Mismatches: {7}<br>'
                'E-value: {13}'.format(asms_covs[asm][i], description, *hits[0][0]))
                asm_summaries.append(top_hit_summary)
            elif type(hits) is list: # Zero hits
                asm_summaries.append(False)
            elif hits is None: # Not searched
                asm_summaries.append(None)
            assert hits or type(hits) is list or hits is None
        asms_summaries[asm] = asm_summaries
    return asms_summaries


def plotly(asms_names, asms_stats):
    cov_max = max(sum([i for i in asms_stats['covs'].values()], []))
    cov_scale_factor = round(cov_max/5000, 1) # For bubble scaling

    traces = []
    for asm_name in asms_names:
        traces.append(
            go.Scatter(
                x=asms_stats['lens'][asm_name],
                y=asms_stats['gc'][asm_name],
                mode='lines+markers',
                name=asm_name,
                text=asms_stats['blast_summary'][asm_name],
                line=dict(shape='spline'),
                marker=dict(
                    opacity=0.5,
                    symbol='circle',
                    sizemode='area',
                    sizeref=cov_scale_factor,
                    size=asms_stats['covs'][asm_name],
                    line=dict(width=1))))

    layout = go.Layout(
        title='Assembly contig length, coverage and GC content',
        xaxis=dict(
            title='Contig length',
            gridcolor='rgb(255, 255, 255)',
            zerolinewidth=1,
            type='log',
            gridwidth=2),
        yaxis=dict(
            title='GC content',
            gridcolor='rgb(255, 255, 255)',
            zerolinewidth=1,
            gridwidth=2),
        paper_bgcolor='rgb(243, 243, 243)',
        plot_bgcolor='rgb(243, 243, 243)')

    fig = go.Figure(data=traces, layout=layout)
    return py.plot(fig)


def report(chart_url, start_time, end_time, params):
    elapsed_time = end_time - start_time
    report_content = 'wall_time\t{}\nchart_url\t{}'.format(elapsed_time, chart_url)
    with open(params['out'] + '/summary.txt', 'w') as report:
        report.write(report_content)
        print(report_content)


def main(
    fwd_fq=None, rev_fq=None, norm_c_list=None, norm_k_list=None, asm_k_list='21,33,55,77',
    untrusted_contigs=False, out_dir='', blast_max_seqs=5, min_len=100, threads=4):
    
    start_time = int(time.time())
    params = {
    'name': name_sample(fwd_fq),
    'out': out_dir + 'sparna_' + name_sample(fwd_fq),
    'pipe': os.path.dirname(os.path.realpath(__file__)),
    'threads': threads}
    for dir in ['raw', 'trim', 'norm', 'asm', 'asm_prune', 'remap', 'eval']:
        os.makedirs(params['out'] + '/' + dir)

    import_reads(fwd_fq, rev_fq, params)
    trim(params)
    norm_perms = normalise(norm_c_list, norm_k_list, params)
    asms_paths_full = assemble(norm_perms, asm_k_list, params)
    asms_paths = prune_assemblies(asms_paths_full, min_len, params)
    
    asms_names = {a: [r.id for r in SeqIO.parse(p, 'fasta')] for a, p in asms_paths.items()}
    asms_lens = {a: [int(n.split('_')[3]) for n in ns] for a, ns in asms_names.items()}
    asms_covs = map_to_assemblies(asms_paths, params)
    blast_results = blast_assemblies(asms_paths, blast_max_seqs, min_len)
    
    asms_stats = {'names': asms_names,
                 'lens': asms_lens,
                 'covs': asms_covs,
                 'blast_summary': blast_summary(blast_results, asms_covs), 
                 'blast_superkingdoms': blast_superkingdoms(blast_results),
                 'gc': gc_content(asms_paths),
                 'cpg': None}


    print(asms_stats)
    # pprint.pprint(blast_results)

    chart_url = plotly(asms_names, asms_stats)

    report(chart_url, start_time, time.time(), params)


argh.dispatch_command(main)