#!/usr/bin/env python

# Author: Bede Constantinides
# Python (2.7+) pipeline for paired-end HepC assembly developed during a placement at PHE
# Please seek permission prior to use or distribution

# TODO
# | Fix final assembly evaluation 
# | Interleaved reads (ONE TRUE FORMAT)
# | Python3
# | stop using envoy
# | accept interleaved input
# | use khmer for deinterleaving (split-paired-reads.py)
# | add minimum similarity threshold for reference selection
# | report on trimming, %remapped
# | increase khmer table size
# | TESTS
# | allow reference to be supplied for QUAST and assembly
# | use roiginal filenames

# DEPENDENCIES
# | python packages:
# |    argh, biopython, envoy, khmer, matplotlib
# | others, expected inside $PATH:
# |    blast, samtools, seqtk, spades, quast, parallel (GNU)
# | others, bundled inside res/ directory:
# |    trimmomatic, fastq_deinterleave
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
    print('Identifying input... ')
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
        # print(fastqs)
    elif paths['in_fwd'] and paths['in_rev']:
        fastqs['f'].append(os.path.basename(paths['in_fwd']))
        fastqs['r'].append(os.path.basename(paths['in_rev']))
        fastqs = zip(fastqs['f'], fastqs['r'])
        fastqs = {os.path.splitext(f[0].replace(fwd_reads_sig,''))[0]: f for f in fastqs}
        # print(fastqs)
    print('\tDone') if fastqs else sys.exit('ERR_READS')
    return fastqs
    
def import_reads(multiple_samples, sample_name, fastq_names, paths, i=1):
    # need to pass in list of files with a good name
    print('Importing reads... (SAMPLE {})'.format(i))
    if multiple_samples:
        fastq_path = paths['in_dir']
    else:
        fastq_path = os.path.split(paths['in_fwd'])[0]
    cmd_vars = {
    'i':str(i),
    'sample_name':sample_name,
    'fastq_path_f':fastq_path + '/' + fastq_names[0],
    'fastq_path_r':fastq_path + '/' + fastq_names[1],
    'path_o':paths['o']}
    cmd_import = (
    'cp {fastq_path_f} {path_o}/merge/{i}.{sample_name}.raw.r1.fastq && '
    'cp {fastq_path_r} {path_o}/merge/{i}.{sample_name}.raw.r2.fastq && '
    'interleave-reads.py {path_o}/merge/{i}.{sample_name}.raw.r1.fastq '
    '{path_o}/merge/{i}.{sample_name}.raw.r2.fastq 2> /dev/null > '
    '{path_o}/merge/{i}.{sample_name}.raw.r12.fastq'
    .format(**cmd_vars))
    cmd_import = os.system(cmd_import)
    sys.exit('ERR_IMPORT') if cmd_import != 0 else None

def count_reads(sample_name, paths, i=1):
    print('Counting reads...')
    cmd_count = ('wc -l {path_o}/merge/{i}.{sample_name}.raw.r12.fastq'
    .format(i=str(i),
            path_o=paths['o'],
            sample_name=sample_name))
    cmd_count = envoy.run(cmd_count)
    n_reads = int(cmd_count.std_out.strip().split(' ')[0])/4
    print('\tDone') if cmd_count.status_code == 0 else sys.exit('ERR_COUNT')
    return n_reads

def sample_reads(n_reads, paths, i=1):
    n_reads_sample = 1e4 if n_reads >= 1e4 else n_reads
    print('Sampling ' + str(int(n_reads_sample)) + ' reads...')
    cmd_sample = (
     'cat {path_o}/merge/{i}.raw.r12.fastq | seqtk sample - '
     '{n_reads_sample} | seqtk seq -a - > {path_o}/sample/{i}.sample.fasta'
     .format(i=str(i),
             path_o=paths['o'],
             n_reads_sample=n_reads_sample))
    cmd_sample = os.system(cmd_sample)
    print('\tDone') if cmd_sample == 0 else sys.exit('ERR_SAMPLE')
    return n_reads_sample

def blast_references(paths, threads, i=1):
    print('BLASTing reference sequences...')
    if not os.path.exists(paths['pipe'] + '/res/hcv_db/db.fasta.nhr'):
        cmd_blastn_index = (
         'makeblastdb -dbtype ngcl -input_type fasta '
         '-in {path_pipe}/res/hcv_db/db.fasta -title db'
         .format(path_pipe=paths['pipe']))
        cmd_blastn_index = os.system(cmd_blastn_index)
    cmd_blastn = NcbiblastnCommandline(
        query = paths['o'] + '/sample/' + str(i) + '.sample.fasta', num_alignments = 1,
        db = paths['pipe'] + '/res/hcv_db/db.fasta', evalue = 1e-4, outfmt = 7,
        out = paths['o'] + '/blast/' + str(i) + '.blast.tsv', 
        num_threads = threads)
    cmd_blastn()
    print('\tDone')

def choose_reference(paths, i=1):
    print('Choosing reference sequence...')
    accession_freqs = {}
    with open(paths['o'] + '/blast/' + str(i) + '.blast.tsv', 'r') as blast_out:
        for line in blast_out:
            if not line.startswith('#'):
                accession = line.split('\t')[1]
                if accession in accession_freqs.keys():
                    accession_freqs[accession] += 1
                else: accession_freqs[accession] = 1
    if accession_freqs:
        ref_found = True
        top_ref_accession = max(accession_freqs, key=accession_freqs.get)
    else:
        ref_found = False
        top_ref_accession = None
        print('\tWARNING: failed to identify a similar reference sequence')
    print('\tDone')
    return ref_found, top_ref_accession

def extract_reference(top_ref_accession, paths, i=1):
    print('Extracting reference ' + top_ref_accession + '...')
    reference = ''
    with open(paths['pipe'] + '/res/hcv_db/db.fasta', 'r') as references_fa:
        inside_best_reference = False
        for line in references_fa:
            if line.startswith('>'):
                if top_ref_accession in line:
                    inside_best_reference = True
                else: inside_best_reference = False
            elif inside_best_reference:
                reference += line.strip()
    ref_len = len(reference)
    ref_path = paths['o'] + '/ref/' + str(i) + '.ref.fasta'
    with open(ref_path, 'w') as reference_fa:
        reference_fa.write('>' + top_ref_accession + '\n' + reference)
    print('\tDone')
    return ref_path, ref_len

def genotype(n_reads, n_reads_sample, paths, i=1):
    print('Genotyping...')
    prop_reads_sample = n_reads_sample/n_reads
    genotype_freqs = {}
    with open(paths['o'] + '/blast/' + str(i) + '.blast.tsv', 'r') as blast_out:
        for line in blast_out:
            if not line.startswith('#'):
                genotype = line.split('\t')[1].split('_')[1].split('.')[0]
                if genotype in genotype_freqs.keys():
                    genotype_freqs[genotype] += 1
                else: genotype_freqs[genotype] = 1
    top_genotype = max(genotype_freqs, key=genotype_freqs.get) if genotype_freqs else None
    genotype_props = {k: genotype_freqs[k]/n_reads*100 for k in genotype_freqs.keys()}
    genotype_props_pc = {k: round(genotype_freqs[k]/n_reads_sample*100, 3) for k in genotype_freqs.keys()}
    genotype_props_pc_sorted = []
    for genotype, proportion in reversed(sorted(genotype_props_pc.items(), key=lambda(k,v):(v,k))):
        record = (genotype + ': ' + str(proportion) + '% (' + str(genotype_freqs[genotype]) + ')')
        genotype_props_pc_sorted.append(record)
        print('\t' + record)
    with open(paths['o'] + '/blast/' + str(i) + '.genotypes.txt', 'w') as genotypes_file:
        for item in genotype_props_pc_sorted:
             genotypes_file.write(item + '\n')
    print('\tDone')
    return top_genotype

def map_reads(ref_path, paths, threads, i=1):
    print('Aligning... ')
    cmd_vars = {
     'i':str(i),
     'path_pipe':paths['pipe'],
     'path_o':paths['o'],
     'ref_path':ref_path,
     'threads':threads}
    cmds_map = [
     'bwa index {ref_path} &> /dev/null',
     'bwa mem -v 0 -p -t {threads} {ref_path} {path_o}/merge/{i}.raw.r12.fastq 2> /dev/null > {path_o}/map/{i}.mapped.sam',
     # '{path_pipe}/res/segemehl/segemehl.x -d {ref_path} -x {ref_path}.idx &> /dev/null',
     # '{path_pipe}/res/segemehl/segemehl.x -d {ref_path} -x {ref_path}.idx -q {path_o}/merge/{i}.raw.r12.fastq --threads {threads} -A 60 2> /dev/null > {path_o}/map/{i}.mapped.sam',
     'samtools view -bS {path_o}/map/{i}.mapped.sam | samtools sort - {path_o}/map/{i}.mapped',
     'samtools index {path_o}/map/{i}.mapped.bam',
     'samtools mpileup -d 1000 -f {ref_path} {path_o}/map/{i}.mapped.bam 2> /dev/null > {path_o}/map/{i}.mapped.pile',
     'samtools mpileup -ud 1000 -f {ref_path} {path_o}/map/{i}.mapped.bam 2> /dev/null | bcftools call -c | vcfutils.pl vcf2fastq | seqtk seq -a - | fasta_formatter -o {path_o}/map/{i}.consensus.fasta']
    for j, cmd in enumerate(cmds_map, start=1):
        cmd_map = os.system(cmd.format(**cmd_vars))
        print('\tDone (' + cmd.split(' ')[0] + ')') if cmd_map == 0 else sys.exit('ERR_MAP')

def assess_coverage(ref_len, paths, i=1):
    print('Identifying low coverage regions... ')
    min_depth = 1
    min_coverage = 0.9
    depths = {}
    uncovered_sites = []
    with open(paths['o'] + '/map/' + str(i) + '.mapped.pile', 'r') as pileup:
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

    print('\tUncovered sites: ' + str(len(uncovered_sites)))
    print('\tUncovered regions: ' + str(len(uncovered_regions)))
    if not uncovered_sites:
        print('\tAll reference bases covered!')
    elif len(uncovered_sites) < (1-min_coverage)*ref_len:
        print('\tReference coverage exceeds threshold')
    else:
        print('\tReference coverage below threshold')
    print('\tDone')

def trim(sample_name, paths, i=1):
    print('Trimming... ')
    cmd_trim = (
     'java -jar {path_pipe}/res/trimmomatic-0.32.jar PE '
     '{path_o}/merge/{i}.{sample_name}.raw.r1.fastq {path_o}/merge/{i}.{sample_name}.raw.r2.fastq '
     '{path_o}/trim/{i}.{sample_name}.trim.r1_pe.fastq {path_o}/trim/{i}.{sample_name}.trim.r1_se.fastq '
     '{path_o}/trim/{i}.{sample_name}.trim.r2_pe.fastq {path_o}/trim/{i}.{sample_name}.trim.r2_se.fastq '
     'ILLUMINACLIP:{path_pipe}/res/illumina_adapters.fa:2:30:10 MINLEN:25'
     .format(i=str(i),
             path_pipe=paths['pipe'],
             path_o=paths['o'],
             sample_name=sample_name))
    cmd_trim_pp = (
     'cat {path_o}/trim/{i}.{sample_name}.trim.r1_se.fastq {path_o}/trim/{i}.{sample_name}.trim.r2_se.fastq > '
     '{path_o}/trim/{i}.{sample_name}.trim.se.fastq '
     '&& interleave-reads.py {path_o}/trim/{i}.{sample_name}.trim.r1_pe.fastq '
     '{path_o}/trim/{i}.{sample_name}.trim.r2_pe.fastq 2> /dev/null > {path_o}/trim/{i}.{sample_name}.trim.r12_pe.fastq'
     .format(i=str(i), 
             path_pipe=paths['pipe'],
             path_o=paths['o'],
             sample_name=sample_name))
    cmd_trim = envoy.run(cmd_trim)
    cmd_trim_stats = ''.join(cmd_trim.std_err).split('\n')[25]
    cmd_trim_pp = os.system(cmd_trim_pp)
    print('\tDone') if cmd_trim.status_code == 0 and cmd_trim_pp == 0 else sys.exit('ERR_TRIM')

def normalise(norm_k_list, norm_cov_list, sample_name, paths, i=1):
    print('Normalising... ')
    ks = norm_k_list.split(',')
    cs = norm_cov_list.split(',')
    norm_perms = [{'k':k, 'c':c} for k in ks for c in cs]
    cmds_norm = []
    for norm_perm in norm_perms:
        cmd_vars = {
         'i':str(i),
         'k':str(norm_perm['k']),
         'c':str(norm_perm['c']),
         'path_pipe':paths['pipe'],
         'path_o':paths['o'],
         'sample_name':sample_name}
        cmd_norm = (
         'normalize-by-median.py -C {c} -k {k} -N 1 -x 1e9 -p '
         '{path_o}/trim/{i}.{sample_name}.trim.r12_pe.fastq '
         '-o {path_o}/norm/{i}.{sample_name}.norm_k{k}c{c}.r12_pe.fastq &> /dev/null '
         '&& normalize-by-median.py -C {c} -k {k} -N 1 -x 1e9 '
         '{path_o}/trim/{i}.{sample_name}.trim.se.fastq '
         '-o {path_o}/norm/{i}.{sample_name}.norm_k{k}c{c}.se.fastq &> /dev/null '
         '&& {path_pipe}/res/fastq_deinterleave '
         '{path_o}/norm/{i}.{sample_name}.norm_k{k}c{c}.r12_pe.fastq '
         '{path_o}/norm/{i}.{sample_name}.norm_k{k}c{c}.r1_pe.fastq '
         '{path_o}/norm/{i}.{sample_name}.norm_k{k}c{c}.r2_pe.fastq '
         '&& cat {path_o}/norm/{i}.{sample_name}.norm_k{k}c{c}.r12_pe.fastq '
         '{path_o}/norm/{i}.{sample_name}.norm_k{k}c{c}.se.fastq > '
         '{path_o}/norm/{i}.{sample_name}.norm_k{k}c{c}.pe_and_se.fastq'
         .format(**cmd_vars))
        cmds_norm.append(cmd_norm)
        print('\tNormalising norm_k={k},norm_c={c}'.format(**cmd_vars))
    with open(os.devnull, 'w') as devnull:
        processes = [subprocess.Popen(cmd, shell=True, stdout=devnull) for cmd in cmds_norm]
        for process in processes:
            process.wait()
            print('\tDone') if process.returncode == 0 else sys.exit('ERR_NORM')
    return norm_perms

def assemble(norm_perms, asm_k_list, asm_using_ref, ref_found, sample_name, paths, threads, i=1):
    print('Assembling...')
    if ref_found and asm_using_ref:
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
         'sample_name':sample_name,
         'threads':threads}
        cmd_asm = (
         'spades.py -m 8 -t {threads} -k {asm_k_list} '
         '--pe1-1 {path_o}/norm/{i}.{sample_name}.norm_k{k}c{c}.r1_pe.fastq '
         '--pe1-2 {path_o}/norm/{i}.{sample_name}.norm_k{k}c{c}.r2_pe.fastq '
         '--s1 {path_o}/norm/{i}.{sample_name}.norm_k{k}c{c}.se.fastq '
         '-o {path_o}/asm/{i}.{sample_name}.norm_k{k}c{c}.asm_k{asm_k_list}.rg{rg} --careful'
         .format(**cmd_vars))
        if asm_perm['rg']:
            cmd_asm += ' --untrusted-contigs ' + paths['o'] + '/ref/' + str(i) + '.ref.fasta'
        cmds_asm.append(cmd_asm)
        print('\tAssembling norm_k={k},norm_c={c},asm_k={asm_k_list},rg={rg}'.format(**cmd_vars))
    with open(os.devnull, 'w') as devnull:
        processes = [subprocess.Popen(cmd, shell=True, stdout=devnull) for cmd in cmds_asm]
        for process in processes:
            process.wait()
            print('\tDone') if process.returncode == 0 else sys.exit('ERR_ASM')

# First arg previosly ref_found
def evaluate_assembly(reference, est_ref_size, paths, threads, i=1):
    print('Comparing assemblies... ')
    asm_dirs = (
     [paths['o'] + '/asm/' + dir + '/contigs.fasta' for dir in
     filter(lambda d: d.startswith(str(i)), os.listdir(paths['o'] + '/asm'))])
    # print(asm_dirs)
    cmd_vars = {
    'i':str(i),
    'asm_dirs':' '.join(asm_dirs),
    'ref_len':est_ref_size,
    'path_ref':paths['ref'],
    'path_o':paths['o'],
    'threads':threads}
    cmd_eval = ('quast.py {asm_dirs} -o {path_o}/eval/{i} --threads {threads}'.format(**cmd_vars))
    # if ref_found:
    #     cmd_eval += ' -R {path_o}/ref/{i}.{sample_name}.ref.fasta'.format(**cmd_vars)
    if reference:
        print('REFERENCE')
        cmd_eval += ' -R {path_ref}'.format(**cmd_vars)
    if est_ref_size:
        cmd_eval += ' --est-ref-size {ref_len}'.format(**cmd_vars)
    cmd_eval += ' &> /dev/null'
    print(cmd_eval)
    cmd_eval = os.system(cmd_eval)
    print('\tDone') if cmd_eval == 0 else sys.exit('ERR_EVAL')

def remap_reads():
    pass

def evaluate_assemblies(paths, i):
    os.makedirs(paths['o'] + '/eval/summary/')
    cmd_vars = {
     'i':str(i),
     'path_o':paths['o']}
    cmd_report = (
     'cp {path_o}/eval/{i}/report.html {path_o}/eval/{i}/transposed_report.tsv '
     '{path_o}/eval/summary/ && cp -R {path_o}/eval/{i}/report_html_aux '
     '{path_o}/eval/summary/'.format(**cmd_vars))
    cmd_report = os.system(cmd_report)
    print('Reporting...')
    print('\tQUAST report: ' + paths['o'] + '/eval/summary/')

def report(start_time, end_time, paths):
    elapsed_time = end_time - start_time
    with open(paths['o'] + '/summary.txt', 'w') as report:
        report.write('wall_time\t{}'.format(elapsed_time))
        print('\tWall time: {0:.{1}f}s'.format(elapsed_time, 1))

def main(fwd_reads=None, rev_reads=None, reads_dir=None, out_dir='output', fwd_reads_sig='_F',
    rev_reads_sig='_R', norm_k_list=None, norm_cov_list=None, asm_k_list=None,
    reference=None, asm_using_ref=False, est_ref_size=None,
    map_before_asm=False, map_after_asm=False,
    hcv=False, threads=1):
   
    multiple_samples = True if reads_dir else False

    print('Run options...')
    print('\tMultiple input samples') if multiple_samples else print('\tSingle sample')
    print('\tPaired read signatures: \'' + fwd_reads_sig + '\', \'' + rev_reads_sig + '\'')
    print('\tVirus agnostic') if not hcv else print('Using HCV-specific features')
    print('\t' + str(threads) + ' threads')
   
    start_time = time.time()
    paths = {
    'in_dir':reads_dir,
    'in_fwd':fwd_reads,
    'in_rev':rev_reads,
    'ref':reference,
    'pipe':os.path.dirname(os.path.realpath(__file__)),
    'o':out_dir + '/run_' + str(int(time.time()))}
    state = {
    'fastqs':None,
    'n_reads':None,
    'n_reads_sample':None,
    'ref_found':None,
    'top_ref_accession':None,
    'top_genotype':None,
    'ref_path':None,
    'ref_len':None}

    job_dirs = ['merge', 'sample', 'blast', 'ref', 'map', 'trim', 'norm', 'asm', 'remap', 'eval']
    for dir in job_dirs:
        os.makedirs(paths['o'] + '/' + dir)

    state['fastqs'] = list_fastqs(fwd_reads_sig, rev_reads_sig, paths)
    i = 0 # counter for file naming
    for sample_name, fastq_names in state['fastqs'].items():
        i += 1
        import_reads(multiple_samples, sample_name, fastq_names, paths, i)
        state['n_reads'] = count_reads(sample_name, paths, i)
        if hcv:
            state['n_reads_sample'] = sample_reads(state['n_reads'], paths, i)
            blast_references(paths, threads, i)
            state['ref_found'], state['top_ref_accession'] = choose_reference(paths, i)
            if state['ref_found']:
                state['ref_path'], state['ref_len'] = extract_reference(state['top_ref_accession'], paths, i)
                state['top_genotype'] = genotype(state['n_reads'], state['n_reads_sample'], paths, i)
                map_reads(state['ref_path'], paths, threads, i)
                assess_coverage(state['ref_len'], paths, i)
        trim(sample_name, paths, i)
        assemble(normalise(norm_k_list, norm_cov_list, sample_name, paths, i), asm_k_list, asm_using_ref,
                 state['ref_found'], sample_name, paths, threads, i)
        evaluate_assembly(reference, paths, threads, i)        
    evaluate_assemblies(paths, i) if multiple_samples else None
    report(start_time, time.time(), paths)

argh.dispatch_command(main)