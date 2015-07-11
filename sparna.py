#!/usr/bin/env python

# SparNA: A pipeline for assembling deep-sequenced viral amplicon reads.
# Copyright 2015 Bede Constantinides. Distributed under the GNU General Public License version 3

# DEPENDENCIES
# | python packages:
# |    argh, biopython, envoy, khmer, matplotlib
# | others, expected inside $PATH:
# |    bwa, bowtie2, blast, samtools, vcftools, bcftools, bedtools, seqtk, spades, quast
# | others, bundled inside res/ directory:
# |    trimmomatic
# | others, bundled and requiring compilation: segemehl (optional)

# TODO
# | DOCUMENTATION
#  | Code comments
#  | Argument parsing 
#  | --help page
# | IMPLEMENTATION
#  | Improve error handling 
#  | Replace envoy.send() and os.system() calls with vanilla subprocess
#  | Migrate to BWA only for mapping... Base reports on samtools flagsts?
#  | Consistently use of r12 / fr pairing terminologuy
#  | Support for gzipped, interleaved input (ONE TRUE FORMAT)
#  | Use khmer's interleave-reads.py once bug in pair handling is fixed in v2.0
#  | Tests
#  | Py3k
#  | PEP8
# | REPORTING
#  | Send only best contigs per sample to QUAST for final eval step
#  | Report on trimming, %remapped to ref vs assembled contig
#  | Send coverage stats to file
# | HCV-SPECIFIC
#  | Improve genotyping rigour, use genome guidance
# | PERFORMANCE
#  | More pipelining to reduce disk I/O associated with mapping

# USAGE (single run): ./sparna.py --fwd-reads reads_F.fastq --rev-reads reads_R.fastq --fwd-reads-sig _F --rev-reads-sig _R --norm-k-list 21,31 --norm-c-list 1,5,20 --asm-k-list 21,33,55,77 --reference ref.fasta --out-dir /path/to/output --threads 12 
# USAGE (batch): ./sparna.py --reads-dir dir/full/of/reads/ --fwd-reads-sig _1 --rev-reads-sig _2 --norm-k-list 21,31 --norm-cov-list 1,2,5,10 --asm-k-list 21,33,55,77 --reference ref.fasta --out-dir /path/to/output --threads 12 

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
    cmd_import = (
    'cp {fastq_path_f} {path_o}/merge/{i}.{sample_name}.raw.r1.fastq && '
    'cp {fastq_path_r} {path_o}/merge/{i}.{sample_name}.raw.r2.fastq && '
    # 'interleave-reads.py {path_o}/merge/{i}.{sample_name}.raw.r1.fastq '
    '{path_pipe}/res/interleave.py {path_o}/merge/{i}.{sample_name}.raw.r1.fastq '
    '{path_o}/merge/{i}.{sample_name}.raw.r2.fastq > '
    '{path_o}/merge/{i}.{sample_name}.raw.r12.fastq'
    .format(**cmd_vars))
    cmd_import = os.system(cmd_import)
    print('\tDone') if cmd_import == 0 else sys.exit('ERR_IMPORT')


def hcv_count_reads(sample_name, paths, i=1):
    print('Counting reads...')
    cmd_count = ('wc -l {path_o}/merge/{i}.{sample_name}.raw.r12.fastq'
    .format(i=str(i),
            path_o=paths['o'],
            sample_name=sample_name))
    cmd_count = envoy.run(cmd_count)
    n_reads = int(cmd_count.std_out.strip().split(' ')[0])/4
    print('\tDone') if cmd_count.status_code == 0 else sys.exit('ERR_COUNT')
    return n_reads


def hcv_sample_reads(n_reads, paths, i=1):
    n_reads_sample = 1e4 if n_reads >= 1e4 else n_reads
    print('Sampling ' + str(int(n_reads_sample)) + ' reads...')
    cmd_sample = (
     'cat {path_o}/merge/{i}.{sample_name}.raw.r12.fastq | seqtk sample - '
     '{n_reads_sample} | seqtk seq -a - > {path_o}/sample/{i}.sample.fasta'
     .format(i=str(i),
             path_o=paths['o'],
             n_reads_sample=n_reads_sample))
    cmd_sample = os.system(cmd_sample)
    print('\tDone') if cmd_sample == 0 else sys.exit('ERR_SAMPLE')
    return n_reads_sample


def hcv_blast_references(paths, threads, i=1):
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


def hcv_choose_reference(paths, i=1):
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
        hcv_ref_found = True
        top_ref_accession = max(accession_freqs, key=accession_freqs.get)
    else:
        hcv_ref_found = False
        top_ref_accession = None
        print('\tWARNING: failed to identify a similar reference sequence')
    print('\tDone')
    return hcv_ref_found, top_ref_accession


def hcv_extract_reference(top_ref_accession, paths, i=1):
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


def hcv_genotype(n_reads, n_reads_sample, paths, i=1):
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


def hcv_map_reads(ref_path, paths, threads, i=1):
    print('Aligning...')
    cmd_vars = {
    'i':str(i),
    'sample_name':sample_name,
    'path_pipe':paths['pipe'],
    'path_o':paths['o'],
    'ref_path':ref_path,
    'threads':threads}
    cmds_map = [
    'bwa index {ref_path} &> /dev/null',
    'bwa mem -v 0 -p -t {threads} {ref_path} {path_o}/merge/{i}.raw.r12.fastq 2> /dev/null > {path_o}/map/{i}.mapped.sam',
    # '{path_pipe}/res/segemehl/segemehl.x -d {ref_path} -x {ref_path}.idx &> /dev/null',
    # '{path_pipe}/res/segemehl/segemehl.x -d {ref_path} -x {ref_path}.idx -q {path_o}/merge/{i}.{sample_name}.raw.r12.fastq --threads {threads} -A 60 2> /dev/null > {path_o}/map/{i}.mapped.sam',
    'samtools view -bS {path_o}/map/{i}.mapped.sam | samtools sort - {path_o}/map/{i}.mapped',
    'samtools index {path_o}/map/{i}.mapped.bam',
    'samtools mpileup -d 1000 -f {ref_path} {path_o}/map/{i}.mapped.bam 2> /dev/null > {path_o}/map/{i}.mapped.pile',
    'samtools mpileup -ud 1000 -f {ref_path} {path_o}/map/{i}.mapped.bam 2> /dev/null | bcftools call -c | vcfutils.pl vcf2fq | seqtk seq -a - > {path_o}/map/{i}.consensus.fasta']
    for j, cmd in enumerate(cmds_map, start=1):
        cmd_map = os.system(cmd.format(**cmd_vars))
        print('\tDone (' + cmd.split(' ')[0] + ')') if cmd_map == 0 else sys.exit('ERR_MAP') # sample_name should be in here


def hcv_assess_coverage(ref_len, paths, i=1):
    print('Identifying low coverage regions...')
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
    print('\tLargest uncovered region: ' + str(largest_uncovered_region))

    # if not uncovered_sites:
    #     print('\tAll reference bases covered!')
    # elif len(uncovered_sites) < (1-min_coverage)*ref_len:
    #     print('\tReference coverage exceeds threshold')
    # else:
    #     print('\tReference coverage below threshold')
    # print('\tDone')


def map_reads(use_segemehl, ref, sample_name, paths, threads, i=1):
    print('Aligning... (Segemehl)') if use_segemehl else print('Aligning... (BWA)')
    cmd_vars = {
     'i':str(i),
     'sample_name':sample_name,
     'path_segemehl':paths['segemehl'],
     'path_pipe':paths['pipe'],
     'path_o':paths['o'],
     'ref':ref,
     'threads':threads}
    cmds_map = ['cp {ref} {path_o}/ref/ref.fasta']
    if use_segemehl and os.system('{path_segemehl} &> /dev/null'.format(**cmd_vars)):
        cmds_map += [
         '{path_segemehl} -d {ref} -x {ref}.idx &> /dev/null',
         '{path_segemehl} -d {ref} -x {ref}.idx -q {path_o}/merge/{i}.{sample_name}.raw.r12.fastq '
         '--threads {threads} -A 60 2> /dev/null > {path_o}/map/{i}.{sample_name}.mapped.sam ']
    else:
        cmds_map += [
         'bwa index {ref} &> /dev/null',
         'bwa mem -v 0 -p -t {threads} {ref} {path_o}/merge/{i}.{sample_name}.raw.r12.fastq '
         '2> /dev/null > {path_o}/map/{i}.{sample_name}.mapped.sam ']
    cmds_map += [
     'samtools view -bS {path_o}/map/{i}.{sample_name}.mapped.sam | '
     'samtools sort - {path_o}/map/{i}.{sample_name}.mapped',
     'samtools index {path_o}/map/{i}.{sample_name}.mapped.bam',
     'samtools mpileup -d 1000 -f {ref} {path_o}/map/{i}.{sample_name}.mapped.bam 2> /dev/null '
     '> {path_o}/map/{i}.{sample_name}.mapped.pile',
     'samtools mpileup -ud 1000 -f {ref} {path_o}/map/{i}.{sample_name}.mapped.bam 2> /dev/null '
     '| bcftools call -c | vcfutils.pl vcf2fq | seqtk seq -a - > {path_o}/map/{i}.{sample_name}.consensus.fasta']
    for j, cmd in enumerate(cmds_map, start=1):
        cmd = cmd.format(**cmd_vars)
        cmd_map = os.system(cmd)
        print('\tDone (' + cmd.split(' ')[0].format(**cmd_vars) + ')') if cmd_map == 0 else sys.exit('ERR_MAP')


def trim(sample_name, paths, i=1):
    print('Trimming...')
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
     # '&& interleave-reads.py {path_o}/trim/{i}.{sample_name}.trim.r1_pe.fastq '
     # '{path_o}/trim/{i}.{sample_name}.trim.r2_pe.fastq 2> /dev/null > {path_o}/trim/{i}.{sample_name}.trim.r12_pe.fastq'
     '&& {path_pipe}/res/interleave.py {path_o}/trim/{i}.{sample_name}.trim.r1_pe.fastq '
     '{path_o}/trim/{i}.{sample_name}.trim.r2_pe.fastq > {path_o}/trim/{i}.{sample_name}.trim.r12_pe.fastq'
     .format(i=str(i), 
             path_pipe=paths['pipe'],
             path_o=paths['o'],
             sample_name=sample_name))
    cmd_trim = envoy.run(cmd_trim)
    cmd_trim_stats = ''.join(cmd_trim.std_err).split('\n')[25]
    print(cmd_trim.std_out, cmd_trim.std_err)
    cmd_trim_pp = os.system(cmd_trim_pp)
    print('\tDone') if cmd_trim.status_code == 0 and cmd_trim_pp == 0 else sys.exit('ERR_TRIM')


def normalise(norm_k_list, norm_cov_list, sample_name, paths, i=1):
    print('Normalising...')
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
         'normalize-by-median.py -C {c} -k {k} -N 4 -x 1e9 -p '
         '{path_o}/trim/{i}.{sample_name}.trim.r12_pe.fastq '
         '-o {path_o}/norm/{i}.{sample_name}.norm_k{k}c{c}.r12_pe.fastq &> /dev/null '
         '&& normalize-by-median.py -C {c} -k {k} -N 1 -x 1e9 '
         '{path_o}/trim/{i}.{sample_name}.trim.se.fastq '
         '-o {path_o}/norm/{i}.{sample_name}.norm_k{k}c{c}.se.fastq &> /dev/null '
         '&& split-paired-reads.py '
         '-1 {path_o}/norm/{i}.{sample_name}.norm_k{k}c{c}.r1_pe.fastq '
         '-2 {path_o}/norm/{i}.{sample_name}.norm_k{k}c{c}.r2_pe.fastq '
         '{path_o}/norm/{i}.{sample_name}.norm_k{k}c{c}.r12_pe.fastq &> /dev/null '
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
         'spades.py -m 8 -t {threads} -k {asm_k_list} '
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


def choose_assembly(est_ref_len, sample_name, paths, threads, i=1):
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
            longest_contig_len = None
            for record in SeqIO.parse(contigs_file, 'fasta'):
                if len(record.seq) > longest_contig_len:
                    longest_contig_len = len(record.seq) 
                    longest_contig_name = record.id
        longest_contigs[asm_name] = (longest_contig_name, longest_contig_len)
        # print(longest_contigs)
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
     'bowtie2-build -q {path_o}/remap/{i}.contig.fasta {path_o}/remap/{i}.contig &> /dev/null',
     'bowtie2 -x {path_o}/remap/{i}.contig -S {path_o}/remap/{i}.sam --no-unal --threads 12 --local '
     '-1 {path_o}/merge/{i}.{sample_name}.raw.r1.fastq '
     '-2 {path_o}/merge/{i}.{sample_name}.raw.r2.fastq '
     '2> {path_o}/remap/{i}.{sample_name}.bt2.stats',
     'echo -e "\t $(tail -n 1 {path_o}/remap/{i}.{sample_name}.bt2.stats)"',
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


def evaluate_assemblies(reference, est_ref_len, sample_name, paths, threads, i=1):
    print('Comparing assemblies...')
    asm_dirs = (
     [paths['o'] + '/asm/' + dir + '/contigs.fasta' for dir in
     filter(lambda d: d.startswith(str(i)), os.listdir(paths['o'] + '/asm'))])
    # print(asm_dirs)
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
     '--gene-finding'.format(**cmd_vars))
    # if reference:
    #     cmd_eval += ' -R {path_o}/ref/{i}.{sample_name}.ref.fasta'.format(**cmd_vars)
    if reference:
        cmd_eval += ' -R {path_ref}'.format(**cmd_vars)
    if est_ref_len:
        cmd_eval += ' --est-ref-size {ref_len}'.format(**cmd_vars)
    cmd_eval += ' --min-contig 50' # Just for testing
    cmd_eval += ' &> /dev/null'
    cmd_eval = os.system(cmd_eval)
    print('\tDone') if cmd_eval == 0 else sys.exit('ERR_EVAL')


def evaluate_all_assemblies(reference, est_ref_len, sample_name, paths, threads, i=1):
    print('Comparing all assemblies...')
    asm_dirs = (
    [paths['o'] + '/asm/' + dir + '/contigs.fasta' for dir in
    filter(lambda d: d[0].isdigit(), os.listdir(paths['o'] + '/asm'))])
    # print(asm_dirs)
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


def main(
    fwd_reads=None, rev_reads=None, reads_dir=None, out_dir='output',
    fwd_reads_sig='_F', rev_reads_sig='_R',
    norm_k_list=None, norm_cov_list=None,
    asm_k_list=None, reference_guided_asm=False, 
    reference=None, est_ref_len=None,
    use_segemehl=False,
    no_map=False, no_remap=False,
    hcv=False,
    threads=1):
   
    multiple_samples = True if reads_dir else False
    if not est_ref_len and reference:
        with open(reference, 'r') as reference_file:
            for record in SeqIO.parse(reference_file, 'fasta'):
                est_ref_len = len(record.seq)

    print('-' * 40)
    print('Run options...')
    print('\tMultiple input samples') if multiple_samples else print('\tSingle sample')
    print('\tPaired read signatures: \'' + fwd_reads_sig + '\', \'' + rev_reads_sig + '\'')
    print('\tVirus agnostic') if not hcv else print('Using HCV-specific features')
    print('\tUsing Segemehl') if use_segemehl else print('\tUsing BWA')
    print('\t' + str(threads) + ' threads available')
   
    start_time = time.time()
    
    paths = {
    'in_dir':reads_dir,
    'in_fwd':fwd_reads,
    'in_rev':rev_reads,
    'ref':reference,
    'pipe':os.path.dirname(os.path.realpath(__file__)),
    'o':out_dir + '/run_' + str(int(time.time())),
    'segemehl':os.path.dirname(os.path.realpath(__file__)) + '/res/segemehl/segemehl.x'}
    
    fastqs = None
    n_reads = None
    n_reads_sample = None
    hcv_ref_found = None
    hcv_ref_len = None
    hcv_ref_path = None
    hcv_top_accession = None
    hcv_top_genotype = None

    job_dirs = ['merge', 'sample', 'blast', 'ref', 'map', 'trim', 'norm', 'asm', 'remap', 'eval']
    for dir in job_dirs:
        os.makedirs(paths['o'] + '/' + dir)

    best_asms = {}
    fastqs = list_fastqs(fwd_reads_sig, rev_reads_sig, paths)
    i = 0 # For sample IDs
    for sample_name, fastq_names in fastqs.items():
        i += 1
        import_reads(multiple_samples, sample_name, fastq_names, paths, i)
        n_reads = hcv_count_reads(sample_name, paths, i)
        if hcv:    
            n_reads_sample = hcv_sample_reads(n_reads, paths, i)
            hcv_blast_references(paths, threads, i)
            hcv_ref_found, hcv_top_accession = hcv_choose_reference(paths, i)
            if hcv_ref_found:
                ref_path, ref_len = hcv_extract_reference(hcv_top_accession, paths, i)
                hcv_top_genotype = hcv_genotype(n_reads, n_reads_sample, paths, i)
                hcv_map_reads(ref_path, paths, threads, i)
                hcv_assess_coverage(ref_len, paths, i)
        elif reference:
            map_reads(use_segemehl, reference, sample_name, paths, threads, i)
        trim(sample_name, paths, i)
        assemble(normalise(norm_k_list, norm_cov_list, sample_name, paths, i),
                 asm_k_list, reference_guided_asm, reference, sample_name, paths, threads, i)
        best_asms[sample_name] = choose_assembly(est_ref_len, sample_name, paths, threads, i)
        map_reads_to_assembly(sample_name, paths, threads, i)
        assess_remap_coverage( best_asms[sample_name][2], paths, i)
        evaluate_assemblies(reference, est_ref_len, sample_name, paths, threads, i)        
    evaluate_all_assemblies(reference, est_ref_len, sample_name, paths, threads, i)
    report(start_time, time.time(), paths)

argh.dispatch_command(main)