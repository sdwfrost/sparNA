#!/usr/bin/env python

# Author: Bede Constantinides
# Python (2.7+) pipeline for paired-end HepC assembly developed during a placement at PHE
# Please seek permission prior to use or distribution

# TODO
# | if no reads are aligned in the blast step, try blasting every single read
# | itop appeasing Insanely Bad Format - keep reads interleaved except as required
# | use khmer for deinterleaving (split-paired-reads.py)
# | sort out assess_coverage() (currently disused)
# | add minimum similarity threshold for reference selection
# | mauve/nucmer integration for when a contiguous assembly is unavailable
# | parellelise normalisation and assembly for speed
# | report on trimming, %remapped
# | increase khmer table size
# | improve command line interface
# | check exit code of every call to os.system?
# | refactor main method
# | pep8

# DEPENDENCIES
# | python packages:
# |   argh, biopython, envoy, khmer
# | others, expected inside $PATH:
# |   blast, samtools, seqtk, spades, quast
# | others, bundled inside res/ directory:
# |   trimmomatic, fastq_deinterleave
# | others, bundled and requiring compilation: segemehl

# USAGE: ./pipeline.py --threads 12 --fwd-reads-sig _F --rev-reads-sig _R --norm-k-list 31 --norm-c-list 5 --asm-k-list 33 --multiple-samples --in-dir /path/to/fastqs --out-dir /path/to/output
# Input fastq filenames should have an extension and a signature to allow identification of forward and reverse reads

# min_cov
# segemehl_sensitivity_pc
# min_depth

from __future__ import division, print_function
import os
import sys
import time
import argh
import envoy
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline

def list_fastqs(fwd_reads_sig, rev_reads_sig, paths):
   print('Identifying input... ')
   fastqs = {'f':[], 'r':[]}
   for fastq in os.listdir(paths['in']):
      if fastq.endswith('.fastq') or fastq.endswith('.fastq'):
         if fwd_reads_sig in fastq:
            fastqs['f'].append(paths['in'] + '/' + fastq)
         elif rev_reads_sig in fastq:
            fastqs['r'].append(paths['in'] + '/' + fastq)
   fastq_pairs = zip(fastqs['f'], fastqs['r'])
   fastq_pairs = {os.path.splitext(p[0].replace(fwd_reads_sig,''))[0]: p for p in fastq_pairs}
   if not fastq_pairs:
      sys.exit('Error reading paired FASTQs')
   # for fastq_pair in fastq_pairs.val:
   #  print(fastq_pair)
   return fastqs, fastq_pairs # return list of filenames, tuples of tuples of filename pairs
   
def merge_fastq_pairs(fastqs, paths, i=1):
   print('Merging reads... ')
   cmd_merge = os.system(''
   + 'cat ' + ' '.join(fastqs['f']) + ' > ' + paths['out'] + '/merge/' + str(i)
   + '.raw.r1.fastq && ' + 'cat ' + ' '.join(fastqs['r']) + ' > ' + paths['out'] + '/merge/'
   + str(i) + '.raw.r2.fastq && ' + 'interleave-reads.py ' + paths['out'] + '/merge/' + str(i)
   + '.raw.r1.fastq ' + paths['out'] + '/merge/' + str(i) + '.raw.r2.fastq 2> /dev/null > ' 
   + paths['out'] + '/merge/' + str(i) + '.raw.r12.fastq')

def copy_fastq_pair(fastq_pair, paths, i=1):
   print('Copying reads for sample ' + str(i) + '... \n\t' + fastq_pair[0] + '\n\t' + fastq_pair[1])
   cmd_copy = os.system(''
   + 'cp ' + fastq_pair[0] + ' ' + paths['out'] + '/merge/' + str(i) + '.raw.r1.fastq && '
   + 'cp ' + fastq_pair[1] + ' ' + paths['out'] + '/merge/' + str(i) + '.raw.r2.fastq && '
   + 'interleave-reads.py ' + paths['out'] + '/merge/' + str(i) + '.raw.r1.fastq ' + paths['out']
   + '/merge/' + str(i) + '.raw.r2.fastq 2> /dev/null > ' + paths['out'] + '/merge/' + str(i)
   + '.raw.r12.fastq')

def count_reads(paths, i=1):
   print('Counting reads...')
   cmd_count = envoy.run('wc -l ' + paths['out'] + '/merge/' + str(i) + '.raw.r12.fastq')
   # print(cmd_count.std_out, cmd_count.std_err)
   n_reads = int(cmd_count.std_out.replace(' ','').split('/')[0])/4
   return n_reads

def sample_reads(n_reads, paths, i=1):
   print('Sampling reads...')
   n_reads_sample = n_reads/100
   cmd_sample = os.system(''
   + 'cat ' + paths['out'] + '/merge/' + str(i) + '.raw.r12.fastq | seqtk sample - '
   + str(n_reads_sample) + ' | seqtk seq -a - > ' + paths['out'] + '/sample/' + str(i)
   + '.sample.fasta')
   return n_reads_sample

def blast_references(paths, threads, i=1):
   print('BLASTing reference sequences...')
   if not os.path.exists(paths['pipe'] + '/res/hcv_db/db.nhr'):
      cmd_blastn_index = envoy.run(''
      + 'makeblastdb -dbtype nucl -input_type fasta -in '
      + paths['pipe'] + '/res/hcv_db/db.fasta -title db')
      # print(cmd_blastn_index.std_out, cmd_blastn_index.std_err)
   cmd_blastn = NcbiblastnCommandline(
   query = paths['out'] + '/sample/' + str(i) + '.sample.fasta',
   db = paths['pipe'] + '/res/hcv_db/db.fasta',
   out = paths['out'] + '/blast/' + str(i) + '.blast.tsv',
   evalue = 1e-4,
   outfmt = 7,
   num_alignments = 1,
   num_threads = threads)
   cmd_blastn()

def choose_reference(paths, i=1):
   print('Choosing reference sequence...')
   accession_freqs = {}
   with open(paths['out'] + '/blast/' + str(i) + '.blast.tsv', 'r') as blast_out:
      for line in blast_out:
         if not line.startswith('#'):
            accession = line.split('\t')[1]
            if accession in accession_freqs.keys():
               accession_freqs[accession] += 1
            else: accession_freqs[accession] = 1
   if accession_freqs:
      reference_found = True
      top_accession = max(accession_freqs, key=accession_freqs.get)
   else:
      reference_found = False
      top_accession = None
      print('\t WARNING: failed to identify a similar reference sequence')
   return reference_found, top_accession

def extract_reference(top_accession, paths, i=1):
   print('\tExtracting ' + top_accession + '...')
   reference = ''
   with open(paths['pipe'] + '/res/hcv_db/db.fasta', 'r') as references_fa:
      inside_best_reference = False
      for line in references_fa:
         if line.startswith('>'):
            if top_accession in line:
               inside_best_reference = True
            else: inside_best_reference = False
         elif inside_best_reference:
            reference += line.strip()
   reference_len = len(reference)
   reference_fa_path = paths['out'] + '/ref/' + str(i) + '.ref.fasta'
   with open(reference_fa_path, 'w') as reference_fa:
      reference_fa.write('>' + top_accession + '\n' + reference)
   return reference_fa_path, reference_len

def genotype(n_reads, paths, i=1):
   genotype_freqs = {}
   with open(paths['out'] + '/blast/' + str(i) + '.blast.tsv', 'r') as blast_out:
      for line in blast_out:
         if not line.startswith('#'):
            genotype = line.split('\t')[1].split('_')[1].split('.')[0]
            if genotype in genotype_freqs.keys():
               genotype_freqs[genotype] += 1
            else: genotype_freqs[genotype] = 1
   top_genotype = max(genotype_freqs, key=genotype_freqs.get) if genotype_freqs else None
   genotype_props = {k: genotype_freqs[k]/n_reads*100 for k in genotype_freqs.keys()}
   genotype_props_pc = {k: round(genotype_freqs[k]/n_reads*1e4, 2) for k in genotype_freqs.keys()}
   genotype_props_pc_sorted = []
   for genotype, proportion in reversed(sorted(genotype_props_pc.items(), key=lambda(k,v):(v,k))):
      record = (genotype + ': ' + str(proportion) + '% (' + str(genotype_freqs[genotype]) + ')')
      genotype_props_pc_sorted.append(record)
      print('\t' + record)
   with open(paths['out'] + '/blast/' + str(i) + '.genotypes.txt', 'w') as genotypes_file:
      for item in genotype_props_pc_sorted:
          genotypes_file.write(item + '\n')
   return top_genotype

def map_reads(reference_fa_path, paths, threads, i=1):
   print('Aligning with Segemehl... ')
   cmds_align = [
      # paths['pipe'] + '/res/segemehl/segemehl.x -d ' + reference_fa_path + ' -x '
      # + reference_fa_path + '.idx &> /dev/null',
      paths['pipe'] + '/res/segemehl/segemehl.x -d ' + reference_fa_path + ' -x '
      + reference_fa_path + '.idx' + ' -q ' + paths['out'] + '/merge/' + str(i) + '.raw.r12.fastq'
      + ' --threads ' + str(threads) + ' -A 60 > ' + paths['out'] + '/map/' + str(i)
      + '.segemehl_mapped.sam',
      'samtools view -bS ' + paths['out'] + '/map/' + str(i) + '.segemehl_mapped.sam '
      + '| samtools sort - ' + paths['out'] + '/map/' + str(i) + '.segemehl_mapped; ',
      'samtools index ' + paths['out'] + '/map/' + str(i) + '.segemehl_mapped.bam',
      'samtools mpileup -d 1000 -f ' + reference_fa_path + ' ' + paths['out']
      + '/map/' + str(i) + '.segemehl_mapped.bam > ' + paths['out'] + '/map/' + str(i)
      + '.segemehl_mapped.pile',
      'samtools mpileup -ud 1000 -f ' + reference_fa_path + ' ' + paths['out'] + '/map/' + str(i)
      + '.segemehl_mapped.bam | ' + 'bcftools call -c | vcfutils.pl vcf2fq | seqtk seq -a - | '
      + 'fasta_formatter -o ' + paths['out'] + '/map/' + str(i) + '.consensus.fasta']
   for cmd_align in cmds_align:
      os.system(cmd_align)

def assess_coverage(reference_len, paths, i=1):
   print('Identifying low coverage regions... ')
   min_depth = 1
   min_coverage = 0.9
   depths = {}
   uncovered_sites = []
   with open(paths['out'] + '/map/' + str(i) + '.segemehl_mapped.pile', 'r') as pileup:
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
      print('All reference bases covered!')
   elif len(uncovered_sites) < (1-min_coverage)*reference_len:
      print('Reference coverage safely above threshold')
   else:
      print('Reference coverage below threshold')

def trim(paths, i=1):
   print('Trimming... ')
   cmd_trim = envoy.run(''
   + 'java -jar ' + paths['pipe'] + '/res/trimmomatic-0.32.jar PE ' + paths['out']
   + '/merge/' + str(i) + '.raw.r1.fastq ' + paths['out'] + '/merge/' + str(i) + '.raw.r2.fastq '
   + paths['out'] + '/trim/' + str(i) + '.trim.r1_pe.fastq ' + paths['out'] + '/trim/' + str(i)
   + '.trim.r1_se.fastq ' + paths['out'] + '/trim/' + str(i) + '.trim.r2_pe.fastq ' + paths['out']
   + '/trim/' + str(i) + '.trim.r2_se.fastq ILLUMINACLIP:' + paths['pipe']
   + '/res/illumina_adapters.fa:2:30:10 MINLEN:25')
   # print(cmd_trim.std_out, cmd_trim.std_err)
   os.system(''
   + 'cat ' + paths['out'] + '/trim/' + str(i) + '.trim.r1_se.fastq ' + paths['out'] + '/trim/'
   + str(i) + '.trim.r2_se.fastq > ' + paths['out'] + '/trim/' + str(i) + '.trim.se.fastq')
   os.system(''
   + 'interleave-reads.py ' + paths['out'] + '/trim/' + str(i) + '.trim.r1_pe.fastq '
   + paths['out'] + '/trim/' + str(i) + '.trim.r2_pe.fastq 2> /dev/null > '
   + paths['out'] + '/trim/' + str(i) + '.trim.r12_pe.fastq')
   cmd_trim_stats = ''.join(cmd_trim.std_err).split('\n')[25]

def normalise(norm_k_list, norm_c_list, paths, i=1):
   print('Normalising... ')
   ks = norm_k_list.split(',')
   cs = norm_c_list.split(',')
   norm_perms = [{'k':k, 'c':c} for k in ks for c in cs]
   for norm_perm in norm_perms:
      k, c = str(norm_perm['k']), str(norm_perm['c'])
      cmd_norm = os.system(''
      + 'normalize-by-median.py -C ' + c + ' -k ' + k + ' -N 1 -x 1e9 -p ' + paths['out'] + '/trim/'
      + str(i) + '.trim.r12_pe.fastq -o ' + paths['out'] + '/norm/' + str(i) + '.norm_k' + k + 'c'
      + c + '.r12_pe.fastq &> /dev/null && normalize-by-median.py -C ' + c + ' -k ' + k
      + ' -N 1 -x 1e9 ' + paths['out'] + '/trim/' + str(i) + '.trim.se.fastq ' + '-o '
      + paths['out'] + '/norm/' + str(i) + '.norm_k' + k + 'c' + c + '.se.fastq &> /dev/null && '
      + paths['pipe'] + '/res/fastq_deinterleave ' + paths['out'] + '/norm/' + str(i) + '.norm_k'
      + k + 'c' + c + '.r12_pe.fastq ' + paths['out'] + '/norm/' + str(i) + '.norm_k' + k + 'c' + c
      + '.r1_pe.fastq ' + paths['out'] + '/norm/' + str(i) + '.norm_k' + k + 'c' + c
      + '.r2_pe.fastq &> /dev/null && ' + 'cat ' + paths['out'] + '/norm/' + str(i) + '.norm_k' + k
      + 'c' + c + '.r12_pe.fastq ' + paths['out'] + '/norm/' + str(i) + '.norm_k' + k + 'c' + c
      + '.se.fastq > ' + paths['out'] + '/norm/' + str(i) + '.norm_k' + k + 'c' + c
      + '.pe_and_se.fastq')
      print('\tDone (k=' + k + ', c=' + c + ')') if cmd_norm == 0 else sys.exit('ERR_NORM')
   return norm_perms

def assemble(norm_perms, asm_k_list, asm_untrusted_contigs, reference_found, paths, threads, i=1):
   print('Assembling... ')
   if reference_found and asm_untrusted_contigs:
      asm_perms = [{'k':p['k'],'c':p['c'],'uc':uc} for p in norm_perms for uc in [1, 0]]
   else:
      asm_perms = [{'k':p['k'],'c':p['c'],'uc':uc} for p in norm_perms for uc in [0]]
   for asm_perm in asm_perms:
      k, c, uc = str(asm_perm['k']), str(asm_perm['c']), str(asm_perm['uc'])
      norm_path_prefix = paths['out'] + '/norm/' + str(i) + '.norm_k' + k + 'c' + c    
      asm_path = paths['out'] + '/asm/' + str(i) + '.norm_k' + k + 'c' + c + '.asm_k' + asm_k_list
      + '.uc' + uc
      cmd_asm = (''
      + 'spades.py -m 8 -t ' + str(threads) + ' -k ' + asm_k_list + ' '
      + '--pe1-1 ' + norm_path_prefix + '.r1_pe.fastq '
      + '--pe1-2 ' + norm_path_prefix + '.r2_pe.fastq '
      + '--s1 ' + norm_path_prefix + '.se.fastq '
      + '-o ' + asm_path + ' --careful')
      if asm_perm['uc']:
         cmd_asm += ' --untrusted-contigs ' + paths['out'] + '/ref/' + str(i) + '.ref.fasta'
      cmd_asm = envoy.run(cmd_asm)
      # print(cmd_asm.std_out, cmd_asm.std_err)
      print('\tDone (k=' + asm_k_list + ')') if cmd_asm.status_code == 0 else sys.exit('ERR_ASM')

def evaluate_assemblies(reference_found, paths, threads, i=1):
   print('Comparing assemblies... ')
   asm_dirs = [paths['out'] + '/asm/' + dir + '/contigs.fasta' for dir in os.listdir(paths['out'] + '/asm') if not dir.startswith('.')]
   eval_cmd = ('quast.py ' + ' '.join(asm_dirs) + ' -o ' + paths['out'] + '/eval/' + str(i)
   + ' --threads ' + str(threads))
   if reference_found:
      eval_cmd += ' -R ' + paths['out'] + '/ref/' + str(i) + '.ref.fasta'
   eval_cmd += ' &> /dev/null'
   os.system(eval_cmd)

def remap_reads():
   pass

def report(paths, i):
   os.makedirs(paths['out'] + '/eval/summary/')
   os.system('cp -R ' + paths['out'] + '/eval/' + str(i) + '/report.html ' + paths['out']
   + '/eval/' + str(i) + '/transposed_report.tsv ' + paths['out'] + '/eval/' + str(i)
   + '/report_html_aux ' + paths['out'] + '/eval/summary/')
   print('\tQUAST report: ' + paths['out'] + '/eval/summary/')

def main(in_dir=None, out_dir=None, fwd_reads_sig=None, rev_reads_sig=None, norm_k_list=None,
   norm_c_list=None, asm_k_list=None, asm_untrusted_contigs=False, multiple_samples=False, threads=1):
   paths = {
      'in':in_dir,
      'pipe':os.path.dirname(os.path.realpath(__file__)),
      'out':out_dir + '/run_' + str(int(time.time()))
   }
   os.makedirs(paths['out'])
   os.mkdir(paths['out'] + '/merge/')
   os.mkdir(paths['out'] + '/sample/')
   os.mkdir(paths['out'] + '/blast/')
   os.mkdir(paths['out'] + '/ref/')
   os.mkdir(paths['out'] + '/map/')
   os.mkdir(paths['out'] + '/trim/')
   os.mkdir(paths['out'] + '/norm/')
   os.mkdir(paths['out'] + '/asm/')
   os.mkdir(paths['out'] + '/remap/')
   os.mkdir(paths['out'] + '/eval/')

   fastqs, fastq_pairs = list_fastqs(fwd_reads_sig, rev_reads_sig, paths)
   if multiple_samples:
      for i, fastq_pair in enumerate(fastq_pairs, start=1):
         copy_fastq_pair(fastq_pairs[fastq_pair], paths, i)
         n_reads = count_reads(paths, i)
         sample_reads(n_reads, paths, i)
         blast_references(paths, threads, i)
         reference_found, top_accession = choose_reference(paths, i)
         if reference_found:
            reference_fa_path, reference_len = extract_reference(top_accession, paths, i)
            top_genotype = genotype(n_reads, paths, i)
            map_reads(reference_fa_path, paths, i)
            assess_coverage(reference_len, paths, i)
         trim(paths, i)
         assemble(normalise(norm_k_list, norm_c_list, paths, i), asm_k_list, asm_untrusted_contigs, reference_found, paths, threads, i)
         evaluate_assemblies(reference_found, paths, threads, i)
         remap_reads()
      report(paths, i)
   else:
      merge_fastq_pairs(fastqs. paths)
      n_reads = count_reads(paths)
      sample_reads(n_reads, paths)
      blast_references(paths, threads)
      reference_found, top_accession = choose_reference(paths)
      reference_fa_path, reference_len = extract_reference(top_accession, paths)
      genotype(n_reads, paths)
      map_reads(reference_fa_path, paths)
      assess_coverage(reference_len, paths)
      trim(paths)
      assemble(normalise(norm_k_list, norm_c_list, paths), asm_k_list, asm_untrusted_contigs, reference_found, paths, threads)
      evaluate_assemblies(reference_found, paths, threads)
      remap_reads()
      report(paths)

argh.dispatch_command(main)