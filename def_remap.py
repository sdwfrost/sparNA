#!/usr/bin/env python3

import os 
import networkx
import matplotlib
from Bio import SeqIO

asm = '1.060-660_r1_Cap1.norm_k25c2.asm_k21.uc0'
paths = {
 'o': '/Users/Bede/Research/Analyses/phe_asm/phe_hcv_pipeline/tmp/run_1445604012_fast_3.6',
 'sample_name':'060-660_r1_Cap1'}
threads = 12

def map_to_subgraphs(paths, threads, i=1):
    for dir in os.listdir(paths['o'] + '/subgraph'):
        if not dir.startswith('.'):
            print('Aligning reads to subgraph of {}... (BWA)'.format(dir))
            cmd_vars = {
             'i':str(i),
             'dir':dir,
             'path_o':paths['o'],
             'ref':paths['o'] + '/subgraph/' + dir + '/contigs.fasta',
             'threads':threads}           
            cmds = [
             'bwa index {path_o}/subgraph/{dir}/contigs.fasta',
             'bwa mem -t {threads} {path_o}/subgraph/{dir}/contigs.fasta '
             '{path_o}/merge/{i}.{sample_name}.raw.r12.fastq '
             '> {path_o}/subgraph/{dir}/contigs.mapped.sam ',
             # count reads mapped and uniquely mapped reads ith flagstats
            cmds = [cmd.format(**cmd_vars) for cmd in cmds]
            for cmd in cmds:
                logger.info(cmd)
                cmd_run = run(cmd)
                logger.info(cmd_run.stdout)
                cmd_prefix = cmd.split(' ')[0]
                print('\tDone (' + cmd_prefix + ')') if cmd_run.returncode == 0 else sys.exit('ERR_PREMAP')
            # Return mapping stats
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

map_to_subgraphs(paths, threads)
