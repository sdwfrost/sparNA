#!/usr/bin/env python3

import networkx
import matplotlib
from Bio import SeqIO

# asm = networkx.Graph()

# with open('/Users/Bede/Research/Analyses/phe_asm/phe_hcv_pipeline/tmp/run_1443791863_phe_sample8/asm/1.Sample8.HIV-generic.1-0..norm_k25c2.asm_k21,33,55,77.rg0/contigs.fastg', 'r') as fastg_file:
#     for node in SeqIO.parse(fastg_file, 'fasta'):
#         node_header = node.id[:-1].replace("'","")
#         if ':' in node_header:
#             node_name, node_neighbors = node_header.split(':')
#             node_neighbors = node_neighbors.split(',')
#             for node_neighbor in node_neighbors:
#                 if (node_name, node_name) not in asm.edges():
#                     asm.add_edge(node_name, node_neighbor)

# starting_node = 'NODE_1_length_6578_cov_450.665_ID_16281'

# sub_asm_nodes = networkx.node_connected_component(asm, starting_node)
# sub_asm = asm.subgraph(sub_asm_nodes)
# sub_asm_nodes = sub_asm.nodes()
# sub_asm_node_lens = [int(node.split('length_')[1].split('_cov')[0]) for node in sub_asm_nodes]
# sub_asm_node_labels = [str(item) for item in sub_asm_node_lens]
# sub_asm_nodes_labels = dict(zip(sub_asm_nodes, sub_asm_node_labels))

# positions = networkx.spring_layout(sub_asm)
# networkx.draw(sub_asm, pos=positions, node_size=sub_asm_node_lens, with_labels=False)
# networkx.draw_networkx_labels(sub_asm, pos=positions, labels=sub_asm_nodes_labels)

# matplotlib.pyplot.show()

paths = {'o': '/Users/Bede/Research/Analyses/phe_asm/phe_hcv_pipeline/tmp/run_1443791863_phe_sample8'}

def fetch_subgraph_contigs(contigs_parent_dir, paths, i=1):
# def fetch_subgraph_contigs(contigs_parent_dir, paths, i=1):
    '''
    Fetch any contigs with connectivity to the longest assembly contig by parsing FASTG output.
    Assumes that the longest contig is the first contig, which it is for SPAdes.
    Canonicalise forward and reverse complement nodes
    Handle cases where longest contig is unconnected
    Returns subgraph of contigs and its number of nodes
    '''
    graph = networkx.Graph()
    longest_contig = None
    longest_contig_unconnected = None
    headers_cnames = {}
    with open(contigs_parent_dir + '/contigs.fastg', 'r') as contigs_fastg_file:
        for record in SeqIO.parse(contigs_fastg_file, 'fasta'): # treat fastg as fasta
            node_name, node_neighbors = None, None
            canonicalised_header = record.id[:-1].replace("'","")
            if ':' in canonicalised_header: # is node connected?
                node_name, node_neighbors = canonicalised_header.split(':')
                if longest_contig is None:
                   longest_contig = node_name
                node_neighbors = node_neighbors.split(',')
                for node_neighbor in node_neighbors:
                    if (node_name, node_neighbor) not in graph.edges():
                        graph.add_edge(node_name, node_neighbor)
            else:
                node_name = canonicalised_header
            if longest_contig is None:
                longest_contig = node_name
                graph.add_node(node_name)

    subgraph = graph.subgraph(networkx.node_connected_component(graph, longest_contig))
    subgraph_nodes = subgraph.nodes()
    subgraph_node_lens = [int(node.split('length_')[1].split('_cov')[0]) for node in subgraph_nodes]
    subgraph_nodes_lens = dict(zip(subgraph_nodes, subgraph_node_lens))

    with open(contigs_parent_dir + '/contigs.fasta', 'r') as contigs_file:
        with open(paths['o'] + '/remap/subgraph_contigs.fasta', 'w') as subgraph_contigs_file:
            for record in SeqIO.parse(contigs_file, 'fasta'):
                if record.id in subgraph_nodes:
                    SeqIO.write(record, subgraph_contigs_file, 'fasta')

    subgraph_node_labels = [str(item) for item in subgraph_node_lens]
    subgraph_nodes_labels = dict(zip(subgraph_nodes, subgraph_node_labels))
    positions = networkx.spring_layout(subgraph)

    # networkx.draw(subgraph, pos=positions, node_size=subgraph_node_lens, with_labels=False)
    # networkx.draw_networkx_labels(subgraph, pos=positions, labels=subgraph_nodes_labels)
    # matplotlib.pyplot.show()

    return subgraph, len(subgraph_nodes)

contigs_parent_dir = '/Users/Bede/Research/Analyses/phe_asm/phe_hcv_pipeline/tmp/run_1443791863_phe_sample8/asm/1.Sample8.HIV-generic.1-0..norm_k25c2.asm_k21,33,55,77.rg0/'
fetch_subgraph_contigs(contigs_parent_dir, paths)
