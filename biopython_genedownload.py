#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 17:06:59 2020

@author: ibarlow
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import Entrez
import pandas as pd
from pathlib import Path
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as plt
import re


Entrez.email = "i.barlow@lms.mrc.ac.uk"  # Always tell NCBI who you are


def parse_gene_list(HGNC_FILE):
    with open(HGNC_FILE, 'r') as fid:
        genes = fid.read()

    genes = set([g for g in genes.split(',')])

    # convert gene symbols to entrez ids
    id_list = []
    for x in genes:
        sterm = 'Homo sapiens[orgn] AND {}[gene]'.format(x)
        handle = Entrez.esearch(db="nucleotide",
                                term=sterm,
                                idtype="acc")
        record = Entrez.read(handle)
        IDArray = record["IdList"]
        # toString = str(IDArray[0])
        id_list.append((x, IDArray))

    # and now find the sequences
    record_df = []
    for i in id_list:
        print('processing {}'.format(i[0]))
        if len(i[1]) > 1:
            print('retrieving multiple id entries')
            for entry in i[1]:
                handle = Entrez.efetch(db="nucleotide",
                                       id=entry,
                                       rettype="gb",
                                       retmode="text")
                record = SeqIO.read(handle, "genbank")
                if i[0] in record.description:
                    record_df.append(pd.Series({
                        'HGNC': i[0],
                        'Sequence': str(record.seq),
                        'Description': record.description,
                        'entrez_id': entry,
                        'sequence_length': len(record.seq)
                        }).to_frame().transpose())
    record_df = pd.concat(record_df).reset_index(drop=True)
    return record_df


def save_record_fasta(record_df, save_to_dir):
    records_grouped = record_df.groupby('HGNC')
    genes = list(records_grouped.groups.keys())

    for g in genes:
        records_grouped.get_group(g)
        out_dir = save_to_dir / g
        if not out_dir.exists():
            out_dir.mkdir()
        out_file = out_dir / '{}_sequences.fasta'.format(g)
        seq_list = []
        for i, r in records_grouped.get_group(g).iterrows():
            seq_list.append(SeqRecord(Seq(r.Sequence,
                                          IUPAC.IUPACAmbiguousDNA()),
                                      id=r.entrez_id,
                                      name=r.HGNC,
                                      description=r.Description))
            SeqIO.write(seq_list, out_file, 'fasta')
    return seq_list


def align_hits(fasta_file, record_df):
    """ Use clustalw to align the fasta files to find the best sequence to use
    for the C. elegans vs Human comparison

    REQUIREMENTS - download clustaw from http://www.clustal.org/download/current/
    and put folder in Applications

    Input : fasta_file

    Output """

    from Bio.Align.Applications import ClustalwCommandline
    from Bio import AlignIO
    from Bio import Phylo

    gene = fasta_file.parent.stem
    print('analysis {}'.format(gene))

    # check if alignment has already been done
    if len(list(fasta_file.parent.rglob('*.aln')))>0:
        print('{} alignment already done, nothing to do here'.format(gene))

        return

    # import information about the gene from dataframe
    records = record_df[record_df.HGNC==gene].copy()

    # now do the alignment
    clustalw_exe = r"/Applications/clustalw-2.1-macosx/clustalw2"
    clustalw_cline = ClustalwCommandline(clustalw_exe,
                                         infile =fasta_file,
                                         stats=fasta_file.parent / 'stats.txt')
    stdout, stderr = clustalw_cline()

    #find alignment files
    align_file = list(fasta_file.parent.rglob('*.aln'))[0]
    tree_file = list(fasta_file.parent.rglob('*.dnd'))[0]

    alignment = AlignIO.read(align_file, "clustal")

    # find consensus sequence
    consensus = re.finditer(r"\*",
                            alignment.column_annotations['clustal_consensus'])
    clist = []
    for c in consensus:
        clist.append(c.span(0))

    if len(clist) > 0:
        consensus = alignment[:, clist[0][0]:clist[-1][1]]
    else:
        consensus=alignment[:, ::]

    gap_count = {}
    for sequence in consensus:
        gap_count[sequence.id] = sequence.seq.count('-')

    records.loc[:, 'alignment_gaps'] = records.entrez_id.map(gap_count)
    records.sort_values(by=['alignment_gaps', 'sequence_length'],
                        ascending=[True, False],
                        inplace=True)
    records.reset_index(drop=True, inplace=True)

    #save top ranked to output file
    top_sequence = SeqRecord(Seq(records.Sequence.loc[0],
                                 IUPAC.IUPACAmbiguousDNA()),
                             id=records.entrez_id.loc[0],
                             name=gene)
    SeqIO.write(top_sequence,
                fasta_file.parent / '{}_sequence.fa'.format(top_sequence.id),
                'fasta')

    tree = Phylo.read(tree_file, "newick")
    tree.ladderize()
    Phylo.draw(tree)
    plt.savefig(tree_file.parent / 'tree.png')
    plt.close('all')

    return


if __name__ == '__main__':
    HUMAN_GENE_LIST = Path('/Volumes/behavgenom$/Ida/CRISPRs/DiseaseModels/' +
                       'NeuroDiseaseModelling/v2/RefiningGenes/Round2/' +
                       'TopGenes32_HGNC.csv')
    HUMAN_GENE_DIR = HUMAN_GENE_LIST.parent.parent / 'HGNCsequences'

    if not HUMAN_GENE_DIR.exists():
        HUMAN_GENE_DIR.mkdir()

    record_df = parse_gene_list(HUMAN_GENE_LIST)
    record_df.to_csv(HUMAN_GENE_DIR / 'sequenceDF.csv', index=False)
    save_record_fasta(record_df, HUMAN_GENE_DIR)

    fasta_files = list(HUMAN_GENE_DIR.rglob('*.fasta'))

    for f in fasta_files:
        align_hits(f, record_df)
    # for file in fasta_files:

