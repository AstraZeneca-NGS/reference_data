import csv
from os.path import join
import sys

from ngs_utils.logger import debug, critical

import bed_annotation.ensembl as ebl


bm_fpath = ebl.biomart_fpath()
if not bm_fpath:
    critical('Error: biomart file is not found')


with open(bm_fpath) as inp, \
     open('ensembl_to_refseq.tsv', 'w') as out:
    out = csv.DictWriter(out, delimiter='\t', fieldnames=['Ensembl', 'RefSeq'])
    for r in csv.DictReader(inp, delimiter='\t'):
        refseq_tx = r['RefSeq ncRNA [e.g. NR_002834]'] or r['RefSeq mRNA [e.g. NM_001195597]']
        if refseq_tx.startswith('N'):
            out.writerow({
                'Ensembl': r['Transcript ID'],
                'RefSeq': refseq_tx,
            })
