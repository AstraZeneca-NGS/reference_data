import csv


with open('biomart.tsv') as inp, \
     open('ensembl_to_refseq.tsv.tx', 'w') as out:
    out = csv.DictWriter(out, delimiter='\t', fieldnames=[
        'Ensembl', 'RefSeq', 'TSL', 'HUGO'])
    out.writeheader()
    for r in csv.DictReader(inp, delimiter='\t'):
        tx = r['RefSeq']
        if tx.startswith('N'):
            out.writerow({
                'Ensembl': r['Ensembl'],
                'RefSeq': tx,
                'TSL': r['Transcript Support Level (TSL)'],
                'HUGO': r['HGNC symbol'],
            })
