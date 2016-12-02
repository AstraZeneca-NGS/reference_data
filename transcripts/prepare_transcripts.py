#!/usr/bin/env python
import os
from collections import defaultdict
from os.path import join, basename, splitext, isdir, isfile, dirname

from ngs_utils.file_utils import verify_dir, safe_mkdir, verify_file, splitext_plus, adjust_path
from ngs_utils.logger import warn, info

REPLACEMENTS = {
    'FGFR3':    'ENST00000440486',
    'CDKN2A':   'ENST00000304494',
    'ESR1':     'ENST00000206249',
    'MET':      'ENST00000397752',
    'PPP2R2A':  'ENST00000380737',
    'RAD54L':   'ENST00000371975',
    'RAD51D':   'ENST00000345365',
    'AKT1':     'ENST00000555528',
    'CD79B':    'ENST00000006750',
    'CHEK1':    'ENST00000534070',
    'BRCA1':    'ENST00000357654',
    'FANCL':    'ENST00000233741',
    'MYD88':    'ENST00000396334',
    'CHEK2':    'ENST00000328354',
    'FGFR1':    'ENST00000447712',
    'FGFR2':    'ENST00000457416',
}

"""
Chihua:
    export PATH=/home/vsaveliev/bcbio/anaconda/bin:$PATH && /home/vsaveliev/bcbio_tools/bin/snpEff -Xms750m -Xmx4g \
        -dataDir /home/vsaveliev/bcbio/genomes/Hsapiens/hg38/snpeff
        -canon GRCh37.75 -d -v 2>&1 | tee snpeff_cancer_verbose_output_grch37.75.txt.txt

    export PATH=/home/vsaveliev/bcbio/anaconda/bin:$PATH && /home/vsaveliev/bcbio_tools/bin/snpEff -Xms750m -Xmx4g \
        -dataDir /home/vsaveliev/bcbio/genomes/Hsapiens/hg38/snpeff
        -canon GRCh38.82 -d -v 2>&1 | tee snpeff_cancer_verbose_output_grch38.82.txt.txt

Local:
    python repare_transcripts.py
"""

def later_tx_version(t1, t2):
    return t1.split('.') > t2.split('.')

def eq_tx(t1, t2):
    return t1.split('.')[0] == t2.split('.')[0]

base_dir = adjust_path(join(dirname(__file__)))

def get_transcripts_from_snpeff_output(snpeff_output_fname):
    tx_per_genev_by_gene = defaultdict(lambda: defaultdict(list))

    snpeff_output_fpath = join(base_dir, snpeff_output_fname)
    with open(snpeff_output_fpath) as f:
        transcripts_started = False
        for l in f:
            if transcripts_started:
                fs = l.strip().split()
                if 'done.' in fs:
                    break
                if len(fs) != 4 or 'geneName' in fs:
                    continue
                else:
                    g = fs[0]
                    gv = fs[1]
                    t = fs[2]
                    is_met = False
                    for prev_t in tx_per_genev_by_gene[g][gv]:
                        if eq_tx(t, prev_t):
                            is_met = True
                            # info(t + ' has the same base tx id as ' + prev_t)
                            if later_tx_version(t, prev_t):
                                # info(t + ' is later than ' + prev_t + ', replacing')
                                tx_per_genev_by_gene[g][gv][tx_per_genev_by_gene[g][gv].index(prev_t)] = t
                                break
                            else:
                                # info(t + ' is earlier than ' + prev_t + ', skipping')
                                break
                    if not is_met:
                        tx_per_genev_by_gene[g][gv].append(t)

            elif 'Canonical transcripts:' in l:
                transcripts_started = True
    return tx_per_genev_by_gene

info('Getting transcripts for hg19')
hg19_tx_per_genev_by_gene = get_transcripts_from_snpeff_output(join(base_dir, 'snpeff_cancer_verbose_output_grch37.75.txt'))
hg19_genes = set(hg19_tx_per_genev_by_gene.keys())

info('Getting transcripts for hg38')
hg38_tx_per_genev_by_gene = get_transcripts_from_snpeff_output(join(base_dir, 'snpeff_cancer_verbose_output_grch38.82.txt'))
hg38_genes = set(hg38_tx_per_genev_by_gene.keys())

info('Genes present only in hg19: ' + str(len(hg19_genes - hg38_genes)))
info('Genes present only in hg38: ' + str(len(hg38_genes - hg19_genes)))

canon_tx_hg19_fpath = join(base_dir, 'cancer_transcripts_hg19.txt')
canon_tx_hg38_fpath = join(base_dir, 'cancer_transcripts_hg38.txt')

not_matching_tr_count = 0
mult_hg19_tx_count = 0
mult_hg38_tx_count = 0
with open(canon_tx_hg19_fpath, 'w') as hg19, open(canon_tx_hg38_fpath, 'w') as hg38:
    for g in hg38_genes:
        hg19_trs_by_gv = hg19_tx_per_genev_by_gene[g]
        hg38_trs_by_gv = hg38_tx_per_genev_by_gene[g]

        if g in REPLACEMENTS:
            repl_t = REPLACEMENTS[g]

            print g + ' in replacement list as ' + repl_t

            print '  hg19: replacing all gene versions: ' + ', '.join(hg19_trs_by_gv.keys())
            for gv, trs in hg19_trs_by_gv.items():
                print '     replacing gene version ' + gv + ': ' + ', '.join(trs) + ' -> ' + repl_t
                hg19_trs_by_gv[gv] = [repl_t]

            print '  hg38:'
            for gv, trs in hg38_trs_by_gv.items():
                print '     replacing gene version ' + gv + ': ' + ', '.join(trs) + ' -> ' + repl_t
                hg38_trs_by_gv[gv] = [repl_t]
                is_met = False
                for t in trs:
                    if eq_tx(t, repl_t):
                        print '    ' + t + ' eq to ' + repl_t
                        is_met = True
                if not is_met:
                    print '  not met in ' + ', '.join(trs)

        if not set(hg19_trs_by_gv) == set(hg38_trs_by_gv):
            # print 'Transcripts do not match:\n  hg19: ' + ', '.join(hg19_trs_by_gene[g]) + '\n  hg38: ' + ', '.join(hg38_trs_by_gene[g])
            not_matching_tr_count += 1

        for gv, trs in hg19_trs_by_gv.items():
            if len(trs) > 1:
                # print 'Multiple hg19 tx for ' + g + ': ' + ', '.join(hg38_trs_by_gene[g])
                mult_hg19_tx_count += 1
            for t in trs:
                hg19.write(t + '\n')

        for gv, trs in hg38_trs_by_gv.items():
            if len(trs) > 1:
                # print 'Multiple hg38 tx for ' + g + ': ' + ', '.join(hg38_trs_by_gene[g])
                mult_hg38_tx_count += 1
            for t in trs:
                hg38.write(t + '\n')

print '---'
print 'total genes:', len(hg38_genes)
print 'not_matching_tr_count:', not_matching_tr_count
print 'mult_hg38_tx_count:', mult_hg38_tx_count
print 'mult_hg19_tx_count:', mult_hg19_tx_count







