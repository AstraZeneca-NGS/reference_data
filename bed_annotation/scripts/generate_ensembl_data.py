#!/usr/bin/env python
import os
from os.path import join, isfile
from optparse import OptionParser
import csv

import ngs_utils.reference_data as ref
from ngs_utils import gtf
from ngs_utils import logger
from ngs_utils.bed_utils import bgzip_and_tabix, sort_bed
from ngs_utils.file_utils import verify_file, add_suffix
from ngs_utils.logger import debug, warn, critical

import ensembl as ebl


def main():
    description = '''
Usage:
    ' + __file__ + ' hg19 [db.gtf]
'''

    options = [
        (['--debug'], dict(dest='debug', action='store_true', default=False)),
    ]
    parser = OptionParser(description=description)
    for args, kwargs in options:
        parser.add_option(*args, **kwargs)
    opts, args = parser.parse_args()
    if len(args) == 0:
        parser.exit(1, 'Please provide genome name as the first argument')
    logger.is_debug = opts.debug

    genome_name = args[0]

    if len(args) > 1:
        gtf_fpath = args[1]
    else:
        gtf_fpath = ebl.ensembl_gtf_fpath(genome_name)
    if not isfile(gtf_fpath):
        if not gtf_fpath.endswith('.gz'):
            gtf_fpath += '.gz'
    gtf_fpath = verify_file(gtf_fpath)
    debug('Reading the GTF database')
    db = gtf.get_gtf_db(gtf_fpath)

    debug('Reading biomart data')
    features_by_ens_id = read_biomart(genome_name)

    chroms = [c for c, l in ref.get_chrom_lengths(genome_name)]
    
    output_fpath = join('Ensembl', genome_name, 'ensembl.bed')
    unsorted_output_fpath = add_suffix(output_fpath, 'unsorted')
    debug('Processing features, writing to ' + unsorted_output_fpath)

    def _get(_rec, _key):
        val = _rec.attributes.get(_key)
        if val is None:
            return None
        assert len(val) == 1, (_key, str(val))
        return val[0]

    num_tx_not_in_biomart = 0
    num_tx_diff_gene_in_biomart = 0
    with open(unsorted_output_fpath, 'w') as out:
        out.write('\t'.join(ebl.BedCols.names[i] for i in ebl.BedCols.cols[:-4]) + '\n')

        for rec in db.all_features(order_by=('seqid', 'start', 'end')):
            if rec.featuretype == 'gene':
                continue
            
            if rec.chrom not in chroms:
                continue
                
            tx_id = _get(rec, 'transcript_id')
            gname = _get(rec, 'gene_name')
            tx_biotype = _get(rec, 'transcript_biotype')
            tsl = _get(rec, 'transcript_support_level')

            biomart_rec = features_by_ens_id.get(tx_id)
            if not biomart_rec:
                if rec.featuretype == 'transcript':
                    num_tx_not_in_biomart += 1
                continue

            bm_gname = biomart_rec['Associated Gene Name']
            bm_tx_biotype = biomart_rec['Transcript type']
            bm_tsl = biomart_rec.get('Transcript Support Level (TSL)')
            hugo_gene = biomart_rec['HGNC symbol']

            if bm_gname != gname:
                if rec.featuretype == 'transcript':
                    num_tx_diff_gene_in_biomart += 1
                continue
            tx_biotype = bm_tx_biotype
            if rec.end - rec.start < 0:
                continue
            tsl = bm_tsl.split()[0].replace('tsl', '') if bm_tsl else None

            fs = [None] * len(ebl.BedCols.cols[:-3])
            if not rec.chrom.startswith('chr'):
                rec.chrom = 'chr' + rec.chrom.replace('MT', 'M')
            fs[:6] = [rec.chrom,
                      str(rec.start - 1),
                      str(rec.end),
                      gname,
                      rec.attributes.get('exon_number', ['.'])[0],
                      rec.strand]
            fs[ebl.BedCols.FEATURE] = rec.featuretype or '.'
            fs[ebl.BedCols.BIOTYPE] = tx_biotype or '.'
            fs[ebl.BedCols.ENSEMBL_ID] = tx_id or '.'
            # fs[ebl.BedCols.REFSEQ_ID] = refseq_id or '.'
            # fs[ebl.BedCols.IS_CANONICAL] = 'canonical' if refseq_id in canonical_transcripts_ids else ''
            fs[ebl.BedCols.TSL] = tsl or '.'
            fs[ebl.BedCols.HUGO] = hugo_gene or '.'
            # fs[ebl.BedCols.names[ensembl.BedCols.GC]] = gc
            out.write('\t'.join(fs) + '\n')

    if num_tx_not_in_biomart:
        warn(str(num_tx_not_in_biomart) + ' transcripts not found in biomart')
    if num_tx_diff_gene_in_biomart:
        warn(str(num_tx_diff_gene_in_biomart) + ' transcripts have a different gene name in biomart')

    debug('Sorting results')
    sort_bed(unsorted_output_fpath, output_fpath, fai_fpath=ref.get_fai(genome_name), genome=genome_name)
    os.remove(unsorted_output_fpath)
    bgzip_and_tabix(output_fpath)

    # with open(output_fpath, 'w') as out:
    #     for feature in db.features_of_type('transcript', order_by=("seqid", "start", "end")):
    #         chrom = feature.chrom
    #         start = feature.start
    #         end = feature.end
    #         attributes = feature.attributes.keys()
    #         strand = feature.strand
    #         name = (feature['gene_name'][0] if 'gene_name' in attributes else
    #                 feature['gene_id'][0])
    #         line = "\t".join([str(x) for x in [chrom, start, end, name, ".",
    #                                            strand]])
    #         out.write(line + "\n")


    # db_bed = gtf.gtf_to_bed(gtf_db_fpath, out_dir)


def read_biomart(genome_name):
    features_by_ens_id = dict()
    bm_fpath = ebl.biomart_fpath(genome_name)
    if not verify_file(bm_fpath): critical('Biomart file not found, and needed for TSL values')
    with open(bm_fpath) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            features_by_ens_id[r['Transcript ID']] = r
    
    # hg38 version has TSL, checking if we can populate some TSL from it
    if not genome_name.startswith('hg38'):
        bm_fpath = ebl.biomart_fpath('hg38')
        if not verify_file(bm_fpath): critical('Biomart for hg38 file not found, and needed for TSL values')
        with open(bm_fpath) as f:
            for r in csv.DictReader(f, delimiter='\t'):
                if r['Transcript ID'] not in features_by_ens_id:
                    features_by_ens_id[r['Transcript ID']] = r
                else:
                    features_by_ens_id[r['Transcript ID']]['Transcript Support Level (TSL)'] = r[
                        'Transcript Support Level (TSL)']
    return features_by_ens_id


if __name__ == '__main__':
    main()
