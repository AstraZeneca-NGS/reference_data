#!/usr/bin/env python

import os
from collections import defaultdict, OrderedDict
from os.path import isfile, join, basename

from ngs_utils.bedtools import BedTool
from ngs_utils import reference_data
from ngs_utils.logger import debug
from ngs_utils.utils import OrderedDefaultDict
from ngs_utils.bed_utils import verify_bed, SortableByChrom, count_bed_cols, sort_bed, clean_bed
from ngs_utils.file_utils import file_transaction, adjust_path, safe_mkdir, verify_file
from ngs_utils.logger import critical, info

import ensembl.utils as eu


def bed_chrom_order(bed_fpath):
    chroms = []
    chroms_set = set()
    with open(bed_fpath) as f:
        for l in f:
            l = l.strip()
            if l and not l.startswith('#'):
                chrom = l.split('\t')[0]
                if chrom not in chroms_set:
                    chroms_set.add(chrom)
                    chroms.append(chrom)
    chr_order = {c: i for i, c in enumerate(chroms)}
    return chr_order


def get_sort_key(chr_order):
    return lambda fs: (
        chr_order[fs[eu.BedCols.CHROM]],
        int(fs[eu.BedCols.START]),
        int(fs[eu.BedCols.END]),
        0 if fs[eu.BedCols.FEATURE] == 'transcript' else 1,
        fs[eu.BedCols.ENSEMBL_ID],
        fs[eu.BedCols.GENE]
    )


def overlap_with_features(input_bed_fpath, output_fpath, work_dir, genome=None, is_debug=False,
             only_canonical=False, extended=False, high_confidence=False, collapse_exons=True,
             reannotate=False, cds_only=False, for_seq2c=False):
    return annotate(input_bed_fpath, output_fpath, work_dir, genome=genome, is_debug=is_debug,
             only_canonical=only_canonical, extended=extended, high_confidence=high_confidence,
             collapse_exons=collapse_exons, short=False, output_features=True, reannotate=reannotate,
             cds_only=cds_only, for_seq2c=for_seq2c)


canon_tx_ids, cancer_tx_ids = set(), set()


def annotate(input_bed_fpath, output_fpath, work_dir, genome=None, is_debug=False,
             only_canonical=False, extended=False, high_confidence=False, collapse_exons=True,
             short=False, output_features=False, reannotate=True, cds_only=False, for_seq2c=False):
    debug('Getting features from storage')
    features_bed = eu.get_all_features(genome)
    if features_bed is None:
        critical('Genome ' + genome + ' is not supported. Supported: ' + ', '.join(eu.SUPPORTED_GENOMES))

    if genome:
        fai_fpath = reference_data.get_fai(genome)
        chr_order = reference_data.get_chrom_order(genome)
    else:
        fai_fpath = None
        chr_order = bed_chrom_order(input_bed_fpath)

    input_bed_fpath = sort_bed(input_bed_fpath, work_dir=work_dir, chr_order=chr_order, genome=genome)

    bed = BedTool(input_bed_fpath)
    col_num = bed.field_count()
    keep_gene_column = False
    if col_num > 3:
        if reannotate:
            bed = BedTool(input_bed_fpath).cut([0, 1, 2])
            keep_gene_column = False
        else:
            if col_num > 4:
                bed = BedTool(input_bed_fpath).cut([0, 1, 2, 3])
            keep_gene_column = True

    global canon_tx_ids, cancer_tx_ids
    canon_tx_ids, cancer_tx_ids = eu.get_canonical_transcripts_ids(genome)

    # features_bed = features_bed.saveas()
    # cols = features_bed.field_count()
    # if cols < 12:
    #     features_bed = features_bed.each(lambda f: f + ['.']*(12-cols))
    if high_confidence:
        features_bed = features_bed.filter(eu.high_confidence_filter)
    if only_canonical:
        features_bed = features_bed.filter(eu.get_only_canonical_filter(canon_tx_ids, cancer_tx_ids))
    if cds_only:
        features_bed = features_bed.filter(eu.protein_coding_filter)
    # unique_tx_by_gene = find_best_tx_by_gene(features_bed)

    info('Extracting features from Ensembl GTF')
    features_bed = features_bed.filter(lambda x:
        x[eu.BedCols.FEATURE] in ['exon', 'CDS', 'stop_codon', 'transcript'])
        # x[eu.BedCols.ENSEMBL_ID] == unique_tx_by_gene[x[eu.BedCols.GENE]])

    info('Overlapping regions with Ensembl data')
    if is_debug:
        bed = bed.saveas(join(work_dir, 'bed.bed'))
        features_bed = features_bed.saveas(join(work_dir, 'features.bed'))
    annotated = _annotate(bed, features_bed, chr_order, fai_fpath, work_dir,
        high_confidence=high_confidence, collapse_exons=collapse_exons, for_seq2c=for_seq2c,
        output_features=output_features, is_debug=is_debug, keep_gene_column=keep_gene_column)

    header = [eu.BedCols.names[i] for i in eu.BedCols.cols]

    info('Saving annotated regions...')
    with file_transaction(work_dir, output_fpath) as tx:
        with open(tx, 'w') as out:
            if short:
                header = header[:4]
            if not extended:
                header = header[:6]
            if extended:
                out.write('## ' + eu.BedCols.names[eu.BedCols.TX_OVERLAP_PERCENTAGE] +
                          ': part of region overlapping with transcripts\n')
                out.write('## ' + eu.BedCols.names[eu.BedCols.EXON_OVERLAPS_PERCENTAGE] +
                          ': part of region overlapping with exons\n')
                out.write('## ' + eu.BedCols.names[eu.BedCols.CDS_OVERLAPS_PERCENTAGE] +
                          ': part of region overlapping with protein coding regions\n')
            out.write('\t'.join(header) + '\n')
            for fields in annotated:
                if short:
                    fields = fields[:4]
                if not extended:
                    fields = fields[:6]
                out.write('\t'.join(map(_format_field, fields)) + '\n')

    return output_fpath


class Region(SortableByChrom):
    def __init__(self, chrom, start, end, ref_chrom_order, gene_symbol=None, exon=None,
                 strand=None, other_fields=None):
        SortableByChrom.__init__(self, chrom, ref_chrom_order)
        self.chrom = chrom
        self.start = start
        self.end = end
        self.symbol = gene_symbol
        self.exon = exon
        self.strand = strand
        self.other_fields = other_fields or []
        self.total_merged = 0

    def __str__(self):
        fs = [
            self.chrom,
            self.start,
            self.end,
            self.symbol,
            self.exon,
            self.strand
        ] + self.other_fields
        fs = map(_format_field, fs)
        return '\t'.join(fs) + '\n'

    def get_key(self):
        return SortableByChrom.get_key(self), self.start, self.end, self.symbol

    def header(self):
        pass


def _format_field(value):
    if isinstance(value, list) or isinstance(value, set):
        return ', '.join(_format_field(v) for v in value) or '.'
    elif isinstance(value, float):
        return '{:.1f}'.format(value)
    else:
        return str(value or '.')


# def find_best_tx_by_gene(features_bed):
#     all_transcripts = features_bed.filter(lambda x: x[eu.BedCols.FEATURE] in ['transcript'])
#
#     tx_by_gene = defaultdict(list)
#     for tx_region in all_transcripts:
#         tx_by_gene[tx_region.name].append(tx_region)
#     unique_tx_by_gene = dict()
#     for g, txs in tx_by_gene.iteritems():
#         tsl1_txs = [tx for tx in txs if tx[eu.BedCols.TSL] in ['1', 'NA']]
#         if tsl1_txs:
#             txs = tsl1_txs[:]
#         hugo_txs = [tx for tx in txs if tx[eu.BedCols.HUGO]]
#         if hugo_txs:
#             txs = hugo_txs[:]
#         best_tx = max(txs, key=lambda tx_: int(tx_[eu.BedCols.END]) - int(tx_[eu.BedCols.START]))
#         unique_tx_by_gene[g] = best_tx[eu.BedCols.ENSEMBL_ID]
#     return unique_tx_by_gene


def tx_priority_sort_key(x):
    cancer_tx_key = 0 if x[eu.BedCols.ENSEMBL_ID] in cancer_tx_ids else 1

    canon_tx_key = 0 if x[eu.BedCols.ENSEMBL_ID] in canon_tx_ids else 1

    biotype_key = 1
    biotype_rank = ['protein_coding', '__default__', 'rna', 'decay', 'sense_', 'antisense', 'translated_', 'transcribed_']
    for key in biotype_rank:
        if key in x[eu.BedCols.BIOTYPE].lower():
            biotype_key = biotype_rank.index(key)
    tsl_key = {'1': 0, '2': 2, '3': 3, '4': 4, '5': 5}.get(x[eu.BedCols.TSL], 1)
    hugo_key = 0 if x[eu.BedCols.HUGO] not in ['.', '', None] else 1

    overlap_key = 0
    if len(x) > eu.BedCols.TX_OVERLAP_PERCENTAGE:
        overlap_key = [(-x[key] if x[key] is not None else 0)
                       for key in [eu.BedCols.TX_OVERLAP_PERCENTAGE,
                                   eu.BedCols.CDS_OVERLAPS_PERCENTAGE,
                                   eu.BedCols.EXON_OVERLAPS_PERCENTAGE]]

    length_key = -(int(x[eu.BedCols.END]) - int(x[eu.BedCols.START]))

    return cancer_tx_key, canon_tx_key, biotype_key, tsl_key, hugo_key, overlap_key, length_key


# def select_best_tx(overlaps_by_tx):
#     tx_overlaps = [(o, size) for os in overlaps_by_tx.values() for (o, size) in os if o[eu.BedCols.FEATURE] == 'transcript']
#     tsl1_overlaps = [(o, size) for (o, size) in tx_overlaps if o[eu.BedCols.TSL] in ['1', 'NA']]
#     if tsl1_overlaps:
#         tx_overlaps = tsl1_overlaps
#     hugo_overlaps = [(o, size) for (o, size) in tx_overlaps if o[eu.BedCols.HUGO]]
#     if hugo_overlaps:
#         tx_overlaps = hugo_overlaps
#     (best_overlap, size) = max(sorted(tx_overlaps, key=lambda (o, size): int(o[eu.BedCols.END]) - int(o[eu.BedCols.START])))
#     tx_id = best_overlap[eu.BedCols.ENSEMBL_ID]
#     return tx_id


def _resolve_ambiguities(overlaps_by_tx_by_gene_by_loc, chrom_order,
         collapse_exons=True, output_features=False, for_seq2c=False):
    # sort transcripts by "quality" and then select which tx id to report in case of several overlaps
    # (will be useful further when interating exons)
    # annotated_by_tx_by_gene_by_loc = OrderedDict([
    #     (loc, OrderedDict([
    #         (gname, sorted(annotated_by_tx.itervalues(), key=tx_sort_key))
    #             for gname, annotated_by_tx in annotated_by_tx_by_gene.iteritems()
    #     ])) for loc, annotated_by_tx_by_gene in annotated_by_tx_by_gene_by_loc.iteritems()
    # ])

    annotated = []
    for (chrom, start, end), overlaps_by_tx_by_gene in overlaps_by_tx_by_gene_by_loc.iteritems():
        features = dict()
        annotation_alternatives = []
        for gname, overlaps_by_tx in overlaps_by_tx_by_gene.iteritems():
            consensus = [None for _ in eu.BedCols.cols]
            consensus[:3] = chrom, start, end
            if gname:
                consensus[3] = gname
            consensus[eu.BedCols.FEATURE] = 'capture'
            # consensus[eu.BedCols.EXON_OVERLAPS_BASES] = 0
            consensus[eu.BedCols.EXON_OVERLAPS_PERCENTAGE] = 0
            consensus[eu.BedCols.CDS_OVERLAPS_PERCENTAGE] = 0
            consensus[eu.BedCols.EXON] = set()

            if not overlaps_by_tx:  # not annotated or already annotated but gene did not match the intersection
                annotated.append(consensus)
                continue

            all_tx = []
            for xx in overlaps_by_tx.itervalues():
                for x, overlap_bp in xx:
                    if x[eu.BedCols.FEATURE] == 'transcript':
                        x[eu.BedCols.TX_OVERLAP_PERCENTAGE] = 100.0 * overlap_bp / (int(end) - int(start))
                        all_tx.append(x)
            # all_tx = (x for xx in overlaps_by_tx.itervalues() for x, overlap_bp in xx
            #           if x[eu.BedCols.FEATURE] == 'transcript')

            tx_sorted_list = [x[eu.BedCols.ENSEMBL_ID] for x in sorted(all_tx, key=tx_priority_sort_key)]
            if not tx_sorted_list:
                annotated.append(consensus)
                continue
            tx_id = tx_sorted_list[0]
            overlaps = overlaps_by_tx[tx_id]

            # # sort transcripts by "quality" and then select which tx id to report in case of several overlaps
            # # (will be useful further when interating exons)
            # annotated_by_tx_by_gene_by_loc = OrderedDict([
            #     (loc, OrderedDict([
            #         (gname, sorted(annotated_by_tx.itervalues(), key=tx_sort_key))
            #             for gname, annotated_by_tx in annotated_by_tx_by_gene.iteritems()
            #     ])) for loc, annotated_by_tx_by_gene in annotated_by_tx_by_gene_by_loc.iteritems()
            # ])

            for fields, overlap_size in overlaps:
                c_overlap_bp = overlap_size
                c_overlap_pct = 100.0 * c_overlap_bp / (int(end) - int(start))

                if output_features:
                    f_start = int(fields[1])
                    f_end = int(fields[2])
                    feature = features.get((f_start, f_end))
                    if feature is None:
                        feature = [None for _ in eu.BedCols.cols]
                        feature[:len(fields)] = fields
                        # feature[eu.BedCols.TX_OVERLAP_BASES] = 0
                        # feature[eu.BedCols.TX_OVERLAP_PERCENTAGE] = 0
                        features[(f_start, f_end)] = feature
                    # feature[eu.BedCols.TX_OVERLAP_BASES] += c_overlap_bp
                    # feature[eu.BedCols.TX_OVERLAP_PERCENTAGE] += 100.0 * c_overlap_bp / (int(f_end) - int(f_start))
                    # TODO: don't forget to merge BED if not

                if fields[eu.BedCols.FEATURE] == 'transcript':
                    consensus[eu.BedCols.GENE] = gname
                    consensus[eu.BedCols.STRAND] = fields[eu.BedCols.STRAND]
                    consensus[eu.BedCols.BIOTYPE] = fields[eu.BedCols.BIOTYPE]
                    consensus[eu.BedCols.ENSEMBL_ID] = fields[eu.BedCols.ENSEMBL_ID]
                    consensus[eu.BedCols.TSL] = fields[eu.BedCols.TSL]
                    consensus[eu.BedCols.HUGO] = fields[eu.BedCols.HUGO]
                    # consensus[eu.BedCols.TX_OVERLAP_BASES] = c_overlap_bp
                    consensus[eu.BedCols.TX_OVERLAP_PERCENTAGE] = c_overlap_pct

                elif fields[eu.BedCols.FEATURE] == 'exon':
                    consensus[eu.BedCols.EXON_OVERLAPS_PERCENTAGE] += c_overlap_pct
                    consensus[eu.BedCols.EXON].add(int(fields[eu.BedCols.EXON]))

                elif fields[eu.BedCols.FEATURE] == 'CDS':
                    consensus[eu.BedCols.CDS_OVERLAPS_PERCENTAGE] += c_overlap_pct

            consensus[eu.BedCols.EXON] = sorted(list(consensus[eu.BedCols.EXON]))

            # annotated.append(consensus)
            annotation_alternatives.append(consensus)

        # TODO: sort annotation_alternatives by score, select best
        if annotation_alternatives and for_seq2c:
            annotation_alternatives.sort(key=tx_priority_sort_key)
            # annotation_alternatives = [a for a in annotation_alternatives if a[eu.BedCols.CDS_OVERLAPS_PERCENTAGE] > 50]
            best_alt = annotation_alternatives[0]
            annotation_alternatives = [a for a in annotation_alternatives if tx_priority_sort_key(a) == tx_priority_sort_key(best_alt)]

        annotated.extend(annotation_alternatives)

        features = sorted(features.values(), key=get_sort_key(chrom_order))
        if output_features:
            annotated.extend(features)
    return annotated


def _annotate(bed, ref_bed, chr_order, fai_fpath, work_dir,
              high_confidence=False, collapse_exons=True, output_features=False,
              is_debug=False, keep_gene_column=False, for_seq2c=False):
    # if genome:
        # genome_fpath = cut(fai_fpath, 2, output_fpath=intermediate_fname(work_dir, fai_fpath, 'cut2'))
        # intersection = bed.intersect(ref_bed, sorted=True, wao=True, g='<(cut -f1,2 ' + fai_fpath + ')')
        # intersection = bed.intersect(ref_bed, sorted=True, wao=True, genome=genome.split('-')[0])
    # else:

    intersection_bed = None
    intersection_fpath = None
    if is_debug:
        intersection_fpath = join(work_dir, 'intersection.bed')
        if isfile(intersection_fpath):
            info('Loading from ' + intersection_fpath)
            intersection_bed = BedTool(intersection_fpath)
    if not intersection_bed:
        if count_bed_cols(fai_fpath) == 2:
            debug('Fai fields size is 2 ' + fai_fpath)
            intersection_bed = bed.intersect(ref_bed, wao=True, sorted=True, g=fai_fpath)
        else:
            debug('Fai fields is ' + str(count_bed_cols(fai_fpath)) + ', not 2')
            intersection_bed = bed.intersect(ref_bed, wao=True)
    if is_debug and not isfile(intersection_fpath):
        intersection_bed.saveas(intersection_fpath)
        debug('Saved intersection to ' + intersection_fpath)

    total_annotated = 0
    total_uniq_annotated = 0
    total_off_target = 0

    met = set()

    overlaps_by_tx_by_gene_by_loc = OrderedDefaultDict(lambda: OrderedDefaultDict(lambda: defaultdict(list)))
    # off_targets = list()

    for intersection_fields in intersection_bed:
        inters_list = list(intersection_fields)
        if len(inters_list) < 3 + len(eu.BedCols.cols) - 3 + 1:
            critical('Cannot parse the reference BED file - unexpected number of lines '
                     '(' + str(len(inters_list)) + ') in ' + str(inters_list) +
                     ' (less than ' + str(3 + len(eu.BedCols.cols) - 2 + 1) + ')')

        a_chr, a_start, a_end = intersection_fields[:3]

        overlap_fields = [None for _ in eu.BedCols.cols]

        if keep_gene_column:
            a_gene = intersection_fields[3]
            overlap_fields[:len(intersection_fields[4:])] = intersection_fields[4:]
        else:
            a_gene = None
            overlap_fields[:len(intersection_fields[3:])] = intersection_fields[3:]

        e_chr = overlap_fields[0]
        overlap_size = int(intersection_fields[-1])
        assert e_chr == '.' or a_chr == e_chr, str((a_chr + ', ' + e_chr))

        # fs = [None for _ in eu.BedCols.cols]
        # fs[:3] = [a_chr, a_start, a_end]
        reg = (a_chr, int(a_start), int(a_end))

        if e_chr == '.':
            total_off_target += 1
            # off_targets.append(fs)
            overlaps_by_tx_by_gene_by_loc[reg][a_gene] = OrderedDefaultDict(list)

        else:
            # fs[3:-1] = db_feature_fields[3:-1]
            total_annotated += 1
            if (a_chr, a_start, a_end) not in met:
                total_uniq_annotated += 1
                met.add((a_chr, a_start, a_end))

            e_gene = overlap_fields[eu.BedCols.GENE] if not high_confidence else overlap_fields[eu.BedCols.HUGO]
            if keep_gene_column and e_gene != a_gene:
                overlaps_by_tx_by_gene_by_loc[reg][a_gene] = OrderedDefaultDict(list)
            else:
                tx = overlap_fields[eu.BedCols.ENSEMBL_ID]
                overlaps_by_tx_by_gene_by_loc[reg][e_gene][tx].append((overlap_fields, overlap_size))

    info('  Total annotated regions: ' + str(total_annotated))
    info('  Total unique annotated regions: ' + str(total_uniq_annotated))
    info('  Total off target regions: ' + str(total_off_target))
    info('Resolving ambiguities...')
    annotated = _resolve_ambiguities(overlaps_by_tx_by_gene_by_loc, chr_order,
        collapse_exons, output_features, for_seq2c)

    return annotated


def _save_regions(regions, fpath):
    with open(fpath, 'w') as off_target_f:
        for r in regions:
            off_target_f.write(str(r))

    return fpath


def _split_reference_by_priority(cnf, features_bed_fpath):
    features = ['CDS', 'Exon', 'Transcript', 'Gene']
    info('Splitting the reference file into ' + ', '.join(features))
    features_and_beds = []
    for f in features:
        features_and_beds.append((f, BedTool(features_bed_fpath).filter(lambda x: x[6] == f)))
    return features_and_beds
