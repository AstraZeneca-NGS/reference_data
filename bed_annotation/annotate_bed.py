#!/usr/bin/env python

import os
import shutil
from optparse import OptionParser, SUPPRESS_HELP
from collections import defaultdict, OrderedDict
from os.path import isfile, join, basename
from tempfile import mkdtemp

from ngs_utils.logger import debug
from ngs_utils.bed_utils import verify_bed, SortableByChrom, count_bed_cols, sort_bed, clean_bed
from ngs_utils.file_utils import file_transaction, adjust_path, safe_mkdir, verify_file
from ngs_utils.logger import critical, info
from ngs_utils import logger

import ensembl.utils as eu
from bed_annotation.bed_annotation import annotate


def main():
    options = [
        (['-o', '--output-file'], dict(
            dest='output_file',
            metavar='FILE',
            help='Output file',
        )),
        (['--output-features'], dict(
            dest='output_features',
            action='store_true',
            default=False,
            help='Also output featues that was used to annotate',
        )),
        (['--reuse'], dict(
            dest='reuse_intermediate',
            action='store_true',
            help='reuse intermediate non-empty files in the work dir from previous run',
        )),
        (['-g', '--genome'], dict(
            dest='genome',
            help='Genome build. Accepted values: ' + ', '.join(eu.SUPPORTED_GENOMES),
        )),
        (['--canonical'], dict(
            dest='only_canonical',
            action='store_true',
            default=False,
            help='Use only features from canonical transcripts to annotate',
        )),
        (['--short'], dict(
            dest='short',
            action='store_true',
            default=False,
            help='Add only "Gene" column (outputa 4-col BED file instead of 6-col)',
        )),
        (['--cds-only'], dict(
            dest='cds_only',
            action='store_true',
            default=False,
            help='Use only CDS to annotate',
        )),
        (['--extended'], dict(
            dest='extended',
            action='store_true',
            default=False,
            help='Output additional columns: transcript, GC, overlap size...',
        )),
        (['--high-confidence'], dict(
            dest='high_confidence',
            action='store_true',
            default=False,
            help='Annotate with only high confidence regions (TSL is 1 or NA, with HUGO symbol, total overlap size > 50%)',
        )),
        (['--seq2c'], dict(
            dest='seq2c',
            action='store_true',
            default=False,
            help='Uses --canonical, and uses only one annotation per gene',  # TODO: prefer consecutive annotations
        )),
        (['--debug'], dict(
            dest='debug',
            action='store_true',
            default=False,
            help=SUPPRESS_HELP,
         )),
        (['--not-collapse-exons'], dict(
            dest='collapse_exons',
            action='store_false',
            default=True,
            help=SUPPRESS_HELP,
         )),
        (['--work-dir'], dict(dest='work_dir', metavar='DIR', help=SUPPRESS_HELP)),
        (['--log-dir'], dict(dest='log_dir', metavar='DIR', help=SUPPRESS_HELP)),
    ]

    parser = OptionParser(description='Annotating BED file based on reference features annotations.')
    for args, kwargs in options:
        parser.add_option(*args, **kwargs)
    opts, args = parser.parse_args()
    logger.is_debug = opts.debug

    if not opts.genome:
        critical('Error: please, specify genome build name with -g (e.g. `-g hg19`)')

    if opts.short:
        if opts.extended:        critical('--short and --extended can\'t be set both')
        if opts.output_features: critical('--short and --output-features can\'t be set both')
    elif opts.output_features or opts.extended:
        opts.extended = True
        opts.short = False
    if opts.seq2c:
        opts.only_canonical = True  #opts.cds_only = True

    if len(args) < 1:
        parser.exit('Usage: ' + __file__ + ' Input_BED_file -g hg19 -o Annotated_BED_file [--canonical]')
    input_bed_fpath = verify_file(args[0], is_critical=True, description='Input BED file for ' + __file__)
    output_fpath = adjust_path(opts.output_file)

    # prev_output_fpath = None
    # if opts.debug:
    #     if isfile(output_fpath):
    #         prev_output_fpath = output_fpath + '_' + datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
    #         os.rename(output_fpath, prev_output_fpath)
    if opts.work_dir:
        work_dir = join(adjust_path(opts.work_dir), os.path.splitext(basename(input_bed_fpath))[0])
        safe_mkdir(work_dir)
        info('Created work directory ' + work_dir)
    else:
        work_dir = mkdtemp('bed_annotate')
        debug('Created temporary work directory ' + work_dir)

    input_bed_fpath = clean_bed(input_bed_fpath, work_dir)
    input_bed_fpath = verify_bed(input_bed_fpath, is_critical=True, description='Input BED file for ' + __file__)

    output_fpath = annotate(
        input_bed_fpath, output_fpath, work_dir, genome=opts.genome, is_debug=opts.debug,
        only_canonical=opts.only_canonical, short=opts.short, extended=opts.extended,
        high_confidence=opts.high_confidence, collapse_exons=opts.collapse_exons,
        output_features=opts.output_features, cds_only=opts.cds_only, for_seq2c=opts.seq2c)

    if not opts.work_dir:
        debug('Removing work directory ' + work_dir)
        shutil.rmtree(work_dir)

    # if opts.debug:
    #     if prev_output_fpath:
    #         os.system('diff ' + prev_output_fpath + ' ' + output_fpath)
    info('Done, saved to ' + output_fpath)


if __name__ == '__main__':
    main()
