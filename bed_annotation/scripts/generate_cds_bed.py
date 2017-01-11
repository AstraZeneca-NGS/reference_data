#!/usr/bin/env python

import os
import shutil
from optparse import OptionParser, SUPPRESS_HELP
from os.path import isfile, join, basename, dirname, pardir

from ngs_utils.logger import debug
from ngs_utils.file_utils import file_transaction, adjust_path, safe_mkdir, verify_file
from ngs_utils.logger import critical, info
from ngs_utils import logger

import ensembl as ebl


def main():
    options = [
        (['-g', '--genome'], dict(
            dest='genome',
            help='Genome build. Accepted values: ' + ', '.join(ebl.SUPPORTED_GENOMES),
        )),
    ]
    parser = OptionParser()
    for args, kwargs in options:
        parser.add_option(*args, **kwargs)
    opts, args = parser.parse_args()

    if not opts.genome:
        critical('Error: please, specify genome build name with -g (e.g. `-g hg19`)')
    genome = opts.genome

    debug('Getting features from storage')
    features_bed = ebl.get_all_features(genome)
    if features_bed is None:
        critical('Genome ' + genome + ' is not supported. Supported: ' + ', '.join(ebl.SUPPORTED_GENOMES))

    info('Extracting features from Ensembl GTF')
    features_bed = features_bed.filter(lambda x: x[ebl.BedCols.FEATURE] == 'CDS')
    features_bed = features_bed.filter(ebl.get_only_canonical_filter(genome))

    info('Saving CDS regions...')
    output_fpath = adjust_path(join(dirname(__file__), pardir, genome, 'bed', 'CDS-canonical.bed'))
    with file_transaction(None, output_fpath) as tx:
        features_bed.cut(range(6)).saveas(tx)
    info('Done, saved to ' + output_fpath)


if __name__ == '__main__':
    main()
