from os.path import isfile, join, splitext, dirname, pardir, abspath
import os
from datetime import datetime


script = abspath(join(dirname(__file__), pardir, 'annotate_bed.py'))
test_dir = dirname(__file__)


def run(args, fname, suf):
    fpath = join(test_dir, fname)
    output_fpath = join(test_dir, splitext(fname)[0] + '.anno' + (('_'+suf) if suf else '') + '.bed')
    prev_output_fpath = None
    if isfile(output_fpath):
        prev_output_fpath = output_fpath + '_' + datetime.now().strftime('%Y_%m_%d_%H_%M_%S')
        os.rename(output_fpath, prev_output_fpath)

    cmdl = script + ' ' + fpath + ' ' + args + ' -o ' + output_fpath
    print cmdl
    os.system(cmdl)
    print ''

    if prev_output_fpath:
        os.system('diff ' + prev_output_fpath + ' ' + output_fpath + ' > ' + output_fpath + '.diff')


BED_FNAME = 'test.bed'  # small.bed


run('-g hg19',                                      BED_FNAME, 'default')
run('-g hg19 --extended',                           BED_FNAME, 'extended')
run('-g hg19 --extended --output-features',         BED_FNAME, 'extended_plus_features')
run('-g hg19 --extended --canonical',               BED_FNAME, 'extended_canonical')
run('-g hg19 --extended --ambiguities best_one',    BED_FNAME, 'extended_canonical_bestone')
run('-g hg19 --extended --ambiguities best_all',    BED_FNAME, 'extended_canonical_bestall')
run('-g hg19 --extended --ambiguities all',         BED_FNAME, 'extended_canonical_all')
run('-g hg19 --short',                              BED_FNAME, 'short')


# TODO:
# 1. check missing other cols info:
# /Users/vlad/vagrant/reference_data/bed_annotation/test/small.anno_extended_plus_features_canonical.bed
# 2. compare results in regular and canonical runs
# 3. check ambiguous annotations in std annotation run
# 4. annotate bed files
# 5. get hg38 bed files if possible from vendor; otherwise map; annotate
