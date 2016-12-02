from os.path import isfile, join, basename, splitext
import os
from datetime import datetime

script = "/Users/vladsaveliev/vagrant/TargQC/GeneAnnotation/annotate_bed.py"
test_dir = "/Users/vladsaveliev/vagrant/TargQC/GeneAnnotation/Ensembl/hg38/test_annotation/"

def run(args, fname, suf):
    fpath = join(test_dir, fname)
    output_fpath = join(test_dir, splitext(fname)[0] + '.anno' + (('_'+suf) if suf else '') + '.bed')
    prev_output_fpath = None
    if isfile(output_fpath):
        prev_output_fpath = output_fpath + '_' + datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
        os.rename(output_fpath, prev_output_fpath)

    cmdl = script + ' ' + fpath + ' ' + args + ' -o ' + output_fpath
    print cmdl
    os.system(cmdl)
    print ''

    if prev_output_fpath:
        os.system('diff ' + prev_output_fpath + ' ' + output_fpath)

BED_FNAME = 'small.bed'

run('-g hg38',                                                          BED_FNAME, 'default')
run('-g hg38 --extended',                                               BED_FNAME, 'extended')
run('-g hg38 --extended --output-features',                             BED_FNAME, 'extended_plus_features')
run('-g hg38 --extended --output-features --high-confidence',           BED_FNAME, 'extended_plus_features_high_confidence')
run('-g hg38 --short',                                                  BED_FNAME, 'short')
run('-g hg38 --seq2c',                                                  BED_FNAME, 'seq2c')
