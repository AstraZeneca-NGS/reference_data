import sys
import os
from os.path import dirname, join, abspath, isfile, pardir

from ngs_utils.bedtools import BedTool
from ngs_utils.bed_utils import bedtools_version
from ngs_utils.file_utils import which, open_gzipsafe, verify_file
from ngs_utils.logger import debug, critical


SUPPORTED_GENOMES = ['hg19', 'hg19-noalt', 'hg38', 'hg38-noalt', 'mm10']

class BedCols:
    CHROM, \
    START, \
    END, \
    GENE, \
    EXON, \
    STRAND, \
    FEATURE, \
    BIOTYPE, \
    ENSEMBL_ID, \
    TSL, \
    HUGO, \
    TX_OVERLAP_PERCENTAGE, \
    EXON_OVERLAPS_PERCENTAGE, \
    CDS_OVERLAPS_PERCENTAGE \
        = cols = range(14)

    names = {
        CHROM: '#Chrom',
        START: 'Start',
        END: 'End',
        GENE: 'Gene',
        EXON: 'Exon',
        STRAND: 'Strand',
        FEATURE: 'Feature',
        BIOTYPE: 'Biotype',
        ENSEMBL_ID: 'Ensembl_ID',
        TSL: 'TSL',
        HUGO: 'HUGO',
        # TX_OVERLAP_BASES: 'Tx_overlap_bp',
        TX_OVERLAP_PERCENTAGE: 'Tx_overlap_%',
        # EXON_OVERLAPS_BASES: 'Exon_overlaps_bp',
        EXON_OVERLAPS_PERCENTAGE: 'Exon_overlaps_%',
        CDS_OVERLAPS_PERCENTAGE: 'CDS_overlaps_%',
    }

def check_genome(genome):
    if genome not in SUPPORTED_GENOMES:
        sys.stdout.write('Genome ' + genome + ' is not supported. Supported genomes: ' + ', '.join(SUPPORTED_GENOMES) + '\n')
        sys.exit(1)

#################
### INTERFACE ###
#################
def get_all_features(genome, high_confidence=False):
    bed = _get_ensembl_file('ensembl.bed', genome)
    if high_confidence:
        bed = bed.filter(lambda x: x[BedCols.HUGO] not in ['', '.', None])
    return bed

def get_merged_cds(genome):
    """
    Returns all CDS merged, used:
    - for TargQC general reports CDS coverage statistics for WGS
    - for Seq2C CNV calling when no capture BED available
    """
    bed = get_all_features(genome)
    return bed\
        .filter(lambda r: r.fields[BedCols.FEATURE] in ['CDS', 'stop_codon'])\
        .filter(high_confidence_filter)\
        .merge()

###############
### ENSEMBL ###
###############
def ensembl_gtf_fpath(genome):
    return _get_ensembl_file(join('gtf', 'ref-transcripts.gtf'), genome.split('-')[0])  # no -alt

def biomart_fpath():
    """ Downloaded from http://www.ensembl.org/biomart/martview/691cdd7b3ea5b7dc64f442ba4ef4a402
    """
    return _get_ensembl_file('mart_export.txt')

def _get_ensembl_file(fname, genome=None):
    if genome:
        return _get(join(genome.split('-')[0], fname), genome)
    else:
        return _get(join(fname))


###################
### TRANSCRIPTS ###
###################
def get_canonical_transcripts_ids(genome):
    short_genome = genome.split('-')[0]
    if short_genome.startswith('GRCh37'):
        short_genome = 'hg19'
    if short_genome.startswith('GRCh38'):
        short_genome = 'hg38'
    check_genome(short_genome)

    canon_fpath = _get(join(pardir, pardir, 'transcripts', 'canon_transcripts_{genome}_ensembl.txt'), genome)
    replacement_fpath = _get(join(pardir, pardir, 'transcripts', 'canon_cancer_replacement.txt'))

    canon_fpath = verify_file(canon_fpath, is_critical=True, description='Canonical transcripts path')
    replacement_fpath = verify_file(replacement_fpath, is_critical=True, description='Canonical cancer transcripts replacement path')

    with open(canon_fpath) as f:
        canon_tx_by_gname = dict(l.strip('\n').split('\t') for l in f)
    with open(replacement_fpath) as f:
        for gname, tx_id in (l.strip('\n').split('\t') for l in f):
            canon_tx_by_gname[gname] = tx_id

    return canon_tx_by_gname


def _get(relative_path, genome=None):
    """
    :param relative_path: relative path of the file inside the repository
    :param genome: genome name. Can contain chromosome name after comma, like hg19-chr20,
                   in case of BED, the returning BedTool will be with added filter.
    :return: BedTools object if it's a BED file, or filepath
    """
    chrom = None
    if genome:
        if '-chr' in genome:
            genome, chrom = genome.split('-')
        check_genome(genome)
        relative_path = relative_path.format(genome=genome)

    path = abspath(join(dirname(__file__), relative_path))
    if not isfile(path) and isfile(path + '.gz'):
        path += '.gz'

    if path.endswith('.bed') or path.endswith('.bed.gz'):
        if path.endswith('.bed.gz'):
            bedtools = which('bedtools')
            if not bedtools:
                critical('bedtools not found in PATH: ' + str(os.environ['PATH']))
            bedtools_v = bedtools_version(bedtools)
            if bedtools_v > 25:
                debug('BED is compressed, creating BedTool')
                bed = BedTool(path)
            else:
                debug('BedTools version is < ' + str(bedtools_v) + ', extracting BED file')
                bed = BedTool(open_gzipsafe(path))
        else:
            debug('BED is uncompressed, creating BedTool')
            bed = BedTool(path)

        if chrom:
            bed = bed.filter(lambda r: r.chrom == chrom)
        return bed
    else:
        return path

def get_hgnc_gene_synonyms():
    return _get('HGNC_gene_synonyms.txt')

def high_confidence_filter(x):
    return x[BedCols.TSL] in ['1', '2', 'NA', '.', None] and x[BedCols.HUGO] not in ['', '.', None]

def get_only_canonical_filter(canon_tx_by_gname):
    return lambda x: x[BedCols.ENSEMBL_ID] == canon_tx_by_gname.get(x[BedCols.HUGO]) or \
                     x[BedCols.ENSEMBL_ID] == canon_tx_by_gname.get(x[BedCols.GENE])

def protein_coding_filter(x):
    return x[BedCols.BIOTYPE] == 'protein_coding'
