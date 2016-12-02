import os
from os.path import dirname, join, abspath, isfile, basename, splitext
import sys

from ngs_utils.bedtools import BedTool
from ngs_utils.bed_utils import bedtools_version
from ngs_utils.file_utils import which, open_gzipsafe
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

##############
### REFSEQ ###
##############
REFSEQ_DIR = 'RefSeq'

def refseq_knowngene_fpath(genome):
    return _get_refseq_file('RefSeq_knownGene.txt', genome.split('-')[0])  # no -alt

def get_refseq_gene(genome):
    return _get_refseq_file('refGene.txt.gz', genome.split('-')[0])  # no -alt

def _get_refseq_file(fname, genome=None):
    return _get(join(REFSEQ_DIR, genome, fname), genome)

def get_refseq_dirpath():
    return abspath(join(dirname(__file__), REFSEQ_DIR))

###############
### ENSEMBL ###
###############
ENSEMBL_DIR = 'Ensembl'

def ensembl_gtf_fpath(genome):
    return _get_ensembl_file(join('gtf', 'ref-transcripts.gtf'), genome.split('-')[0])  # no -alt

def biomart_fpath(genome):
    return _get_ensembl_file('biomart.tsv')

def _get_ensembl_file(fname, genome=None):
    if genome:
        return _get(join(ENSEMBL_DIR, genome.split('-')[0], fname), genome)
    else:
        return _get(join(ENSEMBL_DIR, fname))


'''
This repository is made for storing genomic features coordinates and annotations.

Annotation in TargQC:
  - prepocess: generate BED file from ref-transcripts.gtf.db:
    - annotate with RefSeq IDs and TSL

  - annotate by priority:
    - use all RefSeq exons (TODO: check if need to separate processed_transcript)
    - use all RefSeq transcripts
    - use other exons
    - use other transcripts

Annotation in bcbio:
  - generate BED file from ref-transcripts.gtf.db, containing all features
    - annotate with RefSeq IDs and TSL

  - annotate by priority:



- BED files annotation.
  Stores known features:
   - Feature: Exon, CDS, Transcript, Gene
   - Biotype: protein-coding, RNA
   - TSL
   - Transcript ID in Ensembl and RefSeq

  Annotation priority:
   - known protein_coding CDS|stop_codon|ncRNA_exon
   - known protein_coding UTR
   - known protein_coding transcript (annotate as "intron")
   - predicted protein_coding CDS|stop_codon|ncRNA_exon
   - predicted protein_coding UTR
   - predicted protein_coding transcript (annotate as "intron")

   -

  After annotation, write all features used for annotation and make a file for TargQC coverage reports.

- CDS BED file for Seq2C, if capture panel is not known.
  Ensembl-based, contains
   - Known canonical CDS
'''


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

def protein_coding_filter(x):
    return x[BedCols.BIOTYPE] == 'protein_coding'
