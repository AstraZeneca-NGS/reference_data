# Reference data

## Capture region BED files

Collects commonly used capture region BED files. These are installed and available for
use in [bcbio](https://github.com/chapmanb/bcbio-nextgen) analyses. Includes
files for hg19 (chr1, chr2, chr3... style naming) and GRCh37 (1, 2, 3... style naming).

- capture_regions/Exome-AZ_V2 -- Inclusive file of exome regions, hand curated by AstraZeneca Oncology.
- capture_regions/Exome-IDT_V1 -- [IDT xGen Exome Research Panel v1.0](http://www.idtdna.com/pages/products/nextgen/target-capture/xgen-lockdown-panels/xgen-exome-panel)
- capture_regions/Exome-Agilent_V6 -- [Agilent SureSelect Human All Exon V6](http://www.genomics.agilent.com/article.jsp?crumbAction=push&pageId=9500016)
- capture_regions/Exome-Agilent_V5_Clinical -- [Agilent SureSelect Clinical Research Exome](http://www.genomics.agilent.com/article.jsp?crumbAction=push&pageId=4600005) 
- capture_regions/Exome-MedExome -- [Nimblegen SeqCap EZ MedExome](http://sequencing.roche.com/products/nimblegen-seqcap-target-enrichment/seqcap-ez-system/seqcap-ez-medexome.html)
- capture_regions/Exome-NGv3 -- [Nimblegen SeqCap EZ Exome v3](http://sequencing.roche.com/products/nimblegen-seqcap-target-enrichment/seqcap-ez-system/seqcap-ez-exome-v3.html)

## Canonical transcripts

`transcripts/cancer_transcripts_*_ensembl.txt` contain the lists of canonical (longest) transcript IDs that are used by [SnpEff](http://snpeff.sourceforge.net/) variant prediction tool when it run with the `-canon` flag on (only in Ensembl-based versions of reference databases GRCh37.** and GRCh38.** in SnpEff notation). However, not all IDs in the list represent the most cancer-relevant isoforms, and `transcripts/canon_cancer_replacement.txt` provides a map of transcripts for replacement with the `-canonList` option:
```
java -jar snpEff.jar GRCh37.75 test.vcf -canon -canonList transcripts/canon_cancer_replacement.txt
```

To use canonical transcripts for variant annotation in bcbio, add the following into config:
```	
algorithm:
  effects_transcripts: canon
```
To use cancer transcripts, use the following:
```	
algorithm:
  effects_transcripts: canonical_cancer
```

## BED files annotation

### Installation
```
cd bed_annotation
python setup.py install
```

### Usage
```
bed_annotation/annotate_bed.py regions.bed -g hg19 -o regions.ann.bed
``` 

Script checks each region against Ensembl genomic features database, and writes a BED file in a standardized format with a gene symbol, strand and exon rank in 4-6th columns:

`regions.bed`:
```
chr1    69090   70008
chr1    367658  368597
```

`regions.ann.bed`:
```
chr1    69090   70008   OR4F5   1       +
chr1    367658  368597  OR4F29  1       +
```

The piority for choosing transcripts for annotation is the following:
    - Biotype: protein coding > others > RNA > decay > sense_* > antisense > translated_* > transcribed_*
    - TSL: 1 > NA > others > 2 > 3 > 4 > 5
    - Has HUGO gene symbold
    - Overlap % with transcript
    - Overlap % with CDS
    - Overlap % with exons
    - Is cancer canonical
    - Transcript length

Use `--extended` option to report columns with details on overlapping transcripts, and overlap sizes:
```
annotate_bed.py regions.bed -g hg19 -o regions.ann.bed --extended
```

`regions.ann.bed`:
```
## Tx_overlap_%: part of region overlapping with transcripts
## Exon_overlaps_%: part of region overlapping with exons
## CDS_overlaps_%: part of region overlapping with protein coding regions
#Chrom  Start   End     Gene    Exon    Strand  Feature Biotype Ensembl_ID      TSL     HUGO    Tx_overlap_%    Exon_overlaps_% CDS_overlaps_%
chr1    69090   70008   OR4F5   1       +       capture protein_coding  ENST00000335137 NA      OR4F5   100.0   100.0   99.7
chr1    367658  368597  OR4F29  1       +       capture protein_coding  ENST00000426406 NA      OR4F29  100.0   100.0   99.7
```

Use `--report-all-good-annotations` flag to report all reasonable overlaps for each regions (might produce multiple lines per region, but makes sense for long regions)
```
annotate_bed.py regions.bed -g hg19 -o regions.ann.bed --extended --report-all-good-annotations
```
`regions.ann.bed`:
```
## Tx_overlap_%: part of region overlapping with transcripts
## Exon_overlaps_%: part of region overlapping with exons
## CDS_overlaps_%: part of region overlapping with protein coding regions
#Chrom  Start   End     Gene    Exon    Strand  Feature Biotype Ensembl_ID      TSL     HUGO    Tx_overlap_%    Exon_overlaps_% CDS_overlaps_%
chr1    69090   70008   OR4F5   1       +       capture protein_coding  ENST00000335137 NA      OR4F5   100.0   100.0   99.7
chr1    367658  368597  OR4F29  1       +       capture protein_coding  ENST00000426406 NA      OR4F29  100.0   100.0   99.7
chr1    367658  368597  OR4F29  1       +       capture protein_coding  ENST00000412321 NA      OR4F29  100.0   100.0   82.1
```