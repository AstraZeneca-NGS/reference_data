# Reference data

## Capture region BED files

Collects commonly used capture region BED files. These are installed and available for
use in [bcbio](https://github.com/chapmanb/bcbio-nextgen) analyses. Includes
files for hg19 (chr1, chr2, chr3... style naming) and GRCh37 (1, 2, 3... style naming).

- `capture_regions/Exome-AZ_V2` -- Inclusive file of exome regions, hand curated by AstraZeneca Oncology.
- `capture_regions/Exome-IDT_V1` -- [IDT xGen Exome Research Panel v1.0](http://www.idtdna.com/pages/products/nextgen/target-capture/xgen-lockdown-panels/xgen-exome-panel)
- `capture_regions/Exome-Agilent_V6` -- [Agilent SureSelect Human All Exon V6](http://www.genomics.agilent.com/article.jsp?crumbAction=push&pageId=9500016)
- `capture_regions/Exome-Agilent_V5_Clinical` -- [Agilent SureSelect Clinical Research Exome](http://www.genomics.agilent.com/article.jsp?crumbAction=push&pageId=4600005) 
- `capture_regions/Exome-MedExome` -- [Nimblegen SeqCap EZ MedExome](http://sequencing.roche.com/products/nimblegen-seqcap-target-enrichment/seqcap-ez-system/seqcap-ez-medexome.html)
- `capture_regions/Exome-NGv3` -- [Nimblegen SeqCap EZ Exome v3](http://sequencing.roche.com/products/nimblegen-seqcap-target-enrichment/seqcap-ez-system/seqcap-ez-exome-v3.html)

## Canonical transcripts

Files under `transcripts/cancer_transcripts_*_ensembl.txt` contain IDs of canonical (longest) transcripts that are used by [SnpEff](http://snpeff.sourceforge.net/) variant prediction tool when it run with the `-canon` flag (only in Ensembl-based versions of reference databases GRCh37.** and GRCh38.** in SnpEff notation). Since not all IDs in the list represent the most cancer-relevant isoforms, `transcripts/canon_cancer_replacement.txt` provides a map of transcripts for replacement with the `-canonList` option:
```
java -jar snpEff.jar GRCh37.75 test.vcf -canon -canonList transcripts/canon_cancer_replacement.txt
```

To use the canonical transcripts for variant annotation in bcbio, add the following into your configuration YAML file:
```	
algorithm:
  effects_transcripts: canon
```
To use the cancer transcripts, use the following:
```	
algorithm:
  effects_transcripts: canonical_cancer
```

## BED file annotation

The `bed_annotation` package provides tools for annotation of BED file regions with gene symbols, based on Ensembl data.

### Installation
```
cd bed_annotation
python setup.py install
```

### Usage
```
annotate_bed.py INPUT.bed -g hg19 -o OUTPUT.bed
``` 

Script checks each region against the Ensembl genomic features database, and writes a BED file in a standardized format with a gene symbol, strand and exon rank in 4-6th columns:

`INPUT.bed`:
```
chr1    69090   70008
chr1    367658  368597
```

`OUTPUT.bed`:
```
chr1    69090   70008   OR4F5   1       +
chr1    367658  368597  OR4F29  1       +
```

#### Priority

The piority for choosing transcripts for annotation is the following:
- Overlap % with transcript
- Overlap % with CDS
- Overlap % with exons
- Biotype (`protein_coding` > others > `*RNA` > `*_decay` > `sense_*` > `antisense` > `translated_*` > `transcribed_*`)
- TSL (1 > NA > others > 2 > 3 > 4 > 5)
- Presence of a HUGO gene symbol
- Is cancer canonical
- Transcript size

#### Extended annotation

Use `--extended` option to report extra columns with details on features, biotype, overlapping transcripts and overlap sizes:
```
annotate_bed.py INPUT.bed -g hg19 -o OUTPUT.bed --extended
```

`OUTPUT.bed`:
```
## Tx_overlap_%: part of region overlapping with transcripts
## Exon_overlaps_%: part of region overlapping with exons
## CDS_overlaps_%: part of region overlapping with protein coding regions
#Chrom  Start   End     Gene    Exon    Strand  Feature Biotype Ensembl_ID      TSL     HUGO    Tx_overlap_%    Exon_overlaps_% CDS_overlaps_%
chr1    69090   70008   OR4F5   1       +       capture protein_coding  ENST00000335137 NA      OR4F5   100.0   100.0   99.7
chr1    367658  368597  OR4F29  1       +       capture protein_coding  ENST00000426406 NA      OR4F29  100.0   100.0   99.7
```

#### Ambuguous annotations

Regions may overlap mltiple genes. The `--ambiguities` controls how the script resolves such ambiguities
- `--ambiguities all` -- report all reliable overlaps (in order in the "priority" section, see above)
- `--ambiguities all_ask` -- stop execution and ask user which annotation to pick
- `--ambiguities best_all` -- find the best overlap, and if there are several equally good, report all (in terms of the "priority" above)
- `--ambiguities best_ask` -- find the best overlap, and if there are several equally good, ask user
- `--ambiguities best_one` -- find the best overlap, and if there are several equally good, report any of them

Note that the first 4 options might output multiple lines per region, e.g.:
```
annotate_bed.py INPUT.bed -g hg19 -o OUTPUT.bed --extended --ambiguities best_all
```
`OUTPUT.bed`:
```
## Tx_overlap_%: part of region overlapping with transcripts
## Exon_overlaps_%: part of region overlapping with exons
## CDS_overlaps_%: part of region overlapping with protein coding regions
#Chrom  Start   End     Gene    Exon    Strand  Feature Biotype Ensembl_ID      TSL     HUGO    Tx_overlap_%    Exon_overlaps_% CDS_overlaps_%
chr1    69090   70008   OR4F5   1       +       capture protein_coding  ENST00000335137 NA      OR4F5   100.0   100.0   100.0
chr1    367658  368597  OR4F29  1       +       capture protein_coding  ENST00000426406 NA      OR4F29  100.0   100.0   100.0
chr1    367658  368597  OR4F29  1       +       capture protein_coding  ENST00000412321 NA      OR4F29  100.0   100.0   100.0
```