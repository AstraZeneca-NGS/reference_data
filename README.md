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

The full list of genes with replaced transcripts:
```
AKT1     ENST00000555528
BRCA1    ENST00000357654
CD79B    ENST00000006750
CDKN2A   ENST00000304494
CHEK1    ENST00000534070
CHEK2    ENST00000328354
ESR1     ENST00000206249
FANCL    ENST00000233741
FGFR1    ENST00000447712
FGFR2    ENST00000457416
FGFR3    ENST00000440486
MET      ENST00000397752
MYD88    ENST00000396334
PPP2R2A  ENST00000380737
RAD51D   ENST00000345365
RAD54L   ENST00000371975
GNAS     ENST00000371085
TP53     ENST00000269305
ARID1B   ENST00000350026
TET2     ENST00000380013
CEBPA    ENST00000498907
PIK3C2G  ENST00000538779
```

