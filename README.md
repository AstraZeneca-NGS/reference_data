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

## Transcripts

`transcripts/cancer_transcripts_*_ensembl.txt` contain the lists of canonical (longest) transcript IDs that are used by [SnpEff](http://snpeff.sourceforge.net/) variant prediction tool when it run with the `-canon` flag on (only in Ensembl-based versions of reference databases GRCh37.75 and GRCh38.82 in SnpEff notation). However, not all IDs in the list represent the most cancer-relevant isoforms. `transcripts/canon_cancer_replacement.txt` provides a map of transcripts for replacement with the `-canonList` option:
```
java -jar snpEff.jar GRCh37.75 test.vcf -canon -canonList transcripts/canon_cancer_replacement.txt
```

## BED file annotation

```
bed_annotation/annotate_bed.py file.bed -g hg19 -o file.anno.bed
``` 
Script checks each region against Ensembl GTF file and annotates 
with a gene symbol, strand and exon number.

Priority for choosing transcripts to annotate:
    Canonical

    biotype_key = 1
    biotype_rank = ['protein_coding', '__default__', 'rna', 'decay', 'sense_', 'antisense', 'translated_', 'transcribed_']
    for key in biotype_rank:
        if key in x[ensembl.BedCols.BIOTYPE].lower():
            biotype_key = biotype_rank.index(key)
    tsl_key = {'1': 0, '2': 2, '3': 3, '4': 4, '5': 5}.get(x[ensembl.BedCols.TSL], 1)
    hugo_key = 0 if x[ensembl.BedCols.HUGO] not in ['.', '', None] else 1

    overlap_key = 0
    if len(x) > ensembl.BedCols.TX_OVERLAP_PERCENTAGE:
        overlap_key = [(-x[key] if x[key] is not None else 0)
                       for key in [ensembl.BedCols.TX_OVERLAP_PERCENTAGE,
                                   ensembl.BedCols.CDS_OVERLAPS_PERCENTAGE,
                                   ensembl.BedCols.EXON_OVERLAPS_PERCENTAGE]]

    length_key = -(int(x[ensembl.BedCols.END]) - int(x[ensembl.BedCols.START]))

    return canonical_key, biotype_key, tsl_key, hugo_key, overlap_key, length_key
    
    

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
