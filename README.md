# Variant Filtering Script - R & Bash

A set of scripts for filtering multi-sample VCF files to generate candidate variant lists

## Main command:

```bash
./variant_filtering.sh --project my_project \
--input /data/Alignments/your_merged_vcf.vcf \
--output /home/user/ \
--config variant_filtering.config \
--temp TRUE
```

## 1. variant_filtering_results.tsv

#### **Ouput description table** 
| Column Name | Description  
| :---        |:--- 
| Id | Variant Id Number  
|CHROM | Chromosome number
|POS | Variant loci position
|rsID | rsID or position again if not known
|REF | Reference base
|ALT | Alternate base
|QUAL | QUAL score
|GENE | Gene symbol
|TYPE | Exonic/splicing - what type of region
|AA | AA change for multiple transcripts - look up appropriate one
|CONSEQUENCE | Functional change or consequence (stop gain etc)
|X1000G | Rarity in 1000 genomes
|EXAC | Rarity in ExAC consortium
|CADD | CADD score variant damage prediction
|HET_val | Number of samples called with a het alt call at this position
|HOM_val | Number of samples called with a hom alt call at this position
|MISS_rate | The pecentage rate of sites failing to genotype (missing rate)
|AC_all | Alt allelic counts for this site (post filtering) - homozygous = 2 towards count
|AC_nsyn | Same as AC_all but only for nonsynonymous changes
|AC_trnc | Same as AC_all but only for truncating changes (stop gain, frameshifts)
|AC_other | Same as AC_all but only for all other types of changes
|AggAF_Trunc | Aggregated allele frequency of truncating mutations based on 1000Genomes or ExAc
|AggAF_nsyn | Same as AggAF_Trunc but for nonsynonymous mutations
|RVIS_Pct | Residual variant intolerance score percentile - Higher is more tolerant
|GDIS_Phred | Genetic damage intolerance score - Phred scaled - Higher is more tolerant
|Sample_Id | Columns from this point are genotyping columns - 0 = ref, 1 = het, 2 = hom, -9 = missing


## 2. variant_filtering_GeneAC.tsv

* Seperate file containing only the AC_* columns from variant_filtering_results.tsv
* Used for statistical analysis of case/controls for all genes present in the variant file

## 3. variant_filtering_results_AD.tsv

* Allelic depth ratio matrix
* describes the proportion of reads supporting the REF against ALT
* Calculates as ALT READs / REF reads + ALT reads
* Used as a filtering step in variant_filtering.config file - variant retained if one sample passes the threshold
* Multi-allelic variants are omitted from this filter for ease of parsing the data

## 4. R_log.txt

* Log file of each filtering step in the R filtering script
* Denotes start times and variant counts at each stage
* Declares the number of multi-allelic variants carried forward
* Describes the number of variants in which AC_all was not equal to sum(AC_trunc + AC_nsyn + AC_other)
   Should be a very low proportion of total variants and is due to unusual variant type/consequence combinations

## 5. variantfilter.log

* Log file for the bash portion of the variant filtering script
* Contains the parameters used and declares which thresholds were specified in the variant_filtering.config file
* Mainly just for error checks and confirming the script has completed

##

### Required dependencies and software:
* [Annovar](http://annovar.openbioinformatics.org/en/latest/) - Including various databases
* [vcftools](http://vcftools.sourceforge.net/index.html)
* [bcftools](https://samtools.github.io/bcftools/bcftools.html)
* [R (version 3.3.3)](https://cran.r-project.org/)
# Required R Libraries:
* [stringr](https://cran.r-project.org/web/packages/stringr/index.html)     

##
