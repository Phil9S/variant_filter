rm(list = ls())
args = commandArgs(trailingOnly=TRUE)
###set working directory and import arguments and libraries

setwd(args[1]) ##comment out and set manually if working locally - i.e. non-server side and without bash script
require("stringr")
clock <- as.character(Sys.time())


###default variables & log start
write(clock, file = "R_log.txt", append = FALSE)
write("##Variant Filter Script ## R-script Log - Log Begin (Version 3.0)", file = "R_log.txt", append = TRUE)

##################################################################################################################

###Import data tables from bash script
ad <- read.table("allelicdepth.table", header = TRUE, stringsAsFactors = FALSE)
anno <- read.table("annovar.table", header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote="")
#gq <- read.table("genoqual.table", header = TRUE,stringsAsFactors = FALSE) NOT USED currently
gt <- read.table("genotype.table", header = TRUE, stringsAsFactors = FALSE)
dp <- read.table("sitedepth.table", header = TRUE,stringsAsFactors = FALSE)
config <- read.table("variant_filtering.config", stringsAsFactors = FALSE)
vv <- read.table("variant.table", comment.char = "", header = TRUE, stringsAsFactors = FALSE)
###Read in genetic intolerance lists
rvis <- read.table("RVIS_Unpublished_ExACv2_March2017.tsv", sep="\t", header = TRUE, stringsAsFactors = FALSE, quote="")
gdis <- read.table("GDI_full_10282015.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
###remove columns that aren't needed
ad <- ad[,-(2:5)]
#gq <- gq[,-(2:5)] NOT USED currently
gt <- gt[,-(2:5)]
dp <- dp[,-(2:5)]
vv <- vv[,-(8:10)]
rvis <- rvis[c(1,4)]
gdis <- gdis[c(1,3)]
###rename rs id col to something other than "ID" for safety - Rename cols in genetic intolerance data
names(vv)[4] <- "rsID"
names(rvis) <- c("GENE","RVIS_Pct")
names(gdis) <- c("GENE","GDIS_Phred")
###log starting number of variants - start
varcount <- paste("##Variant Filter Script ## R-script Log - Starting variant count:",nrow(vv))
write(varcount, file = "R_log.txt", append = TRUE)

###Add annotation cols to variant file
vv$GENE <- anno$Gene.refGene
vv$TYPE <- anno$Func.refGene
vv$AA <- anno$AAChange.refGene
vv$CONSEQUENCE <- anno$ExonicFunc.refGene
vv$X1000G <- anno$X1000g2015aug_all
vv$EXAC <- anno$ExAC_ALL
vv$CADD <- anno$CADD_phred
vv$SIFT <- anno$SIFT_pred
vv$POLYPHEN <- anno$Polyphen2_HVAR_pred
vv$CLINVAR <- anno$CLNSIG
###AF - making novel results = 0 and not -9 for freq calculations
vv$X1000G[vv$X1000G == "-9"] <- 0
vv$EXAC[vv$EXAC == "-9"] <- 0
##################################################################################################################

###Damage filters
###import config file settings
X1000g <- config$V2[config$V1 == "G1000"]
exac <- config$V2[config$V1 == "EXAC"]
CADD <- config$V2[config$V1 == "CADD"]
HET <- config$V2[config$V1 == "HET"]
ADEPTH <- config$V2[config$V1 == "ADEPTH"]
QUAL <- config$V2[config$V1 == "QUAL"]
SAMP_NUM <- config$V2[config$V1 == "SAMP_NUM"]

##################################################################################################################

###filter variants to those occuring in exonic regions and splice sites
exonic.ft <- as.data.frame(anno$ID[grepl("^exonic$",anno$Func.refGene)
                                 | grepl("^splicing$",anno$Func.refGene)
                                 | grepl("^exonic;splicing$",anno$Func.refGene)])
###rename col to match other ID cols
names(exonic.ft)[1] <- "ID"
###filter by functional consequence
vv.ft <- vv[vv$ID %in% exonic.ft$ID,]

###log number of variants - func
varcount <- paste("##Variant Filter Script ## R-script Log - Variants falling in exon/splice sites:",nrow(exonic.ft))
write(varcount, file = "R_log.txt", append = TRUE)

##################################################################################################################

###filtering on functional consequnce - ExonicFunction.refGene
func.ft <- as.data.frame(vv.ft$ID[!grepl("^synonymous SNV$",vv.ft$CONSEQUENCE)
                                  & !grepl("^unknown$",vv.ft$CONSEQUENCE)])
                                 
###rename col to match other ID cols
names(func.ft)[1] <- "ID"
###filter by functional consequence
vv.ft <- vv.ft[vv.ft$ID %in% func.ft$ID,]

###log number of variants - func
varcount <- paste("##Variant Filter Script ## R-script Log - Variants with amino acid altering consequence:",nrow(func.ft))
write(varcount, file = "R_log.txt", append = TRUE)

##################################################################################################################

###filter by qual
qual.ft <- as.data.frame(vv.ft$ID[vv.ft$QUAL > QUAL])
names(qual.ft)[1] <- "ID"
vv.ft <- vv.ft[vv.ft$ID %in% qual.ft$ID,]
###log number of variants - func
varcount <- paste("##Variant Filter Script ## R-script Log - Variants matching QUAL:",nrow(qual.ft))
write(varcount, file = "R_log.txt", append = TRUE)

##################################################################################################################

###corce gt data.frame into a matrix and replace non-numeric format with numeric genotypes
###multi-allelic sites are replaced with -2 and missing sites with -9
gt <- gt[gt$ID %in% vv.ft$ID,]
gtm <- as.matrix(gt)
#additional handling of multiallelic sites
#gtm[!gtm == "0/0" | !gtm == "0/1" | !gtm == "1/1" | !gtm == "./."] <- -2
gtm[gtm == "0/0"] <- 0
gtm[gtm == "0/1"] <- 1
gtm[gtm == "1/1"] <- 2
gtm[gtm == "./."] <- -9
gt <- as.data.frame(gtm)
###apply function to sum the number of missing, het, hom and ref sites
refHOM <- apply(gtm,1, function(x) sum(x == 0))
HETp <- apply(gtm,1, function(x) sum(x == 1))
nonHOM <- apply(gtm, 1, function(x) sum(x == 2))
miss <- apply(gtm, 1, function(x) sum(x == -9))
###addition of raw counts of HET/HOM
vv.ft$HET_val <- HETp
vv.ft$HOM_val <- nonHOM
###form matrix for pct calculations
calc <- cbind(refHOM,HETp,nonHOM,miss)

###calculate hetpct & hompct (excluding missing sites) and missingness over all sites
hetpct <- ((calc[,2] / (calc[,1] + calc[,2] + calc[,3]))*100)
hompct <- ((calc[,3] / (calc[,1] + calc[,2] + calc[,3]))*100)
misspct <- ((calc[,4] / length(gt[1,-1])*100))
###append values to new columns in vv.ft 
###same purpose
vv.ft$HET_rate <- hetpct
vv.ft$HOM_rate <- hompct
vv.ft$MISS_rate <- misspct
rm(hetpct, hompct, misspct, calc, HETp, nonHOM, refHOM, miss, gtm)

##################################################################################################################

###filter variant ids that are below values for both 1000g & exac_all
rarity.ft <- as.data.frame(anno$ID[anno$X1000g2015aug_all < X1000g & anno$ExAC_ALL < exac])
names(rarity.ft)[1] <- "ID"

###filter variant ids that are above cadd
cadd.ft <- as.data.frame(anno$ID[anno$CADD_phred > CADD | anno$CADD_phred < 0])
names(cadd.ft)[1] <- "ID"

### filter by rarity
vv.ft <- vv.ft[vv.ft$ID %in% rarity.ft$ID,]

###log number of variants - rarity
varcount <- paste("##Variant Filter Script ## R-script Log - Variants matching rarity threshold:",nrow(vv.ft))
write(varcount, file = "R_log.txt", append = TRUE)

###filter by cadd score
vv.ft <- vv.ft[vv.ft$ID %in% cadd.ft$ID,]

### log number of variant - cadd score
varcount <- paste("##Variant Filter Script ## R-script Log - Variants matching CADD threshold:",nrow(vv.ft))
write(varcount, file = "R_log.txt", append = TRUE)

###tidy variables
rm(rarity.ft,func.ft,qual.ft)

##################################################################################################################

###filter by het/hom ratio and het rate (het > 0 & no het rate in cohort above 15%)
hethom.ft <- as.data.frame(vv.ft$ID[vv.ft$HET_rate < HET])
names(hethom.ft)[1] <- "ID" 
#log number of variants - Het/Hom
varcount <- paste("##Variant Filter Script ## R-script Log - Variants matching Het/Hom thresholds:",nrow(hethom.ft))
write(varcount, file = "R_log.txt", append = TRUE)

###Extract variants based on filtered list
vv.ft <- vv.ft[vv.ft$ID %in% hethom.ft$ID,]
###tidy variables and tables
rm(hethom.ft)

##################################################################################################################

###performing allelic depth transformation to allele percent
###make copy of allelicdepth(ad)
af <- ad[ad$ID %in% vv.ft$ID,]
af[af == "."] <- NA

###Indexing and generation of percent allelic depth info
af_index <- af[1]
af_mat1 <- as.data.frame(apply(af[2:ncol(af)], c(1,2), FUN = function(x) str_split_fixed(x, ",",2)[,1]))
af_mat2 <- as.data.frame(apply(af[2:ncol(af)], c(1,2), FUN = function(x) str_split_fixed(x, ",",2)[,2]))
af_mat1[af_mat1 == ""] <- NA
af_mat2[af_mat2 == ""] <- NA

### conversion to matrix and perform matrix arithematic
af_mat1 <- as.matrix(apply(af_mat1,2,function(x) as.numeric(x)))
af_mat2 <- as.matrix(apply(af_mat2,2,function(x) as.numeric(x)))
ad_pct <- af_mat2 / (af_mat1 + af_mat2)
af <- cbind(af_index, ad_pct)

rm(af_index, ad_pct, af_mat1, af_mat2)

###filter on variants with no af rate above threshold
af.ft <- data.frame(x=rep(0,nrow(af)))
for(i in 1:nrow(af)){
    if(max(af[i,2:ncol(af)], na.rm = TRUE) > ADEPTH){
        af.ft[i,1] <- af[i,1]}
    else{
        af.ft[i,1] <- NA
    }
}
names(af.ft)[1] <- "ID"
af.ft <- subset(af.ft, (!is.na(af.ft[,1])))

vv.ft <- vv.ft[vv.ft$ID %in% af.ft$ID,]
rm(af.ft,i)

###Variants filtered on no variants af above than 0.3
varcount <- paste("##Variant Filter Script ## R-script Log - Variants with at least single alt AD > threshold:",nrow(vv.ft))
write(varcount, file = "R_log.txt", append = TRUE)

##################################################################################################################

###aggregate mutation types and af counts - ONLY PERFORMED ON HET CALLS!
###all HET call gene sums - after functional exonic only - QUAL filtered
agg_all_HET <- aggregate(vv.ft$HET_val, by=list(GENE=vv.ft$GENE), drop = FALSE, FUN=sum)
agg_all_HOM <- aggregate(vv.ft$HOM_val*2, by=list(GENE=vv.ft$GENE), drop = FALSE, FUN=sum)

###Merge HET and HOMO values to single allele count - hom * 2
agg_sum <- merge(agg_all_HET, agg_all_HOM, by="GENE", all = TRUE)
agg_sum[4] <- agg_sum[2] + agg_sum[3]
agg_sum <- agg_sum[,c(-2,-3)]
agg_all <- agg_sum
rm(agg_all_HET,agg_all_HOM,agg_sum)

###Counts for different types of variants - nonsynonymous (excluding splicing) class by either CONSEQUENCE OR TYPE
agg_nsyn_HET <- aggregate(vv.ft$HET_val[vv.ft$CONSEQUENCE == "nonsynonymous SNV" & vv.ft$TYPE == "exonic" ], 
                            by=list(GENE=vv.ft$GENE[vv.ft$CONSEQUENCE == "nonsynonymous SNV" & vv.ft$TYPE == "exonic"]), drop = FALSE, FUN=sum)

agg_nsyn_HOM <- aggregate(vv.ft$HOM_val[vv.ft$CONSEQUENCE == "nonsynonymous SNV" & vv.ft$TYPE == "exonic" ]*2, 
                            by=list(GENE=vv.ft$GENE[vv.ft$CONSEQUENCE == "nonsynonymous SNV" & vv.ft$TYPE == "exonic"]), drop = FALSE, FUN=sum)
###Merge HET and HOMO values to single allele count - hom * 2
agg_nsyn_sum <- merge(agg_nsyn_HET, agg_nsyn_HOM, by="GENE", all = TRUE)
agg_nsyn_sum[4] <- agg_nsyn_sum[2] + agg_nsyn_sum[3]
agg_nsyn_sum <- agg_nsyn_sum[,c(-2,-3)]
agg_nsyn <- agg_nsyn_sum
rm(agg_nsyn_HET,agg_nsyn_HOM,agg_nsyn_sum)

###Counts for different types of variants - Truncating (including stop loss) class by either CONSEQUENCE OR TYPE
agg_trunc_HET <- aggregate(vv.ft$HET_val[vv.ft$CONSEQUENCE == "stopgain" | vv.ft$CONSEQUENCE == "stoploss" 
                                           | vv.ft$CONSEQUENCE == "frameshift deletion" 
                                           | vv.ft$CONSEQUENCE == "frameshift insertion"],
                             by=list(GENE=vv.ft$GENE[vv.ft$CONSEQUENCE == "stopgain" | vv.ft$CONSEQUENCE == "stoploss" 
                                           | vv.ft$CONSEQUENCE == "frameshift deletion" 
                                           | vv.ft$CONSEQUENCE == "frameshift insertion"]),drop = FALSE, FUN=sum)

agg_trunc_HOM <- aggregate(vv.ft$HOM_val[vv.ft$CONSEQUENCE == "stopgain" | vv.ft$CONSEQUENCE == "frameshift deletion" 
                                           | vv.ft$CONSEQUENCE == "frameshift insertion"]*2,
                             by=list(GENE=vv.ft$GENE[vv.ft$CONSEQUENCE == "stopgain" | vv.ft$CONSEQUENCE == "frameshift deletion" 
                                                     | vv.ft$CONSEQUENCE == "frameshift insertion"]), drop = FALSE, FUN=sum)
###Merge HET and HOMO values to single allele count - hom * 2
agg_trunc_sum <- merge(agg_trunc_HET, agg_trunc_HOM, by="GENE", all = TRUE)
agg_trunc_sum[4] <- agg_trunc_sum[2] + agg_trunc_sum[3]
agg_trunc_sum <- agg_trunc_sum[,c(-2,-3)]
agg_trunc <- agg_trunc_sum
rm(agg_trunc_HET,agg_trunc_HOM,agg_trunc_sum)

###Counts for different types of variants - OTHER - NON frameshifting/stoploss
agg_other_HET <- aggregate(vv.ft$HET_val[vv.ft$CONSEQUENCE == "stoploss" | vv.ft$CONSEQUENCE == "nonframeshift deletion" 
                                        | vv.ft$CONSEQUENCE == "nonframeshift insertion"],
                          by=list(GENE=vv.ft$GENE[vv.ft$CONSEQUENCE == "stoploss" | vv.ft$CONSEQUENCE == "nonframeshift deletion" 
                                                  | vv.ft$CONSEQUENCE == "nonframeshift insertion"]),drop = FALSE, FUN=sum)

agg_other_HOM <- aggregate(vv.ft$HOM_val[vv.ft$CONSEQUENCE == "stoploss" | vv.ft$CONSEQUENCE == "nonframeshift deletion" 
                                         | vv.ft$CONSEQUENCE == "nonframeshift insertion"]*2,
                           by=list(GENE=vv.ft$GENE[vv.ft$CONSEQUENCE == "stoploss" | vv.ft$CONSEQUENCE == "nonframeshift deletion" 
                                                   | vv.ft$CONSEQUENCE == "nonframeshift insertion"]), drop = FALSE, FUN=sum)
###Merge HET and HOMO values to single allele count - hom * 2
agg_other_sum <- merge(agg_other_HET, agg_other_HOM, by="GENE", all = TRUE)
agg_other_sum[4] <- agg_other_sum[2] + agg_other_sum[3]
agg_other_sum <- agg_other_sum[,c(-2,-3)]
agg_other <- agg_other_sum
rm(agg_other_HET,agg_other_HOM,agg_other_sum)

##Counts for different types of variants - splicing class by either CONSEQUENCE OR TYPE
if("-9" %in% vv.ft$CONSEQUENCE){
  agg_splc_HET <- aggregate(vv.ft$HET_val[vv.ft$CONSEQUENCE == "-9" | vv.ft$TYPE == "splicing" | vv.ft$TYPE == "exonic;splicing" ], 
                    by=list(GENE=vv.ft$GENE[vv.ft$CONSEQUENCE == "-9" | vv.ft$TYPE == "splicing" | vv.ft$TYPE == "exonic;splicing" ]), 
                    drop = FALSE, FUN=sum)

  agg_splc_HOM <- aggregate(vv.ft$HOM_val[vv.ft$CONSEQUENCE == "-9" | vv.ft$TYPE == "splicing" | vv.ft$TYPE == "exonic;splicing" ]*2, 
                            by=list(GENE=vv.ft$GENE[vv.ft$CONSEQUENCE == "-9" | vv.ft$TYPE == "splicing" | vv.ft$TYPE == "exonic;splicing" ]), 
                            drop = FALSE, FUN=sum)
  ###Merge HET and HOMO values to single allele count - hom * 2
  agg_splc_sum <- merge(agg_splc_HET, agg_splc_HOM, by="GENE", all = TRUE)
  agg_splc_sum[4] <- agg_splc_sum[2] + agg_splc_sum[3]
  agg_splc_sum <- agg_splc_sum[,c(-2,-3)]
  agg_splc <- agg_splc_sum
  rm(agg_splc_HET,agg_splc_HOM,agg_splc_sum)
}

###merge splicing and OTHER counts - provided splice sites were found
if(exists("agg_splc")){
  agg_splcother_sum <- merge(agg_splc, agg_other, by="GENE", all = TRUE)
  agg_splcother_sum[4] <- agg_splcother_sum[2] + agg_splcother_sum[3]
  agg_splcother_sum <- agg_splcother_sum[,c(-2,-3)]
  agg_other <- agg_splcother_sum
  rm(agg_splcother_sum)
}

###Gene count merging to single data.frame - replace na with 0
agg_allnsyn <- merge(agg_all, agg_nsyn, by = "GENE", all = TRUE)
names(agg_allnsyn)[2] <- "AC_all"
names(agg_allnsyn)[3] <- "AC_nsyn"
agg_alltruncnsyn <- merge(agg_allnsyn, agg_trunc, by = "GENE", all = TRUE)
names(agg_alltruncnsyn)[4] <- "AC_trnc"
agg_alltruncnsyn <- merge(agg_alltruncnsyn, agg_other, by = "GENE", all = TRUE)
names(agg_alltruncnsyn)[5] <- "AC_other"
agg_alltruncnsyn[is.na(agg_alltruncnsyn)] <- 0

rm(agg_all,agg_allnsyn,agg_trunc,agg_nsyn,agg_other)

##report number of lines in which the variant counter doesn't report a matching value between all and the sub columns
allelecount_false_sum <- length(grep("FALSE", agg_alltruncnsyn[2] == agg_alltruncnsyn[3] + agg_alltruncnsyn[4] + agg_alltruncnsyn[5]))
varcount <- paste("##Variant Filter Script ## R-script Log - Number of rows with non-matching allele counts:",allelecount_false_sum)
write(varcount, file = "R_log.txt", append = TRUE)

###merge the gene counts into the original vv.ft file - adding 4 columns - write table as output - rm variables
###totals of subsets may not all sum to all due to certain combinations (e.g. exonic;splicing stoploss)
vv.ft <- merge(vv.ft, agg_alltruncnsyn, by = 'GENE', sort = FALSE)[,union(names(vv.ft), names(agg_alltruncnsyn))]
write.table(agg_alltruncnsyn,file = "variant_filtering_GeneAC.tsv", quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
rm(agg_alltruncnsyn)

##################################################################################################################
### aggregate AF per gene using rare variants - truncating and nonsynonymous
### TRUNCATING aggregation - X1000G
AF_gene_truncX1000 <- aggregate(vv.ft$X1000G[vv.ft$CONSEQUENCE == "stopgain" | vv.ft$CONSEQUENCE == "frameshift deletion" 
                        | vv.ft$CONSEQUENCE == "frameshift insertion"],
                   by=list(GENE=vv.ft$GENE[vv.ft$CONSEQUENCE == "stopgain" | vv.ft$CONSEQUENCE == "frameshift deletion" 
                        | vv.ft$CONSEQUENCE == "frameshift insertion"]),drop = FALSE, FUN=sum)
### TRUNCATING aggregation - EXAC
AF_gene_truncEXAC <- aggregate(vv.ft$EXAC[vv.ft$CONSEQUENCE == "stopgain" | vv.ft$CONSEQUENCE == "frameshift deletion" 
                                             | vv.ft$CONSEQUENCE == "frameshift insertion"],
                                by=list(GENE=vv.ft$GENE[vv.ft$CONSEQUENCE == "stopgain" | vv.ft$CONSEQUENCE == "frameshift deletion" 
                                                        | vv.ft$CONSEQUENCE == "frameshift insertion"]),drop = FALSE, FUN=sum)
###Merge and take max of EXAC & 1000G
AF_gene_trunc <- merge(AF_gene_truncX1000, AF_gene_truncEXAC, by = "GENE", all = TRUE)
AF_gene_trunc$AggAF_Trunc <- apply(AF_gene_trunc[2:3], 1, function(x) max(x))
AF_gene_trunc <- AF_gene_trunc[-(2:3)]


### NONSYNONYMOUS aggregation - X1000G
AF_gene_nsynX1000 <- aggregate(vv.ft$X1000G[vv.ft$CONSEQUENCE == "nonsynonymous SNV"],
                                by=list(GENE=vv.ft$GENE[vv.ft$CONSEQUENCE == "nonsynonymous SNV"]),drop = FALSE, FUN=sum)
### NONSYNONYMOUS aggregation - EXAC
AF_gene_nsynEXAC <- aggregate(vv.ft$EXAC[vv.ft$CONSEQUENCE == "nonsynonymous SNV"],
                               by=list(GENE=vv.ft$GENE[vv.ft$CONSEQUENCE == "nonsynonymous SNV"]),drop = FALSE, FUN=sum)
###Merge and take max of EXAC & 1000G
AF_gene_nsyn <- merge(AF_gene_nsynX1000, AF_gene_nsynEXAC, by = "GENE", all = TRUE)
AF_gene_nsyn$AggAF_nsyn <- apply(AF_gene_nsyn[2:3], 1, function(x) max(x))
AF_gene_nsyn <- AF_gene_nsyn[-(2:3)]

###Merging to main dataframe and removing NA values
AF_gene <- merge(AF_gene_trunc, AF_gene_nsyn, by = "GENE", all = TRUE)
AF_gene[is.na(AF_gene)] <- -9
vv.ft <- merge(vv.ft, AF_gene, by = "GENE", all.x = TRUE)[, union(names(vv.ft), names(AF_gene))]
rm(AF_gene, AF_gene_nsyn, AF_gene_nsynEXAC, AF_gene_nsynX1000, AF_gene_trunc, AF_gene_truncEXAC, AF_gene_truncX1000)

##################################################################################################################

###log number of variants - Final
clock <- as.character(Sys.time())
write(clock, file = "R_log.txt", append = TRUE)
varcount <- paste("##Variant Filter Script ## R-script Log - Final variants passing filters:",nrow(vv.ft))
write(varcount, file = "R_log.txt", append = TRUE)

###Adding genetic intolerance scores
vv.ft <- merge(vv.ft, rvis, by = "GENE", all.x = TRUE)[, union(names(vv.ft), names(rvis))]
vv.ft <- merge(vv.ft, gdis, by = "GENE", all.x = TRUE)[, union(names(vv.ft), names(gdis))]
vv.ft[is.na(vv.ft)] <- -9
###Add genotype information for remaining variants
if(ncol(gt) > SAMP_NUM){
  gt.ft <- gt[gt$ID %in% vv.ft$ID,]
  vvgt <- vv.ft
  clock <- as.character(Sys.time())
  write(clock, file = "R_log.txt", append = TRUE)
  varcount <- paste("##Variant Filter Script ## R-script Log - Seperate genotype file generated",nrow(vv.ft))
  write(varcount, file = "R_log.txt", append = TRUE)
  names(gt.ft)[1] <- "Id"
  write.table(gt.ft,file = "variant_filtering_results_GT.tsv",sep = "\t",row.names = FALSE, quote = FALSE)
} else {
  gt.ft <- gt[gt$ID %in% vv.ft$ID,]
  vvgt <- merge(vv.ft,gt.ft, sort = FALSE)
}
###rename ID col - issues with opening files in excel with "ID" as the first value
names(vvgt)[1] <- "Id"
names(af)[1] <- "Id"

##################################################################################################################
### Col trimming for final tables
drop_col <- c("HET_rate","HOM_rate")
vvgt <- vvgt[ , !(names(vvgt) %in% drop_col)]

##################################################################################################################
###Generating biallelic - potential compound het calls
#form empty vector
bi_allelic_var <- character()
##for each unique gene, identify samples where number of non-ref calls > 1 i.e. two or more het calls for one gene in one sample
for(i in unique(vvgt$GENE)){
  bi_gene <- vvgt[vvgt$GENE == i,]
  bi_gene[26:ncol(bi_gene)] <- sapply(bi_gene[26:ncol(bi_gene)],FUN = function(x) as.numeric(as.character(x)))
  bi_gene[bi_gene == "-9"] <- NA
  bi_count <- as.data.frame(apply(bi_gene[26:ncol(bi_gene)],2, function(x) sum(as.numeric(x),na.rm = TRUE)))
  bi_count$names <- row.names(bi_count)
  sample_list <- bi_count$names[bi_count[1] > 1]
  ## above generates a list of samples with two or more het calls per gene
  ## secondary loop finds the variant Id numbers for the variants called in each sample in sample_list per gene
  for(s in sample_list){
    p <- bi_gene$Id[bi_gene[s] == 1]
    p <- p[!is.na(p)]
    bi_allelic_var <- append(bi_allelic_var,p)
  ## appends each successive loop of variant Ids that are potentially biallelic/compound het
  }
}
### select filtered ids from the main output dataframe
vvgt_bi <- vvgt[vvgt$Id %in% bi_allelic_var,]

###write filtered table out
write.table(vvgt,file = "variant_filtering_results.tsv",sep = "\t",row.names = FALSE,quote = FALSE)
write.table(af,file = "variant_filtering_results_AD.tsv",sep = "\t",row.names = FALSE,quote = FALSE)
write.table(vvgt_bi,file = "variant_comphet-biallelic_results.tsv",sep = "\t",row.names = FALSE,quote = FALSE)
