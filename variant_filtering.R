rm(list = ls())
args = commandArgs(trailingOnly=TRUE)
###set working directory and import arguments and libraries

setwd(args[1]) ##comment out and set manually if working locally - i.e. non-server side and without bash script
require("stringr")
clock <- as.character(Sys.time())


###default variables & log start
write(clock, file = "R_log.txt", append = FALSE)
write("##Variant Filter Script ## R-script Log - Log Begin (Version 2.2)", file = "R_log.txt", append = TRUE)

##################################################################################################################

###Import data tables from bash script
ad <- read.table("allelicdepth.table", header = TRUE,stringsAsFactors = FALSE)
anno <- read.table("annovar.table", header = TRUE, sep = "\t",stringsAsFactors = FALSE, quote="")
#gq <- read.table("genoqual.table", header = TRUE,stringsAsFactors = FALSE) NOT USED currently
gt <- read.table("genotype.table", header = TRUE,stringsAsFactors = FALSE)
dp <- read.table("sitedepth.table", header = TRUE,stringsAsFactors = FALSE)
config <- read.table("variant_filtering.config",stringsAsFactors = FALSE)
vv <- read.table("variant.table",comment.char = "",header = TRUE,stringsAsFactors = FALSE)
###remove columns that aren't needed
ad <- ad[,-(2:5)]
#gq <- gq[,-(2:5)] NOT USED currently
gt <- gt[,-(2:5)]
dp <- dp[,-(2:5)]
vv <- vv[,-(8:10)]
###rename rs id col to something other than "ID" for safety
names(vv)[4] <- "rsID"

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

##################################################################################################################

###Damage filters
###import config file settings
X1000g <- config$V2[config$V1 == "1000G"]
exac <- config$V2[config$V1 == "EXAC"]
CADD <- config$V2[config$V1 == "CADD"]
HET <- config$V2[config$V1 == "HET"]
ADEPTH <- config$V2[config$V1 == "ADEPTH"]
QUAL <- config$V2[config$V1 == "QUAL"]

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
varcount <- paste("##Variant Filter Script ## R-script Log - Variants matching functional consequence:",nrow(exonic.ft))
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
varcount <- paste("##Variant Filter Script ## R-script Log - Variants matching functional consequence:",nrow(func.ft))
write(varcount, file = "R_log.txt", append = TRUE)

##################################################################################################################

###filter by qual
qual.ft <- as.data.frame(vv.ft$ID[vv.ft$QUAL > QUAL])
names(qual.ft)[1] <- "ID"
vv.ft <- vv.ft[vv.ft$ID %in% qual.ft$ID,]
###log number of variants - func
varcount <- paste("##Variant Filter Script ## R-script Log - Variants matching QUAL > 100:",nrow(qual.ft))
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
vv.ft$HET_rate <- hetpct
vv.ft$HOM_rate <- hompct
vv.ft$MISS_rate <- misspct
rm(hetpct, hompct, misspct, calc, HETp, nonHOM, refHOM, miss, gtm)

##################################################################################################################

###aggregate mutation types and af counts - ONLY PERFORMED ON HET CALLS!
###all HET call gene sums - after functional exonic only - QUAL filtered
aggAF_all_HET <- aggregate(vv.ft$HET_val, by=list(GENE=vv.ft$GENE), drop = FALSE, FUN=sum)
#
##Counts for different types of variants - splicing class by either CONSEQUENCE OR TYPE
aggAF_splc_HET <- aggregate(vv.ft$HET_val[vv.ft$CONSEQUENCE == "-9" | vv.ft$TYPE == "splicing" | vv.ft$TYPE == "exonic;splicing" ], 
     by=list(GENE=vv.ft$GENE[vv.ft$CONSEQUENCE == "-9" | vv.ft$TYPE == "splicing" | vv.ft$TYPE == "exonic;splicing" ]), drop = FALSE, FUN=sum)

###Counts for different types of variants - nonsynonymous (excluding splicing) class by either CONSEQUENCE OR TYPE
aggAF_nsyn_HET <- aggregate(vv.ft$HET_val[vv.ft$CONSEQUENCE == "nonsynonymous SNV" & vv.ft$TYPE == "exonic" ], 
                        by=list(GENE=vv.ft$GENE[vv.ft$CONSEQUENCE == "nonsynonymous SNV" & vv.ft$TYPE == "exonic"]), drop = FALSE, FUN=sum)

###Counts for different types of variants - Truncating (including stop loss) class by either CONSEQUENCE OR TYPE
aggAF_trunc_HET <- aggregate(vv.ft$HET_val[vv.ft$CONSEQUENCE == "stopgain" | vv.ft$CONSEQUENCE == "stoploss" 
                                                                           | vv.ft$CONSEQUENCE == "frameshift deletion" 
                                                                           | vv.ft$CONSEQUENCE == "frameshift insertion"
                                                                           | vv.ft$CONSEQUENCE == "nonframeshift insertion"
                                                                           | vv.ft$CONSEQUENCE == "nonframeshift deletion"],
                   by=list(GENE=vv.ft$GENE[vv.ft$CONSEQUENCE == "stopgain" | vv.ft$CONSEQUENCE == "stoploss" 
                                                                           | vv.ft$CONSEQUENCE == "frameshift deletion" 
                                                                           | vv.ft$CONSEQUENCE == "frameshift insertion"
                                                                           | vv.ft$CONSEQUENCE == "nonframeshift insertion"
                                                                           | vv.ft$CONSEQUENCE == "nonframeshift deletion"]),
                                                                           drop = FALSE, FUN=sum)
###Gene count merging to single data.frame - replace na with 0
aggAF_HET <- merge(aggAF_all_HET, aggAF_nsyn_HET, by = "GENE", all = TRUE)
names(aggAF_HET)[2] <- "C_all"
names(aggAF_HET)[3] <- "C_nsyn"
aggAF_HET <- merge(aggAF_HET, aggAF_splc_HET, by = "GENE", all = TRUE)
aggAF_HET <- merge(aggAF_HET, aggAF_trunc_HET, by = "GENE", all = TRUE)
names(aggAF_HET)[4] <- "C_splc"
names(aggAF_HET)[5] <- "C_trunc"
aggAF_HET[is.na(aggAF_HET)] <- 0

###merge the gene counts into the original vv.ft file - adding 4 columns - write table as output - rm variables
###totals of subsets may not all sum to all due to certain combinations (e.g. exonic;splicing stoploss)
vv.ft <- merge(vv.ft, aggAF_HET, by = 'GENE', sort = FALSE)[,union(names(vv.ft), names(aggAF_HET))]
write.table(aggAF_HET,file = "variant_filtering_GeneVarCounts.tsv", quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
rm(aggAF_all_HET, aggAF_nsyn_HET, aggAF_splc_HET, aggAF_trunc_HET, aggAF_HET)

##################################################################################################################

###filter variant ids that are below values for both 1000g & exac_all
rarity.ft <- as.data.frame(anno$ID[anno$X1000g2015aug_all < X1000g & anno$ExAC_ALL < exac])
names(rarity.ft)[1] <- "ID"
###filter variant ids that are above cadd
cadd.ft <- as.data.frame(anno$ID[anno$CADD_phred > CADD | anno$CADD_phred < 0])
names(cadd.ft)[1] <- "ID"

###filter vv by rarity and cadd score
raritycadd.ft <- merge(rarity.ft, cadd.ft,sort = FALSE)
vv.ft <- vv.ft[vv.ft$ID %in% raritycadd.ft$ID,]

###log number of variants - rarity
varcount <- paste("##Variant Filter Script ## R-script Log - Variants matching rarity & CADD threshold:",nrow(vv.ft))
write(varcount, file = "R_log.txt", append = TRUE)
###tidy variables
rm(rarity.ft,raritycadd.ft,func.ft,qual.ft)

##################################################################################################################

###filter by het/hom ratio and het rate (het > hom & no het rate in cohort above 15%)
hethom.ft <- as.data.frame(vv.ft$ID[vv.ft$HET_rate < HET & vv.ft$HET_rate >= vv.ft$HOM_rate])
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
###loop over each sample column
for(i in 2:ncol(af)){
  ###split ad values on ,  to generate a maj allelic read depth and min allelic read depth
  ###replace empty cells in second col with "." to match first
  alfe <- as.data.frame(str_split(af[,i],",",simplify = TRUE),stringsAsFactors = FALSE)
  alfe$V2[alfe$V2 %in% ""] <- "."
  d <- data.frame(x=rep(0,nrow(af)))
  ###loop over new maj/min read depth dataframe & convert to float or replace with na
  for(r in 1:nrow(alfe)){
    if(alfe[r,2] == "."){d[r,1] <- NA}
    else{
      afmaj <- as.numeric(alfe[r,1])
      afmin <- as.numeric(alfe[r,2])
      aftotal <- afmaj + afmin
      afperct <- (afmin/aftotal)
      d[r,1] <- afperct}
  }
  ###replace column in af with new float/na values
  af[,i] <- as.data.frame(d)
  rm(d)
}
###tidy variables and tables
rm(alfe,afmaj,afmin,afperct,aftotal,i,r)
 
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
af.ft <- subset(af.ft, (af.ft[,1] != "NaN"))

vv.ft <- vv.ft[vv.ft$ID %in% af.ft$ID,]
rm(af.ft,i)

###Variants filtered on no variants af above than 0.3
varcount <- paste("##Variant Filter Script ## R-script Log - Variants with at least single AF > threshold:",nrow(vv.ft))
write(varcount, file = "R_log.txt", append = TRUE)

##################################################################################################################

###log number of variants - Final
clock <- as.character(Sys.time())
write(clock, file = "R_log.txt", append = TRUE)
varcount <- paste("##Variant Filter Script ## R-script Log - Final variants passing filters:",nrow(vv.ft))
write(varcount, file = "R_log.txt", append = TRUE)

###Add genotype information for remaining variants
gt.ft <- gt[gt$ID %in% vv.ft$ID,]
vvgt <- merge(vv.ft,gt.ft, sort = FALSE)
###rename ID col - issues with opening files in excel with "ID" as the first value
names(vvgt)[1] <- "Id"
names(af)[1] <- "Id"

##################################################################################################################

###write filtered table out
write.table(vvgt,file = "variant_filtering_results.tsv",sep = "\t",row.names = FALSE, quote = FALSE)
write.table(af,file = "variant_filtering_results_AD.tsv",sep = "\t",row.names = FALSE, quote = FALSE)

