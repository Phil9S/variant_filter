rm(list = ls())
args = commandArgs(trailingOnly=TRUE)
###set working directory and import arguments and libraries
#setwd(args[1])
require("stringr")
clock <- as.character(Sys.time())


###default variables & log start
write(clock, file = "R_log.txt", append = FALSE)
write("##Variant Filter Script ## R-script Log - Log Begin (Version 2)", file = "R_log.txt", append = TRUE)


###Import data tables from bash script
ad <- read.table("allelicdepth.table", header = TRUE,stringsAsFactors = FALSE)
anno <- read.table("annovar.table", header = TRUE, sep = "\t",stringsAsFactors = FALSE, quote="")
gq <- read.table("genoqual.table", header = TRUE,stringsAsFactors = FALSE)
gt <- read.table("genotype.table", header = TRUE,stringsAsFactors = FALSE)
dp <- read.table("sitedepth.table", header = TRUE,stringsAsFactors = FALSE)
config <- read.table("variant_filtering.config",stringsAsFactors = FALSE)
vv <- read.table("variant.table",comment.char = "",header = TRUE,stringsAsFactors = FALSE)
###remove columns that aren't needed
ad <- ad[,-(2:5)]
gq <- gq[,-(2:5)]
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
vv$AA <- anno$AAChange.refGene
vv$CONSEQUENCE <- anno$ExonicFunc.refGene
vv$X1000G <- anno$X1000g2015aug_all
vv$EXAC <- anno$ExAC_ALL
vv$CADD <- anno$CADD_phred



###Rarity & damage filters
###import config file settings
X1000g <- config$V2[config$V1 == "1000G"]
exac <- config$V2[config$V1 == "EXAC"]
CADD <- config$V2[config$V1 == "CADD"]
###filter variant ids that are below values for both 1000g & exac_all
rarity.ft <- as.data.frame(anno$ID[anno$X1000g2015aug_all < X1000g & anno$ExAC_ALL < exac])
names(rarity.ft)[1] <- "ID"
###filter variant ids that are above cadd
cadd.ft <- as.data.frame(anno$ID[anno$CADD_phred > CADD | anno$CADD_phred < 0])
names(cadd.ft)[1] <- "ID"

###log number of variants - rarity
varcount <- paste("##Variant Filter Script ## R-script Log - Variants matching rarity threshold:",nrow(rarity.ft))
write(varcount, file = "R_log.txt", append = TRUE)
###log number of variants - cadd
varcount <- paste("##Variant Filter Script ## R-script Log - Variants matching CADD threshold:",nrow(cadd.ft))
write(varcount, file = "R_log.txt", append = TRUE)

###filter vv by rarity and cadd score
raritycadd.ft <- merge(rarity.ft, cadd.ft,sort = FALSE)
vv.ft <- vv[vv$ID %in% raritycadd.ft$ID,]


###filtering on functional consequnce - ExonicFunction.refGene
func.ft <- as.data.frame(anno$ID[grepl("nonsynonymous SNV",anno$ExonicFunc.refGene,fixed = TRUE)
                                 | grepl("stopgain",anno$ExonicFunc.refGene,fixed = TRUE)
                                 | grepl("stoploss",anno$ExonicFunc.refGene,fixed = TRUE)
                                 | grepl("frameshift insertion",anno$ExonicFunc.refGene,fixed = TRUE)
                                 | grepl("frameshift deletion",anno$ExonicFunc.refGene,fixed = TRUE)
                                 | grepl("nonframeshift insertion",anno$ExonicFunc.refGene,fixed = TRUE)
                                 | grepl("nonframeshift deletion",anno$ExonicFunc.refGene,fixed = TRUE)])
###rename col to match other ID cols
names(func.ft)[1] <- "ID"
###filter by functional consequence
vv.ft <- vv.ft[vv.ft$ID %in% func.ft$ID,]

###log number of variants - func
varcount <- paste("##Variant Filter Script ## R-script Log - Variants matching functional consequence:",nrow(func.ft))
write(varcount, file = "R_log.txt", append = TRUE)

###filter by qual
qual.ft <- as.data.frame(vv.ft$ID[vv.ft$QUAL > 100])
names(qual.ft)[1] <- "ID"
vv.ft <- vv.ft[vv.ft$ID %in% qual.ft$ID,]
###log number of variants - func
varcount <- paste("##Variant Filter Script ## R-script Log - Variants matching QUAL > 100:",nrow(qual.ft))
write(varcount, file = "R_log.txt", append = TRUE)


###Add het calls & cummulative genotyping info
gthet <- data.frame(x=rep(0,nrow(vv.ft)))
gthom <- data.frame(x=rep(0,nrow(vv.ft)))
colnames(gthet) <- c("Het_percentage")
colnames(gthom) <- c("Hom_percentage")
gt <- gt[gt$ID %in% vv.ft$ID,]
######################
gtm <- as.matrix(gt[,-1])
gtm[gtm == "0/0"] <- 0
gtm[gtm == "0/1"] <- 1
gtm[gtm == "1/1"] <- 2
#gtm[gtm == "./."] <- NA
gt <- as.data.frame(gtm)
#######################
refHOM2 <- apply(gtm,1, function(x) sum(x == 0))
HET2 <- apply(gtm,1, function(x) sum(x == 1))
nonHOM2 <- apply(gtm, 1, function(x) sum(x == 2))
calc <- cbind(refHOM2,HET2,nonHOM2)

for(r in 2:nrow(gt)){
  refHOM <- as.numeric(sum(gt[r,2:ncol(gt)] == 0))
  HET <- as.numeric(sum(gt[r,2:ncol(gt)] == 1))
  nonHOM <- as.numeric(sum(gt[r,2:ncol(gt)] == 2))
  hetpct <- (HET / (HET + nonHOM + refHOM))*100
  hompct <- (nonHOM / (HET + nonHOM + refHOM))*100
  gthet[r,1] <- hetpct
  gthom[r,1] <- hompct
}
vv.ft$HET_rate <- gthet$Het_percentage
vv.ft$HOM_rate <- gthom$Hom_percentage




















###tidy variables and tables
rm(HET,hetpct,nonHOM,refHOM,r,hompct)


###filter by het/hom ratio and het rate (het > hom & no het rate in cohort above 15%)
hethom.ft <- as.data.frame(vv.ft$ID[vv.ft$HET_rate < 15 & vv.ft$HET_rate >= vv.ft$HOM_rate])
names(hethom.ft)[1] <- "ID" 
#log number of variants - Het/Hom
varcount <- paste("##Variant Filter Script ## R-script Log - Variants matching Het/Hom thresholds:",nrow(hethom.ft))
write(varcount, file = "R_log.txt", append = TRUE)


###Extract variants based on filtered list
vv.ft <- vv.ft[vv.ft$ID %in% hethom.ft$ID,]
rm(gthet,gthom,func.ft,cadd.ft,hethom.ft)



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
    if(max(af[i,2:ncol(af)], na.rm = TRUE) > 0.3){
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


###log number of variants - Final
clock <- as.character(Sys.time())
write(clock, file = "R_log.txt", append = TRUE)
varcount <- paste("##Variant Filter Script ## R-script Log - Final variants passing filters:",nrow(vv.ft))
write(varcount, file = "R_log.txt", append = TRUE)


###Add genotype information for remaining variants
gt.ft <- gt[gt$ID %in% vv.ft$ID,]
vvgt <- merge(vv.ft,gt.ft, sort = FALSE)
###rename ID col
names(vvgt)[1] <- "Id"

###write filtered table out
write.table(vvgt,file = "variant_filtering_results.tsv",sep = "\t",row.names = FALSE, quote = FALSE)
write.table(af,file = "variant_filtering_results_AF.table",sep = "\t",row.names = FALSE, quote = FALSE)

