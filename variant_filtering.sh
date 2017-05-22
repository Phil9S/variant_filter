#!/bin/bash

###Static variables
BUILD="hg38"
ANNO="/data/Resources/Software/"


###Command-line variables and script start
INPUT="NULL"
CONFIG="NULL"
PROJECT="Unknown_project"
OUTPUT="NULL"
TEMP="FALSE"

###Help block - provides details on arguments
for arg in "$@"; do
	if [[ "$arg" == "--help" ]]; then
		echo -e "
##Variant Filtering Script ## - Help Documentation

Options:

--project		Specifies a folder to generate results in - by default the folder is named 'Unknown_project_variantfiltering'.
			Folder is generated where the script is located.

--input			Full path to the VCF file you want to filter - REQUIRED

--output		Full path to the desired output folder - REQUIRED

--config		Full path to the provided variant_filtering.config file provided with the script - should contain 2 column tab
			separated file with names and values for each filtering parameter. Current parameters are;
			
			NAME		TYPE		Effect
			MEANDP		Interger	Filters on mean read-depth across all samples provided on provided value
			MAF		Float		Filters sites in found in provided proportion of all samples in VCF
			GQ		Interger	Filters on genotype quality at a given site across all samples
			MISSING		Float		Filters sites in which the provided proportion of genotypes are missing
			1000G		Float		Filters in 1000G data for given rarity (R-Script AND logic with ExAC)
			EXAC		Float		Filters in ExAC data for given rarity (RScript - AND logic with 1000G)
			CADD		Interger	Filters CADD score for each site at the provided value

--TEMP			Only option TRUE - automatically retains the temporary .table files and asociated intermediate files when script
			is run using this option

Example usage			

./variant_filtering.sh --input /home/user/myvariants.vcf --output /home/user/ --project myvariants --config /home/user/variant_filtering.config --TEMP TRUE 
		"			               
		exit
	fi
done


###Default output - No arguments
if [[ $# -eq 0 ]]; then
	echo -e `date +[%D-%R]` "## Variant Filter Script ## - You need to provide at least SOME arguments! Try using --help for documentation and examples!"
	exit
fi


###Command-line argument passing
while [[ $# > 1 ]]
	do
	key="$1"
	case $key in
		--project)
		PROJECT=${2}
		shift
		;;
		--input)
		INPUT=${2}
		shift
		;;
		--output)
                OUTPUT=${2}
                shift
                ;;
		--config)
		CONFIG=${2}
		shift
		;;
		--TEMP)
		TEMP=${2}
		shift
		;;
	esac
	shift
done


###ERROR CATCHING AND ARGUMENT CHECKS
###Input
if  [[ ${INPUT} == "NULL" ]]; then
	echo -e "## Variant Filter Script ## - ERROR - You must provide an input file - Please use --input arguement or see the Help options (--help)"
	exit
fi
if [[ ! -f ${INPUT} ]]; then
	echo -e "## Variant Filter Script ## - ERROR - The input provided does not exist - Please confirm the input file path"	
	exit
fi
###OUTPUT
if  [[ ${OUTPUT} == "NULL" ]]; then
        echo -e "## Variant Filter Script ## - ERROR - You must provide an output directory - Please use --output arguement or see the Help options (--help)"
        exit
fi
if [[ ! -d ${OUTPUT} ]]; then
        echo -e "## Variant Filter Script ## - ERROR - The output provided does not exist - Please confirm the output file path"
        exit
fi
###config
if  [[ ${CONFIG} == "NULL" ]]; then
        echo -e "## Variant Filter Script ## - ERROR - You must provide an config file path - Please use --config arguement or see the Help options (--help)"
        exit
fi
if [[ ! -f ${CONFIG} ]]; then
        echo -e "## Variant Filter Script ## - ERROR - The config file provided does not exist - Please confirm the config file path"
        exit
fi
###TEMP status
if [[ ${TEMP} -ne "FALSE" ]] && [[ ${TEMP} -ne "TRUE" ]]; then
        echo -e "## Variant Filter Script ## - ERROR - Unknown TEMP file string used - Please provide TRUE if you wish to retain temporary files"
        exit
fi


###Log start
echo -e "\n"
echo -e `date +[%D-%R]` "## Variant Filter Script ## - Script started" | tee -a variantfilter.log
echo -e `date +[%D-%R]` "## Variant Filter Script ## - Project name set to ${2}" | tee -a variantfilter.log
echo -e `date +[%D-%R]` "## Variant Filter Script ## - Input VCF is ${2}" | tee -a variantfilter.log
echo -e `date +[%D-%R]` "## Variant Filter Script ## - Output folder ${2}" | tee -a variantfilter.log
echo -e `date +[%D-%R]` "## Variant Filter Script ## - Using config file ${2}" | tee -a variantfilter.log
echo -e `date +[%D-%R]` "## Variant Filter Script ## - Temporary file retention set to ${2}" | tee -a variantfilter.log

###config arugment settings
MEANDP=$(grep MEANDP ${CONFIG} | cut -f2)
echo -e `date +[%D-%R]` "## Variant Filter Script ## - Minimum Mean Read Depth set to ${MEANDP}" | tee -a variantfilter.log
MAF=$(grep MAF ${CONFIG} | cut -f2)
echo -e `date +[%D-%R]` "## Variant Filter Script ## - Maximum cohort Minor Allele Frequency set to ${MAF}" | tee -a variantfilter.log
GQ=$(grep GQ ${CONFIG} | cut -f2)
echo -e `date +[%D-%R]` "## Variant Filter Script ## - Minumum Genotype Quality set to ${GQ}" | tee -a variantfilter.log
MISSING=$(grep MISSING ${CONFIG} | cut -f2)
echo -e `date +[%D-%R]` "## Variant Filter Script ## - Maximum number of missing values at site to ${MISSING}" | tee -a variantfilter.log

###Progress checkpoint for variable check
while true; do
    read -p "## Variant Filter Script ## - Do you wish to use the variables set above? (y/n) " yn
    case $yn in
        [Yy]* ) break;;
        [Nn]* ) exit;;
        * ) echo "Please answer yes (y) or no (n)";;
    esac
done

###Run folder setup and config file availability
if [[ ! -d ${OUTPUT}${PROJECT}_variantfiltering ]]; then
	mkdir ${OUTPUT}${PROJECT}_variantfiltering
else
	echo -e `date +[%D-%R]` "## Variant Filter Script ## - Project folder already exists - Overwritting Content" | tee -a variantfilter.log
fi

###Migrate required ref files, config and logs to working directory
cp ${CONFIG} ${OUTPUT}${PROJECT}_variantfiltering/variant_filtering.config
cp RVIS_Unpublished_ExACv2_March2017.tsv ${OUTPUT}${PROJECT}_variantfiltering/RVIS_Unpublished_ExACv2_March2017.tsv
cp GDI_full_10282015.tsv ${OUTPUT}${PROJECT}_variantfiltering/GDI_full_10282015.tsv
mv variantfilter.log ${OUTPUT}${PROJECT}_variantfiltering/variantfilter.log
cd ${OUTPUT}${PROJECT}_variantfiltering


###filter all sites containing ref/ref for all positions & on provided filters
echo -en `date +[%D-%R]` "## Variant Filter Script ## - Filtering VCF on specified values..." | tee -a variantfilter.log
cp ${INPUT} variant_orig.vcf
vcftools --vcf variant_orig.vcf --min-meanDP ${MEANDP} --max-maf ${MAF} --minGQ ${GQ} --max-missing ${MISSING} --recode --out variant_orig > /dev/null 2>&1
mv variant_orig.recode.vcf variant_filtered.vcf
echo -e "\r"`date +[%D-%R]` "## Variant Filter Script ## - Filtering VCF on specified values...Done" | tee -a variantfilter.log


###removing header command & generating intermediate files with bcftools
echo -en `date +[%D-%R]` "## Variant Filter Script ## - Generating intermediate files..." | tee -a variantfilter.log
vcftools --vcf variant_filtered.vcf --max-indv 0 --recode --out annotate > /dev/null 2>&1
sed -n '/#CHROM/,${p}' annotate.recode.vcf  > variant.table
echo -e "\r"`date +[%D-%R]` "## Variant Filter Script ## - Generating intermediate files...Done" | tee -a variantfilter.log


echo -en `date +[%D-%R]` "## Variant Filter Script ## - Generating INFO field tables..." | tee -a variantfilter.log 
###Use bcftools to extract depth/genotype/INFO_field information
bcftools query --print-header -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' -o genotype.table variant_filtered.vcf
sed -i 's/\[[0-9]\+\]//g' genotype.table
bcftools query --print-header -f '%CHROM\t%POS\t%REF\t%ALT[\t%DP]\n' -o sitedepth.table variant_filtered.vcf
sed -i 's/\[[0-9]\+\]//g' sitedepth.table
bcftools query --print-header -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' -o allelicdepth.table variant_filtered.vcf
sed -i 's/\[[0-9]\+\]//g' allelicdepth.table
bcftools query --print-header -f '%CHROM\t%POS\t%REF\t%ALT[\t%GQ]\n' -o genoqual.table variant_filtered.vcf
sed -i 's/\[[0-9]\+\]//g' genoqual.table
###removing incorrect #CHROM to CHROM for R input
sed -i 's/# CHROM/CHROM/' genotype.table
sed -i 's/# CHROM/CHROM/' sitedepth.table
sed -i 's/# CHROM/CHROM/' allelicdepth.table
sed -i 's/# CHROM/CHROM/' genoqual.table
sed -i 's/#CHROM/CHROM/' variant.table
###Removing bcftools tags
sed -i 's/:GT//g' genotype.table
sed -i 's/:DP//g' sitedepth.table
sed -i 's/:AD//g' allelicdepth.table
sed -i 's/:GQ//g' genoqual.table

echo -e "\r"`date +[%D-%R]` "## Variant Filter Script ## - Generating INFO field tables...Done" | tee -a variantfilter.log


###Generating annovar-annotation file for use as table - inc gMAF and damage predictions
echo -en `date +[%D-%R]` "## Variant Filter Script ## - Generating Annovar annotation file..." | tee -a variantfilter.log
${ANNO}annovar/convert2annovar.pl -format vcf4old annotate.recode.vcf --outfile annovarform > /dev/null 2>&1
${ANNO}annovar/table_annovar.pl annovarform ${ANNO}annovar/humandb/ -buildver ${BUILD} -out annotated -remove -protocol refGene,1000g2015aug_all,exac03,avsnp144,dbnsfp30a,clinvar_20160302,cosmic70,nci60,dbscsnv11 -operation g,f,f,f,f,f,f,f,f -nastring -9 > /dev/null 2>&1
mv annotated.hg38_multianno.txt annovar.table
echo -e "\r"`date +[%D-%R]` "## Variant Filter Script ## - Generating Annovar annotation file...Done" | tee -a variantfilter.log


###Check file length for the correct number of rows - all .tables & annovar files should have the same number of rows = to each vcf record
VAR=$(cat variant.table | wc -l)
GENOT=$(cat genotype.table | wc -l)
GENOQ=$(cat genoqual.table | wc -l)
ALLELIC=$(cat allelicdepth.table| wc -l)
SITE=$(cat sitedepth.table | wc -l)
ANNOTATE=$(cat annovar.table | wc -l)

if [ "$GENOT" -ne "$VAR" ] || [ "$GENOQ" -ne "$VAR" ] || [ "$ALLELIC" -ne "$VAR" ] || [ "$SITE" -ne "$VAR" ] || [ "$ANNOTATE" -ne "$VAR" ]; then
	echo -e `date +[%D-%R]` "## Variant Filter Script ## - Mismatched rows - The data tables had the following row counts:" | tee -a variantfilter.log
	echo -e `date +[%D-%R]` "Genotype.table - ${GENOT}" | tee -a variantfilter.log
        echo -e `date +[%D-%R]` "Genoqual.table - ${GENOQ}" | tee -a variantfilter.log
        echo -e `date +[%D-%R]` "Allelicdepth.table - ${ALLELIC}" | tee -a variantfilter.log
	echo -e `date +[%D-%R]` "Sitedepth.table - ${SITE}" | tee -a variantfilter.log
	echo -e `date +[%D-%R]` "Annovar.table - ${ANNOTATE}" | tee -a variantfilter.log
	echo -e `date +[%D-%R]` "Reference VCF - ${VAR}" | tee -a variantfilter.log
	echo -e `date +[%D-%R]` "## Variant Filter Script ## - Exiting now" | tee -a variantfilter.log
	exit
else
	echo -e `date +[%D-%R]` "## Variant Filter Script ## - All data tables contain matching rows" | tee -a variantfilter.log
fi

###index all the files with header ID followed by var1-var(nrows-1)
echo -en `date +[%D-%R]` "## Variant Filter Script ## - Indexing data tables..." | tee -a variantfilter.log
awk -F'\t' -v OFS='\t' 'NR == 1 {print "ID", $0; next} {print "Var"(NR-1), $0}' variant.table > awk.table
mv awk.table variant.table
awk -F'\t' -v OFS='\t' 'NR == 1 {print "ID", $0; next} {print "Var"(NR-1), $0}' genotype.table > awk.table
mv awk.table genotype.table
awk -F'\t' -v OFS='\t' 'NR == 1 {print "ID", $0; next} {print "Var"(NR-1), $0}' genoqual.table > awk.table
mv awk.table genoqual.table
awk -F'\t' -v OFS='\t' 'NR == 1 {print "ID", $0; next} {print "Var"(NR-1), $0}' allelicdepth.table > awk.table
mv awk.table allelicdepth.table
awk -F'\t' -v OFS='\t' 'NR == 1 {print "ID", $0; next} {print "Var"(NR-1), $0}' sitedepth.table > awk.table
mv awk.table sitedepth.table
awk -F'\t' -v OFS='\t' 'NR == 1 {print "ID", $0; next} {print "Var"(NR-1), $0}' annovar.table > awk.table
mv awk.table annovar.table
echo -e "\r"`date +[%D-%R]` "## Variant Filter Script ## - Indexing data tables...Done" | tee -a variantfilter.log

wd=`pwd`
###File clean up
rm annovarform
rm annotate.recode.vcf


###Run R script to filter variants	
echo -e `date +[%D-%R]` "## Variant Filter Script ## - R Script Started" | tee -a variantfilter.log	
echo -en `date +[%D-%R]` "## Variant Filter Script ## - Completing filtering on consequence, allele frequency & rarity..." | tee -a variantfilter.log
Rscript ../variant_filtering.R ${wd} > /dev/null 2>&1
echo -e "\r"`date +[%D-%R]` "## Variant Filter Script ## - Completing filtering consequence, allele frequency & rarity...Done" | tee -a variantfilter.log

###Clean up temporary files
if [ "$TEMP" == "FALSE" ]; then
	echo -en `date +[%D-%R]` "## Variant Filter Script ## - Removing temporary files..." | tee -a variantfilter.log
	rm allelicdepth.table
	rm annovar.table
	rm genoqual.table
	rm genotype.table
	rm sitedepth.table
	rm variant.table
	rm variant_orig.vcf
	echo -e "\r"`date +[%D-%R]` "## Variant Filter Script ## - Removing temporary files...Done" | tee -a variantfilter.log 
fi
echo -e `date +[%D-%R]` "## Variant Filter Script ## - Variant Filter Script Finsihed" | tee -a variantfilter.log
