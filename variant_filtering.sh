#!/bin/bash

###Static variables
BUILD="hg38"
ANNO="/home/pss41/resources/annovar/"
GATK="/data/Resources/Software/Javas/GenomeAnalysisTK.jar"
REF="/data/Resources/References/hg38.bwa/hg38.bwa.fa"

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

--temp			Only option TRUE - automatically retains the temporary .table files and asociated intermediate files when script
			is run using this option

--help-more		Print extended help documentation for altering default reference files etc.

Example usage			

./variant_filtering.sh --input /home/user/myvariants.vcf --output /home/user/ --project myvariants --config /home/user/variant_filtering.config --TEMP TRUE 
		"			               
		exit
	fi
done

for arg in "$@"; do
        if [[ "$arg" == "--help-more" ]]; then
                echo -e "
##Variant Filtering Script ## - Help Documentation

Options:

--project               Specifies a folder to generate results in - by default the folder is named 'Unknown_project_variantfiltering'.
                        Folder is generated where the script is located.

--input                 Full path to the VCF file you want to filter - REQUIRED

--output                Full path to the desired output folder - REQUIRED

--config                Full path to the provided variant_filtering.config file provided with the script - should contain 2 column tab
                        separated file with names and values for each filtering parameter. Current parameters are;

                        NAME            TYPE            Effect
                        MEANDP          Interger        Filters on mean read-depth across all samples provided on provided value
                        MAF             Float           Filters sites in found in provided proportion of all samples in VCF
                        GQ              Interger        Filters on genotype quality at a given site across all samples
                        MISSING         Float           Filters sites in which the provided proportion of genotypes are missing
                        1000G           Float           Filters in 1000G data for given rarity (R-Script AND logic with ExAC)
                        EXAC            Float           Filters in ExAC data for given rarity (RScript - AND logic with 1000G)
                        CADD            Interger        Filters CADD score for each site at the provided value

--temp                  Only option TRUE - automatically retains the temporary .table files and asociated intermediate files when script
                        is run using this option

--help-more             Print extended help documentation for altering default reference files etc.

Default parameters:

--anno			The full path of the annovar directory - e.g. /resources/sofware/annovar/

--gatk			The full path of the GenomeAnalysisTK.jar file - e.g /resources/software/GATK/GenomeAnalysisTK.jar

--build			Required by annovar to deploy the correct annotation information - Only human builds "hg19" and "hg38" are 
			currently accepted for annotation

--ref			The full path for a reference .fasta file for the genome build of your data

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
		--temp)
		TEMP=${2}
		shift
		;;
		--build)
                BUILD=${2}
                shift
                ;;
		--anno)
                ANNO=${2}
                shift
                ;;
		--gatk)
                GATK=${2}
                shift
                ;;
		--ref)
                REF=${2}
                shift
                ;;
	esac
	shift
done


###ERROR CATCHING AND ARGUMENT CHECKS
###Input
if  [[ ${INPUT} == "NULL" ]]; then
	echo -e "## Variant Filter Script ## - ERROR - You must provide an input file - Please use --input arguement or see the help options (--help)"
	exit
fi
if [[ ! -f ${INPUT} ]]; then
	echo -e "## Variant Filter Script ## - ERROR - The input provided does not exist - Please confirm the input file path"	
	exit
fi
###OUTPUT
if  [[ ${OUTPUT} == "NULL" ]]; then
        echo -e "## Variant Filter Script ## - ERROR - You must provide an output directory - Please use --output arguement or see the help options (--help)"
        exit
fi
if [[ ! -d ${OUTPUT} ]]; then
        echo -e "## Variant Filter Script ## - ERROR - The output provided does not exist - Please confirm the output file path"
        exit
fi
###config
if  [[ ${CONFIG} == "NULL" ]]; then
        echo -e "## Variant Filter Script ## - ERROR - You must provide an config file path - Please use --config arguement or see the help options (--help)"
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

###annovar directory
if [[ ! -d ${ANNO} ]]; then
        echo -e "## Variant Filter Script ## - ERROR - Default Annovar directory not found - Please provide with --anno"
        exit
fi

###gatk jar
if [[ ! -f ${GATK} ]]; then
        echo -e "## Variant Filter Script ## - ERROR - Default GATK .jar file not found - Please provide with --gatk"
        exit
fi

###reference genome
if [[ ! -f ${GATK} ]]; then
        echo -e "## Variant Filter Script ## - ERROR - Default reference file not found - Please provide with --ref"
        exit
fi

###Build version - for annovar
if [[ ${BUILD} -ne "hg38" ]] && [[ ${BUILD} -ne "hg19" ]]; then
        echo -e "## Variant Filter Script ## - ERROR - Unknown or Non-human build ID provided - Please provide hg19/hg38 builds with --build"
        exit
fi


###Log start
echo -e "\n"
echo -e `date +[%D-%R]` "## Variant Filter Script ## - Script started" | tee -a variantfilter.log
echo -e `date +[%D-%R]` "## Variant Filter Script ## - Project name set to ${PROJECT}" | tee -a variantfilter.log
echo -e `date +[%D-%R]` "## Variant Filter Script ## - Input VCF is ${INPUT}" | tee -a variantfilter.log
echo -e `date +[%D-%R]` "## Variant Filter Script ## - Output folder ${OUTPUT}" | tee -a variantfilter.log
echo -e `date +[%D-%R]` "## Variant Filter Script ## - Using config file ${CONFIG}" | tee -a variantfilter.log
echo -e `date +[%D-%R]` "## Variant Filter Script ## - Temporary file retention set to ${TEMP}" | tee -a variantfilter.log

###config arugment settings
MEANDP=$(grep MEANDP ${CONFIG} | cut -f2)
echo -e `date +[%D-%R]` "## Variant Filter Script ## - Minimum Mean Read Depth set to ${MEANDP}" | tee -a variantfilter.log
MAF=$(grep MAF ${CONFIG} | cut -f2)
echo -e `date +[%D-%R]` "## Variant Filter Script ## - Maximum cohort Minor Allele Frequency set to ${MAF}" | tee -a variantfilter.log
GQ=$(grep GQ ${CONFIG} | cut -f2)
echo -e `date +[%D-%R]` "## Variant Filter Script ## - Minumum Genotype Quality set to ${GQ}" | tee -a variantfilter.log
MISSING=$(grep MISSING ${CONFIG} | cut -f2)
echo -e `date +[%D-%R]` "## Variant Filter Script ## - Maximum number of missing values at site to ${MISSING}" | tee -a variantfilter.log
G1000=$(grep G1000 ${CONFIG} | cut -f2)
echo -e `date +[%D-%R]` "## Variant Filter Script ## - 1000 genomes rarity filter value set to ${G1000}" | tee -a variantfilter.log
EXAC=$(grep EXAC ${CONFIG} | cut -f2)
echo -e `date +[%D-%R]` "## Variant Filter Script ## - ExAC consortium rarity filter value set to ${EXAC}" | tee -a variantfilter.log
CADD=$(grep CADD ${CONFIG} | cut -f2)
echo -e `date +[%D-%R]` "## Variant Filter Script ## - CADD variant damage prediction threshold set to ${CADD}" | tee -a variantfilter.log
HET=$(grep HET ${CONFIG} | cut -f2)
echo -e `date +[%D-%R]` "## Variant Filter Script ## - Maximum number of heterozygous calls/site over all samples set to ${HET}%" | tee -a variantfilter.log
ADEPTH=$(grep ADEPTH ${CONFIG} | cut -f2)
echo -e `date +[%D-%R]` "## Variant Filter Script ## - Miniumum alt allelic depth (%) for each site across all samples set to ${ADEPTH}" | tee -a variantfilter.log
QUAL=$(grep QUAL ${CONFIG} | cut -f2)
echo -e `date +[%D-%R]` "## Variant Filter Script ## - Miniumum site QUAL threshold set to ${QUAL}" | tee -a variantfilter.log
SAMP_NUM=$(grep SAMP_NUM ${CONFIG} | cut -f2)
echo -e `date +[%D-%R]` "## Variant Filter Script ## - Seperate genotype file generated if sample number exceeds ${SAMP_NUM}" | tee -a variantfilter.log

###Progress checkpoint for variable check
while true; do
    read -p "[# User Input #] ## Variant Filter Script ## - Do you wish to use the variables set above? (y/n) " yn
    case $yn in
        [Yy]* ) break;;
        [Nn]* ) exit;;
        * ) echo "[# User Input #] ## Variant Filter Script ## - Please answer yes (y) or no (n)";;
    esac
done

###Run folder setup and config file availability
if [[ ! -d ${OUTPUT}${PROJECT}_variantfiltering ]]; then
	mkdir ${OUTPUT}${PROJECT}_variantfiltering
else
	echo -e `date +[%D-%R]` "## Variant Filter Script ## - Project folder already exists - Overwritting content" | tee -a variantfilter.log
fi

###Migrate required ref files, config and logs to working directory
cp ${CONFIG} ${OUTPUT}${PROJECT}_variantfiltering/variant_filtering.config
cp variant_filtering.R ${OUTPUT}${PROJECT}_variantfiltering/variant_filtering.R
cp RVIS_Unpublished_ExACv2_March2017.tsv ${OUTPUT}${PROJECT}_variantfiltering/RVIS_Unpublished_ExACv2_March2017.tsv
cp GDI_full_10282015.tsv ${OUTPUT}${PROJECT}_variantfiltering/GDI_full_10282015.tsv
mv variantfilter.log ${OUTPUT}${PROJECT}_variantfiltering/variantfilter.log
cd ${OUTPUT}${PROJECT}_variantfiltering


####filter all sites containing ref/ref for all positions & on provided filters
echo -en `date +[%D-%R]` "## Variant Filter Script ## - Filtering VCF on specified values..." | tee -a variantfilter.log
cp ${INPUT} variant_orig.vcf
vcftools --vcf variant_orig.vcf --non-ref-ac-any 1 --min-meanDP ${MEANDP} --max-maf ${MAF} --minGQ ${GQ} --max-missing ${MISSING} --recode --out variant_orig 
mv variant_orig.recode.vcf variant_filtered.vcf
echo -e "\r"`date +[%D-%R]` "## Variant Filter Script ## - Filtering VCF on specified values...Done" | tee -a variantfilter.log

###Spliting of multiallelic sites
echo -en `date +[%D-%R]` "## Variant Filter Script ## - Splitting multiallelic sites..." | tee -a variantfilter.log
java -jar ${GATK} -T LeftAlignAndTrimVariants -R ${REF} --variant variant_filtered.vcf -o variant_filtered.bi.vcf --splitMultiallelics > /dev/null 2>&1
echo -e "\r"`date +[%D-%R]` "## Variant Filter Script ## - Splitting multiallelic sites...Done" | tee -a variantfilter.log

###removing header command & generating intermediate files with bcftools
echo -en `date +[%D-%R]` "## Variant Filter Script ## - Generating intermediate files..." | tee -a variantfilter.log
vcftools --vcf variant_filtered.bi.vcf --max-indv 0 --recode --out annotate > /dev/null 2>&1
sed -n '/#CHROM/,${p}' annotate.recode.vcf  > variant.table
echo -e "\r"`date +[%D-%R]` "## Variant Filter Script ## - Generating intermediate files...Done" | tee -a variantfilter.log


echo -en `date +[%D-%R]` "## Variant Filter Script ## - Generating vcf info field tables..." | tee -a variantfilter.log 
###Use bcftools to extract depth/genotype/INFO_field information
bcftools query --print-header -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' -o genotype.table variant_filtered.bi.vcf
sed -i 's/\[[0-9]\+\]//g' genotype.table
bcftools query --print-header -f '%CHROM\t%POS\t%REF\t%ALT[\t%DP]\n' -o sitedepth.table variant_filtered.bi.vcf
sed -i 's/\[[0-9]\+\]//g' sitedepth.table
bcftools query --print-header -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' -o allelicdepth.table variant_filtered.bi.vcf
sed -i 's/\[[0-9]\+\]//g' allelicdepth.table
bcftools query --print-header -f '%CHROM\t%POS\t%REF\t%ALT[\t%GQ]\n' -o genoqual.table variant_filtered.bi.vcf
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

echo -e "\r"`date +[%D-%R]` "## Variant Filter Script ## - Generating vcf info field tables...Done" | tee -a variantfilter.log


###Generating annovar-annotation file for use as table - inc gMAF and damage predictions
echo -en `date +[%D-%R]` "## Variant Filter Script ## - Generating annovar annotation file..." | tee -a variantfilter.log
${ANNO}convert2annovar.pl -format vcf4old annotate.recode.vcf --outfile annovarform > /dev/null 2>&1
${ANNO}table_annovar.pl annovarform ${ANNO}humandb/ -buildver ${BUILD} -out annotated -remove -protocol refGene,1000g2015aug_all,exac03,avsnp150,dbnsfp35a,clinvar_20180603,cosmic70,nci60,dbscsnv11 -operation g,f,f,f,f,f,f,f,f -nastring -9 > /dev/null 2>&1
mv annotated.${BUILD}_multianno.txt annovar.table
echo -e "\r"`date +[%D-%R]` "## Variant Filter Script ## - Generating annovar annotation file...Done" | tee -a variantfilter.log


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
rm *.idx

###Run R script to filter variants	
echo -e `date +[%D-%R]` "## Variant Filter Script ## - R Script started" | tee -a variantfilter.log	
echo -en `date +[%D-%R]` "## Variant Filter Script ## - Completing filtering on consequence, allele frequency & rarity..." | tee -a variantfilter.log
Rscript variant_filtering.R ${wd} > /dev/null 2>&1
echo -e "\r"`date +[%D-%R]` "## Variant Filter Script ## - Completing filtering consequence, allele frequency & rarity...Done" | tee -a variantfilter.log
echo -e `date +[%D-%R]` "## Variant Filter Script ## - R Script completed" | tee -a variantfilter.log
###Clean up temporary files
if [ "$TEMP" == "FALSE" ]; then
	echo -en `date +[%D-%R]` "## Variant Filter Script ## - Removing temporary files..." | tee -a variantfilter.log
	rm allelicdepth.table
	rm annovar.table
	rm genoqual.table
	rm genotype.table
	rm sitedepth.table
	rm variant.table
	rm variant_filtered.vcf
	rm variant_filtered.bi.vcf
	rm GDI_full_10282015.tsv
	rm RVIS_Unpublished_ExACv2_March2017.tsv
	rm variant_orig.vcf
	rm variant_filtering.R
	rm variant_filtering.config
	echo -e "\r"`date +[%D-%R]` "## Variant Filter Script ## - Removing temporary files...Done" | tee -a variantfilter.log 
fi
echo -e `date +[%D-%R]` "## Variant Filter Script ## - Script finished - Results file: ${OUTPUT}${PROJECT}_variantfiltering/variant_filtering_results.tsv" | tee -a variantfilter.log
