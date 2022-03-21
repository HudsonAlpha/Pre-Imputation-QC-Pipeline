#!/bin/bash

working_dir=$PWD/
wraynor_script=/cluster/home/jtaylor/scripts/Pre-Imputation-QC-Pipeline/dependencies/HRC-1000G-check-bim.pl
snp_sift=/cluster/home/jtaylor/scripts/Pre-Imputation-QC-Pipeline/dependencies/SnpSift.jar
r_script=/cluster/home/jtaylor/scripts/Pre-Imputation-QC-Pipeline/dependencies/qc_imputation_r_script.R

### set variables for options
input_prefix=""
reference_flag=38
keep_files=no
combined_file=yes
het_check=no
mend_check=no
hardy=no

### parse arguments for the script
for arg in "$@"; do
	case "$arg" in
  		"-i" | "--input" ) input_prefix="$2"; shift ;;
		"-r" | "--reference" ) reference_flag=19; shift ;;
		"--keep-files" ) keep_files=yes; shift ;;
		"--combine-vcfs" ) combined_file=no; shift ;;
		"--het" ) het_check=yes; shift ;;
		"--mend" ) mend_check=yes; shift ;;
		"--hardy" ) hardy=yes; shift ;;
	esac
done

original_prefix=${input_prefix} #set the original prefix the user inputs

### set two reference files to hg19 or hg38 depending on user's input	
if [[ $reference_flag -eq 38 ]]
then
	freeze8_1=/cluster/home/jtaylor/scripts/Pre-Imputation-QC-Pipeline/dependencies/chrALL.BRAVO_TOPMed_Freeze_8_hg38.tab.gz
	freeze8_2=/cluster/home/jtaylor/scripts/Pre-Imputation-QC-Pipeline/dependencies/PASS.Variantschr.BRAVO_TOPMed_Freeze_8_hg38.tab.gz	

elif [[ $reference_flag -eq 19 ]]
then
	freeze8_1=/cluster/home/jtaylor/scripts/Pre-Imputation-QC-Pipeline/dependencies/chrALL.BRAVO_TOPMed_Freeze_8_hg19.tab.gz
	freeze8_2=/cluster/home/jtaylor/scripts/Pre-Imputation-QC-Pipeline/dependencies/PASS.Variantschr.BRAVO_TOPMed_Freeze_8_hg19.tab.gz

fi

### load modules for SspSift and plink
module load cluster/samtools/1.9
module load cluster/plink/1.90
module load cluster/R/4.1.1

### make a directory for the output
mkdir qc_output

### filters heterozygostic outliers
if [ "$het_check" = "yes" ]
then
	
	plink \
		--bfile ${input_prefix} \
		--het \
		--out ${input_prefix}_het
	
	Rscript --vanilla $r_script ${input_prefix}_het
	
	plink \
		--bfile ${input_prefix} \
		--keep ${input_prefix}_het.valid.sample \
		--make-bed \
		--out ${input_prefix}_het
	
	
	input_prefix=${input_prefix}_het #set the new input prefix
	
	### move .het file to output directory
	mv ${input_prefix}.het ${working_dir}qc_output/
	mv ${input_prefix}.valid.sample ${working_dir}qc_output/
	
	### make a new directory for flip exclusions
	mkdir ${input_prefix}

	### move the files to the new directory
	mv ${input_prefix}.* ${input_prefix}

	### change to then new directory
	cd ${input_prefix}
	
fi

### apply a 5% missingness by person and 1% missingness by geno filter, convert the X chromosome to 23 for Wraynor script, and subset to chr1-23 for the TOPMed database
plink \
	--bfile ${input_prefix} \
	--mind 0.05 \
	--geno 0.01 \
	--make-bed \
	--out ${input_prefix}_mind-0.05_geno-0.01_plink1_chr1-23 \
	--chr 1-23

input_prefix=${input_prefix}_mind-0.05_geno-0.01_plink1_chr1-23 #set the new input prefix

### make a directory for new files
mkdir ${input_prefix}

### move the files to this directory
mv ${input_prefix}.* ${input_prefix}

### change to the new directory
cd ${input_prefix}

### applies hardy weinburg filter with a p-value of 0.001
if [ "$hardy" = "yes" ]
then
	
	plink \
		--bfile ${input_prefix} \
		--hardy \
		--hwe 0.000001 \
		--make-bed \
		--out ${input_prefix}_hardy
	
	mv *.hwe ${working_dir}qc_output/
	
	input_prefix=${input_prefix}_hardy #set the new input prefix
	
	### make a new directory for flip exclusions
	mkdir ${input_prefix}

	### move the files to the new directory
	mv ${input_prefix}.* ${input_prefix}

	### change to then new directory
	cd ${input_prefix}
	
fi

### mendelian inconsistencies filter
if [ "$mend_check" = "yes" ]
then
	plink \
		--bfile ${input_prefix} \
		--set-me-missing \
		--mendel-duos \
		--mendel-multigen \
		--make-bed \
		--out ${input_prefix}_temp_1
	
	plink \
		--bfile ${input_prefix}_temp_1 \
		--mind 0.05 \
		--geno 0.01 \
		--make-bed \
		--out ${input_prefix}_mendel \
		--chr 1-23
	
	input_prefix=${input_prefix}_mendel #set the new input prefix
	
	### make a new directory for flip exclusions
	mkdir ${input_prefix}

	### move the files to the new directory
	mv ${input_prefix}.* ${input_prefix}

	### change to then new directory
	cd ${input_prefix}
	
fi

### Generate .frq file.
plink \
	--freq \
	--bfile ${input_prefix} \
	--out ${input_prefix}

### run the Wraynor script to compare to ALL TOPMed freeze 8 variants to be comprehensive at this stage
perl $wraynor_script \
	-b ${input_prefix}.bim \
	-f ${input_prefix}.frq \
	-r $freeze8_1 \
	-h \
	--verbose

### The excluded genotypes are flipped to rescue more variants.
plink \
	--bfile ${input_prefix} \
	--flip Exclude-${input_prefix}-HRC.txt \
	--make-bed \
	--out ${input_prefix}_FlipExclusions

input_prefix=${input_prefix}_FlipExclusions

### filter this filter to geno 0.01 to only retain very high call rate variants
plink \
	--bfile ${input_prefix} \
	--geno 0.01 \
	--make-bed \
	--out ${input_prefix}_geno-0.01

input_prefix=${input_prefix}_geno-0.01 #set the new input prefix

### make a new directory for flip exclusions
mkdir ${input_prefix}

### move the files to the new directory
mv ${input_prefix}.* ${input_prefix}

### change to then new directory
cd ${input_prefix}

### generate the .frq file for this set
plink \
	--freq \
	--bfile ${input_prefix} \
	--out ${input_prefix}


### compare to PASS filter variants for TOPMed to ensure highest quality to go into imputation
perl $wraynor_script \
	-b ${input_prefix}.bim \
	-f ${input_prefix}.frq \
	-r $freeze8_2 \
	-h \
	--verbose

### ensure the run_plink script is executable
chmod 755 Run-plink.sh

### run the Run-plink.sh script
sh Run-plink.sh

### rename 23 to X for input to the TOPMed server
for i in ${input_prefix}; do sed -e "s/##contig=<ID=23/##contig=<ID=X/" ${i}-updated-chr23.vcf | awk -F $'\t' 'BEGIN {OFS = FS} { gsub("23","X",$1); print $0 }' > ${i}-updated-chrX.vcf && rm ${i}-updated-chr23.vcf; done

### makes a combined vcf based on user input
if [ "$combined_file" = "yes" ]
then
	for i in ${input_prefix}; do java -Xmx64G -jar $snp_sift split -j ${i}-updated-chr1.vcf ${i}-updated-chr2.vcf ${i}-updated-chr3.vcf ${i}-updated-chr4.vcf ${i}-updated-chr5.vcf ${i}-updated-chr6.vcf ${i}-updated-chr7.vcf ${i}-updated-chr8.vcf ${i}-updated-chr9.vcf ${i}-updated-chr10.vcf ${i}-updated-chr11.vcf ${i}-updated-chr12.vcf ${i}-updated-chr13.vcf ${i}-updated-chr14.vcf ${i}-updated-chr15.vcf ${i}-updated-chr16.vcf ${i}-updated-chr17.vcf ${i}-updated-chr18.vcf ${i}-updated-chr19.vcf ${i}-updated-chr20.vcf ${i}-updated-chr21.vcf ${i}-updated-chr22.vcf ${i}-updated-chrX.vcf > ${i}-updated-Combined.vcf; done
fi

### zip up all VCF files
for i in *.vcf; do bgzip "$i"; done

### Make a new folder to contain only the input files for imputation
mkdir ZippedVCFsForImputationInput

mv *-updated-chr*.vcf.gz ZippedVCFsForImputationInput/

if [ "$combined_file" = "yes" ]
then
	mv ${input_prefix}-updated-Combined.vcf.gz ${working_dir}qc_output/	
fi

### move the directory with vcfs to the output directory	
mv ZippedVCFsForImputationInput ${working_dir}qc_output/
cd $working_dir

### clean up all intermediate files unless the user decides to keep all files
if [ "$keep_files" = "no" ]
then
	if [ "$het_check" = "yes" ]
	then
		rm -r ${original_prefix}_het/
	else
		rm -r ${original_prefix}_mind-0.05_plink1_chr1-23/
	fi
fi
