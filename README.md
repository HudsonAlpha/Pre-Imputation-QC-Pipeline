Title / Pre-Imputation QC Pipeline

Author / Jared Taylor

About / This pipeline prepares array data for upload to TopMed imputation server. It takes a Plink dataset as input (.bed, .bim, .fam) and has option for different qc functions for the data. The pipeline will will output the vcfs in format ready for upload to the imputation server. 

Options / "-i" or "--input" <input prefix> - takes the prefix for the input files
	"-r" or "--reference" - use this flag to change the reference files to hg19
	"--keep-files" - keeps all files generated throughout the different steps of the pipeline
	"--combine-vcfs" - makes a combined vcf
	"--het" - filters heterozygostic outliers
	"--mend" - filters mendelian inconsistencies 
	"--hardy" - hardy Weinberg filter with p-value set to 0.000001

Files / qc_imputation_prep.sh
	dependencies/chrALL.BRAVO_TOPMed_Freeze_8_hg19.tab.gz
	dependencies/chrALL.BRAVO_TOPMed_Freeze_8_hg38.tab.gz
	dependencies/HRC-1000G-check-bim.pl
	dependencies/PASS.Variantschr.BRAVO_TOPMed_Freeze_8_hg19.tab.gz
	dependencies/PASS.Variantschr.BRAVO_TOPMed_Freeze_8_hg38.tab.gz
	dependencies/qc_imputation_r_script.R
	dependencies/SnpSift.jar

Modules / plink 1.90, R 4.1.1, htslib 1.9

Notes / To run, make sure to update all hard coded paths to the actual path that you have access to. Also make sure you are in the directory that contains the array data. Use the "-i" or "--input" flags to input the the prefix for the input files and use any other flags you may want to use. The pipeline will use hg38 reference files unless you use the "r" or "--reference" flags. The pipeline will also remove all intermediate files unless you use the "--keep-files" flag. All output files will be in a new directory named "qc_output".
	sample command "./qc_imputation_prep.sh -i input_array_data --het --hardy"

Help / For any questions, email jtaylor@hudsonalpha.org

