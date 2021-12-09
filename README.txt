SMAP V1.0.0
----------------------
Contents of this file
----------------------

 * 1. Introduction
 * 2. Software Requirements
 * 3. Hardware Requirements
 * 4. Installation
 * 5. Command Line Arguments
 * 6. Demo and Program Testing
 * 7. Output
 * 8. Maintainers

----------------------
1. Introduction
----------------------

SMAP is a pipeline to validate and correct sample identity based on a combination of concordance and specificity scores. SMAP first detects variant peptides from multiplexed isobaric labeling-based quantitative proteomics data using the proteogenomics approach, and then infers allelic information for each sample based on its expression level of the variant peptides.

You can visit the Shiny platform https://smap.shinyapps.io/smap/
More information could be found in https://sites.google.com/view/smapwanglab/home

----------------------
2. Software Requirements
---------------------- 

The program a cross-platform application that can be run on Windows, Linux, and macOS. SMAP has both standalone and Rshiny versions. The standalone version supports all 64-bit operating systems. The program is written by a combination of Perl and R. The minimum required Perl version should be Perl 5.6 or R 3.6.0.

For the standalone version, no perl dependency is required. 

For Rshiny version, the following packages are needed:
R (version 3.6 or above)
Perl (version 5.6 or above)

R libraries:
	shiny
	shinydashboard
	ggplot2
	DT
	vcfR
	shinycssloaders
	tidyr
	dplyr
	ggrepel

----------------------
3. Hardware Requirements
---------------------- 

The program can be run on any computers, ranging from a single personal computer to high performance computing systems  
   
----------------------
4. Installation
---------------------- 

For the standalone version:

	The source code and example data can be download from https://github.com/UND-Wanglab/SMAP.

	You can download it using the clone command:
	git clone https://github.com/UND-Wanglab/SMAP.git

	or directly download using the following link and unzip it:

	https://github.com/UND-Wanglab/SMAP/archive/refs/heads/main.zip

	After downloading the source code, you can put it in any working directory (e.g. /home/User/Software/SMAP) and ready to execute the program.  


For Rshiny version:

you can run it from Command Line:

	Download the SMAP_github.zip file and unzip and run the following command:

		Rscript -e "shiny::runApp('/path/to/SMAP_github', launch.browser = TRUE)"

or run it from R Studio (must also have R Studio installed):

	Download the SMAP_github.zip file and unzip
	Open either server.R or ui.R in R Studio
	Click the Run App button:


----------------------
5. Command Line Arguments
----------------------
you can use help for usage:

perl SMAP.pl -h

Command line:

perl SMAP.pl -vf variant_peptide_table.txt[file] -g genotype_table.vcf[file] -o result.txt[file]

	--variant_peptide,-vf (A file containing quantitative values of variant peptides; required)
	--genotype, -g (A genotype file used sample verification; required )
	--output, -o (An output filename; required)
	--plex, -p (Multiplex number of the isobaric labeling approach)
	--fold_change, -fc (Signal to Noise ratio (optional; default is 3))
	--noise_level, -nl (The upper threshold of a noise level)
	--version, -h (Print version)
	--help, -h (Print help)
	--license, -l (Print license)

----------------------
6. Demo and Program Testing
----------------------

If you download the standalone program under the folder of /home/User/Software/SMAP, you can test the program using the following command:

perl src/SMAP.pl -g data/genotype_table.vcf -vf data/variant_peptide_table.txt -o output.txt -fc 5

The program takes two inputs: 

1.Variant peptide quantification table. The table contains peptide id, gene/protein name, peptide spectrum match(PSM), SNP id and the quantification data for each sample. 
2.Genotype file in VCF format.


----------------------
7. Output
----------------------

a. Result table (e.g. output.txt)
This table lists the final report containing Sample ID, Inferred ID, CScore and DeltaCScore.

	An example of output:
	Sample ID	Inferred ID	Csore	DeltaCscore
	S2015_1341	S2015_1341	3.31	0.48
	S2015_737	S2015_737	2.75	0.47
	S2015_804	S2015_804	1.90	0.24
	S2015_42	S2015_37	3.25	0.44
	S2015_1555	S2015_1555	2.60	0.43
	S2015_244	S2015_244	2.31	0.41
	S2015_735	S2015_735	3.83	0.51
	S2014_2200	S2015_857	2.43	0.41
	S2016_958	S2016_958	2.00	0.37
	S2016_965	S2016_965	1.96	0.31
	Internal_standard	S2015_1339	1.87	0.01

b. Intermediate results files (More information can be found in our manual document)
There are four intermediate files including inferred_genotype.txt, sample_specific_genotype.snp, sample_specific_genotype.vcf, Score.txt.

c. SMAP.log
This table lists the steps for monitoring.

----------------------
8. Maintainers
----------------------

To submit bug reports and feature suggestions, or to track changes, please contact:

Xusheng Wang (xusheng.wang@und.edu)
 
 
