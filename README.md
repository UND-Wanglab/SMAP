# SMAP
README.md
Welcome to SMAP (Sample Matching for Large-scale Protomomics)!

SMAP is a pipeline to validate and correct sample identity based on a combination of concordance and specificity scores. SMAP first detects variant peptides from multiplexed isobaric labeling-based quantitative proteomics data using the proteogenomics approach, and then infers allelic information for each sample based on its expression level of the variant peptides.

Highlights

use mutant peptide to get the peptide's chromosome and start position

use the peptide's chromosome and start position to extract the corresponding reference peptides in the file of

extract quantification data from both reference and mutation peptides from the file of quantity results of peptides.

Sex information is not necessary but having both genetics-based and reported sexes will help identify true IDs.

You can visit the Shiny platform https://smap.shinyapps.io/smap/

You can visit the github https://github.com/XWangLab/SMAP
