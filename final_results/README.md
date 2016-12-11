# Folder

This folder contains the source code used to generate results in the project.

# Contents

* __subset_genes__: These scripts were used conjointly with bedtools calls for subsetting lincRNAs into categories based on promoter/enhancer overlap. Bedtools calls are not included, but subsetted gene lists are available in the data folder.
    + enhanc_assoc.R: Used to produce bedfiles containing promoter regions of genes. Used conjointly with bedtools intersect calls to categorize lincRNAs based on enhancer overlap, or promoter overlap.
    + all_combinations.R: Used conjointly with bedtools to produce categories of lincRNAs based on promoter AND enhancer overlap, as shown in the report.
* __TAD_processing__: Scripts used to filter TADs or splitting them into bins.
    + TAD2bins.R: Code used to split TADs into bins of equally sized length.
    + small_TADs.R: Code used to filter out large encompassing small_TADs.
* __HiC_processing__: Scripts used to normalize and perform operations on HiC matrices.
    + script_norm_HiC.R: Code used to normalize Hi-C data
    + int_TAD_TAD.R: Code used to compute mean interactions per TAD
    + twosteps_normhiC2boundaries.R: Code used to compute boundaries from normalized Hi-C matrices.
* __gene_characterization__: Scipts used to compute different features of lincRNAs and protein coding genes and store them into compact tables containing gene categories (promoter/enhancer overlaps).
    + expression.R: Used to generate a compact table of median expression values for all genes.
    + seq_conservation.R: Code used to generate a compact table of averaged phastCons scored through mammalian and primate evolution for all genes.
    + tissue_spec.R: Code used to compute tissue specificity and generate a compact table.
* __GAT__: contains files required to perform enrichment tests
    + GAT_run_template.sh: Custom template for testing enrichment of many different segment in an annotation. Segments, annotations and workspace can be changed according to the desired test. Parameters are all set to what has been used in the report.
    + enrichment_analysis.R: Script used for processing GAT output. Allows to take all output files in a folder into a table, provided the filenames are consistent.
* __report_figures.R__: Code used to generate figures in the report


