#!/usr/bin/python3

# This script is a wrapper for the scRNAbox command line interface.
# It is used to run the scRNAbox pipeline from the command line.

# Based on a command line that contain all the parameters available in the scRNAbox pipeline,
# this script will generate the step1_*.txt file and the step2_*.txt file that will be used to run the pipeline.

# The script will also generate the command line that will be used to run the pipeline.

import argparse
import subprocess
import sys
import os
import time

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# Utility methods                  #
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

def is_dir(path):
    if os.path.isdir(path):
        return path
    else:
        raise argparse.ArgumentTypeError(f"{path} is not a valid path")

def is_file(path):
    if os.path.isfile(path):
        return path
    else:
        raise argparse.ArgumentTypeError(f"{path} is not a valid file")


parser = argparse.ArgumentParser(description='scRNAbox command line wrapper')

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#                                Pipeline parameters                             #
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

############################
# Configuration parameters #
############################

configuration_options = parser.add_argument_group('Configuration parameters')

# -R  (--R_path)  =  Path to R executable
configuration_options.add_argument('--R_path', type=is_dir, help='Path to R executable')

# job_mode
configuration_options.add_argument('--job_mode', type=str, help='Job mode (local, slurm)')

######################
# General parameters #
######################

general_options = parser.add_argument_group('General parameters')

# -d  (--dir)  = Working directory (where all the outputs will be printed) (give full path)
general_options.add_argument('-d', '--dir', type=str, help='Working directory (where all the outputs will be printed)')

# --steps  =  Specify what steps, e.g., 2 to run step 2. 2-6, run steps 2 through 6
general_options.add_argument('--steps', type=str, help='Specify what steps, e.g., 2 to run step 2. 2-6, run steps 2 through 6')

# --method  = Select your preferred method: HTO and SCRNA for hashtag, and Standard scRNA, respectively.
general_options.add_argument('--method', type=str, help='Select your preferred method: HTO and SCRNA for hashtag, and Standard scRNA, respectively')

# --msd  = You can get the hashtag labels by running the following code (HTO Step 4).
general_options.add_argument('--msd', action='store_true', help='You can get the hashtag labels by running the following code (HTO Step 4)')

# --markergsea  = Identify marker genes for each cluster and run marker gene set enrichment analysis (GSEA) using EnrichR libraries (Step 7).
general_options.add_argument('--markergsea', action='store_true', help='Identify marker genes for each cluster and run marker gene set enrichment analysis (GSEA) using EnrichR libraries (Step 7)')

# --knownmarkers  = Profile the individual or aggregated expression of known marker genes.
general_options.add_argument('--knownmarkers', action='store_true', help='Profile the individual or aggregated expression of known marker genes')

# --referenceannotation  = Generate annotation predictions based on the annotations of a reference Seurat object (Step 7).
general_options.add_argument('--referenceannotation', action='store_true', help='Generate annotation predictions based on the annotations of a reference Seurat object (Step 7)')

# --annotate  = Add clustering annotations to Seurat object metadata (Step 7).
general_options.add_argument('--annotate', action='store_true', help='Add clustering annotations to Seurat object metadata (Step 7)')

# --addmeta  = Add metadata columns to the Seurat object (Step 8).
general_options.add_argument('--addmeta', action='store_true', help='Add metadata columns to the Seurat object (Step 8)')

# --rundge  = Perform differential gene expression contrasts (Step 8).
general_options.add_argument('--rundge', action='store_true', help='Perform differential gene expression contrasts (Step 8)')

# --seulist  = You can directly call the list of Seurat objects to the pipeline.
general_options.add_argument('--seulist', action='store_true', help='You can directly call the list of Seurat objects to the pipeline')

# --rcheck  = You can identify which libraries are not installed.
general_options.add_argument('--rcheck', action='store_true', help='You can identify which libraries are not installed')

#####################
# Parser for step 1 #
#####################

group1_options = parser.add_argument_group('Step 1')

## Do you want to perform automated library prep?
#par_automated_library_prep= "no"
group1_options.add_argument('--par_automated_library_prep', action='store_true', help='If you want to perform automated library preparation')

## Path to the directory containing the FASTQ files for the experiment. This folder should only contain the FASTQ files for the experiment.
#par_fastq_directory= "/path/to/fastqs/directory"
group1_options.add_argument('--par_fastq_directory', type=is_dir, help='Directory containing the FASTQ files for the experiment')

## List the sample names used in the FASTQ nomenclature
#par_sample_names= c("Sample1", "Sample2", "Sample3")
group1_options.add_argument('--par_sample_names_step1', type=str, help='List the sample names used in the FASTQ nomenclature')

## If you want to rename the samples, set par_rename_samples to yes.
#par_rename_samples= "yes"
group1_options.add_argument('--par_rename_samples', action='store_true', help='If you want to rename the samples')

## If you want to renames the samples (i.e. par_rename_samples= "yes"), list the new sample names in the same order as the old labels are listed in  par_sample_names.
#par_new_sample_names= c("NewSample1", "NewSample2", "NewSample3")
group1_options.add_argument('--par_new_sample_names', type=str, help='List the new sample names in the same order as the old labels are listed in par_sample_names')

## If your sequencing is paired-end, set the following to TRUE. Otherwise set it as FALSE.
#par_paired_end_seq=TRUE
group1_options.add_argument('--par_paired_end_seq', action='store_true', help='If your sequencing is paired-end')

# CellRanger counts pipeline parameters.
## Path to reference genome
#par_ref_dir_grch='/path/to/CellRanger/reference/genome'
group1_options.add_argument('--par_ref_dir_grch', type=is_dir, help='Reference genome')

## Minimum number of bases to retain for R1 sequence of gene expression assay. If you want to use this parameter uncomment the line below and define your par_r1_length.
#par_r1_length=20
group1_options.add_argument('--par_r1_length', type=int, help='Minimum number of bases to retain for R1 sequence of gene expression assay')

## Minimum number of bases to retain for R2 sequence of gene expression assay. If you want to use this parameter uncomment the line below and define your par_r2_length.
#par_r2_length=20
group1_options.add_argument('--par_r2_length', type=int, help='Minimum number of bases to retain for R2 sequence of gene expression assay')

## For clusters whose job managers do not support memory requests, it is possible to request memory in the form of cores. This option will scale up the number of threads requested via the MRO_THREADS variable according to how much memory a stage requires when given to the ratio of memory on your nodes.
#par_mempercode=30
group1_options.add_argument('--par_mempercode', type=int, help='For clusters whose job managers do not support memory requests, it is possible to request memory in the form of cores')

## If you want CellRnager to include introns when producing the gene expression matrices set the following parameter to "yes", otherwise keep the default as "no".
#par_include_introns="no"
group1_options.add_argument('--par_include_introns', action='store_true', help='If you want CellRnager to include introns when producing the gene expression matrices')

## If you want to turn off CellRanger's target UMI filtering subpipeline uncomment the parameter below.
#par_no_target_umi_filter="no"
group1_options.add_argument('--par_no_target_umi_filter', action='store_true', help='If you want to turn off CellRanger\'s target UMI filtering subpipeline')

## If you want to specify the number of expected cells, uncomment the parameter below and enter the value. By default, CellRanger's auto-estimate algorithm will be used.
#par_expect_cells=6000
group1_options.add_argument('--par_expect_cells', type=int, help='If you want to specify the number of expected cells')

## If you want to force the CellRanger count pipeline to use a certain number of cells, uncomment the parameter below and enter the number of cells
#par_force_cells=6000
group1_options.add_argument('--par_force_cells', type=int, help='If you want to force the CellRanger count pipeline to use a certain number of cells')

## If you want to skip the bam file generation, uncomment the parameter below.
#par_no_bam="no"
group1_options.add_argument('--par_no_bam', action='store_true', help='If you want to skip the bam file generation')

#####################
# Parser for step 1 #
#        HTO        #
#####################

group1_options_hto = parser.add_argument_group('Step 1 - HTO')

# par_RNA_run_names	NULL	The names of the sequencing runs for the RNA assay
group1_options_hto.add_argument('--par_RNA_run_names', type=str, help='The names of the sequencing runs for the RNA assay')

# par_HTO_run_names	NULL	The names of the sequencing runs for the HTO assay
group1_options_hto.add_argument('--par_HTO_run_names', type=str, help='The names of the sequencing runs for the HTO assay')

# par_seq_run_names	NULL	The user-selected name for the sequencing run. These names will be used to identify cells in the Seurat objects
group1_options_hto.add_argument('--par_seq_run_names', type=str, help='The user-selected name for the sequencing run')

# id	NULL	Barcode ID which will be used to track the feature counts
group1_options_hto.add_argument('--id', type=str, help='Barcode ID which will be used to track the feature counts')

# name	NULL	The user-selected name for the barcode identifier
group1_options_hto.add_argument('--name', type=str, help='The user-selected name for the barcode identifier')

# read	R2	Which RNA sequencing read contains the barcode sequence. This value Will be either R1 or R2.
group1_options_hto.add_argument('--read', type=str, help='Which RNA sequencing read contains the barcode sequence. This value Will be either R1 or R2')

# pattern	NULL	The pattern of the barcode identifiers
group1_options_hto.add_argument('--pattern', type=str, help='The pattern of the barcode identifiers')

# sequence	NULL	The nucleotide sequence associated with the barcode identifier
group1_options_hto.add_argument('--sequence', type=str, help='The nucleotide sequence associated with the barcode identifier')

#####################
# Parser for step 2 #
#####################

group2_options = parser.add_argument_group('Step 2')

# If you want to save an RNA expression matrix and metadata dataframe set the following to "yes"
############################################################################
#step2_par_save_RNA= "yes"
group2_options.add_argument('--par_save_RNA_step2', action='store_true', help='If you want to save an RNA expression matrix')
#step2_par_save_metadata= "yes"
group2_options.add_argument('--par_save_metadata_step2', action='store_true', help='If you want to save a metadata dataframe')

# Ambient RNA removal
############################################################################
## If you want to remove ambient RNA from the expression matrix, keep the default as “yes”. If you do not want values changed to remove ambient RNA, change to “no”.
# par_ambient_RNA= "yes"
group2_options.add_argument('--par_ambient_RNA', action='store_true', help='If you want to remove ambient RNA from the expression matrix')

# Filtering parameters
############################################################################
## Only retain genes that are present in at least a specified number of cells.
# par_min.cells_L= 3
group2_options.add_argument('--par_min_cells_L', type=int, help='Only retain genes that are present in at least a specified number of cells')

# Normalization and scaling parameters for individual Seurat object
## Note: normalization and scaling must be performed prior to cell cycle scoring
############################################################################
## Normalization method
# par_normalization.method= "LogNormalize"
group2_options.add_argument('--par_normalization_method_step2', type=str, help='Normalization method')

## Scale factor
# par_scale.factor= 10000
group2_options.add_argument('--par_scale_factor_step2', type=int, help='Scale factor')

## Method for choosing the top variable features. vst, mean.var.plot (mvp), dispersion (disp).
# par_selection.method= "vst"
group2_options.add_argument('--par_selection_method_step2', type=str, help='Method for choosing the top variable features')

## Number of features to select as top variable features
# par_nfeatures= 2500
group2_options.add_argument('--par_nfeatures_step2', type=int, help='Number of features to select as top variable features')

#####################
# Parser for step 3 #
#####################

group3_options = parser.add_argument_group('Step 3')

# If you want to save an RNA expression matrix and metadata data frame set the following to "yes"
############################################################################
# par_save_RNA= "yes"
group3_options.add_argument('--par_save_RNA_step3', action='store_true', help='If you want to save an RNA expression matrix')
# par_save_metadata= "yes"
group3_options.add_argument('--par_save_metadata_step3', action='store_true', help='If you want to save a metadata dataframe')

# If you already have a processed Seurat RDS object, and did not perform Step 2 of scRNAbox,
# use parameter this to add the path to the directory containing you Seurat object(s).
# Make sure that no other files/objects are present in the directory besides Seurat RDS objects.
# Uncomment the line to activate the parameter
############################################################################
#par_seurat_object= "/path/to/directory/containing/seurat/object"
group3_options.add_argument('--par_seurat_object_step3', type=is_dir, help='Directory containing Seurat object(s)')

# Quality control parameters
# Uncomment the line to activate the parameter and add the desired value. Cells will be filtered out accordingly.
# L = lower bound threshold
# R = upper bound threshold
############################################################################
## Minimum number of unique RNA transcripts
# par_nFeature_RNA_L= 300
group3_options.add_argument('--par_nFeature_RNA_L', type=int, help='Minimum number of unique RNA transcripts')

## Maximum number of unique RNA transcripts
# par_nFeature_RNA_U= 10000
group3_options.add_argument('--par_nFeature_RNA_U', type=int, help='Maximum number of unique RNA transcripts')

## Minimum number of total RNA transcripts
# par_nCount_RNA_L= 300
group3_options.add_argument('--par_nCount_RNA_L', type=int, help='Minimum number of total RNA transcripts')

## Maximum number of total RNA transcripts
# par_nCount_RNA_U= 20000
group3_options.add_argument('--par_nCount_RNA_U', type=int, help='Maximum number of total RNA transcripts')

## Minimum mitochondrial RNA percentage
# par_mitochondria_percent_L= 0
group3_options.add_argument('--par_mitochondria_percent_L', type=int, help='Minimum mitochondrial RNA percentage')

## Maximum mitochondrial RNA percentage
# par_mitochondria_percent_U= 20
group3_options.add_argument('--par_mitochondria_percent_U', type=int, help='Maximum mitochondrial RNA percentage')

## Minimum ribosomal RNA percentage
# par_ribosomal_percent_L= 0
group3_options.add_argument('--par_ribosomal_percent_L', type=int, help='Minimum ribosomal RNA percentage')

## Maximum ribosomal RNA percentage
# par_ribosomal_percent_U= 100
group3_options.add_argument('--par_ribosomal_percent_U', type=int, help='Maximum ribosomal RNA percentage')

# Parameters to filter out genes
############################################################################
## If you want to filter out mitochondrial and ribosomal genes set the following parameters to "yes". If you do not want to remove them keep the default as "no".
# par_remove_mitochondrial_genes= "no"
group3_options.add_argument('--par_remove_mitochondrial_genes', action='store_true', help='If you want to filter out mitochondrial genes')
# par_remove_ribosomal_genes= "no"
group3_options.add_argument('--par_remove_ribosomal_genes', action='store_true', help='If you want to filter out ribosomal genes')

## If you have specific genes that you want to remove, enter a vector of the genes. Uncomment the line to activate the parameter.
#par_remove_genes= c("gene1", "gene2")
group3_options.add_argument('--par_remove_genes', type=str, help='List of genes to remove')

# Regress genes
############################################################################
## If you want to regress cell cycle genes, set the following parameters to "yes". If you do not want to regress them, keep the default as "no". Note: if you are using your own Seurat object (i.e. not from Step 2), you can only regress cell cycle genes if your Seurat object has the cell cycle score computed.
# par_regress_cell_cycle_genes= "no"
group3_options.add_argument('--par_regress_cell_cycle_genes', action='store_true', help='If you want to regress cell cycle genes')

## If you want to regress a custom list of genes, set the following parameters to "yes". If you do not want to regress a custom list, keep the default as "no".
# par_regress_custom_genes= "no"
group3_options.add_argument('--par_regress_custom_genes', action='store_true', help='If you want to regress a custom list of genes')

## Enter the genes that you want to regress in the list below.
# par_regress_genes= c("gene1", "gene2")
group3_options.add_argument('--par_regress_genes', type=str, help='List of genes to regress')

# Parameters for normalization and scaling after quality control
############################################################################
## Normalization method
# par_normalization.method= "LogNormalize"
group3_options.add_argument('--par_normalization_method_step3', type=str, help='Normalization method')

## Scale factor
# par_scale.factor= 10000
group3_options.add_argument('--par_scale_factor_step3', type=int, help='Scale factor')

## Method for choosing the top variable features. vst, mean.var.plot (mvp), dispersion (disp).
# par_selection.method= "vst"
group3_options.add_argument('--par_selection_method_step3', type=str, help='Method for choosing the top variable features')

## Number of features to select as top variable features
# par_nfeatures= 2500
group3_options.add_argument('--par_nfeatures_step3', type=int, help='Number of features to select as top variable features')

## Number of most variable features to be reported in csv file
# par_top= 10
group3_options.add_argument('--par_top', type=int, help='Number of most variable features to be reported in csv file')

## Total Number of PCs to compute and store for RunPCA
# par_npcs_pca= 30
group3_options.add_argument('--par_npcs_pca', type=int, help='Total Number of PCs to compute and store for RunPCA')

#####################
# Parser for step 4 #
#####################

group4_options = parser.add_argument_group('Step 4')

# If you want to save an RNA expression matrix and metadata data frame set the following to "yes"
############################################################################
# par_save_RNA= "yes"
group4_options.add_argument('--par_save_RNA_step4', action='store_true', help='If you want to save an RNA expression matrix')
# par_save_metadata= "yes"
group4_options.add_argument('--par_save_metadata_step4', action='store_true', help='If you want to save a metadata dataframe')

# If you already have a processed Seurat RDS object, and did not perform Step 3 of scRNAbox use this to add the path to the directory containing you Seurat object(s).
# Make sure that no other files/objects are present in the directory besides Seurat RDS objects.
# Uncomment the line to activate the parameter
############################################################################
#par_seurat_object= "/path/to/directory/containing/seurat/object"
group4_options.add_argument('--par_seurat_object_step4', type=is_dir, help='Directory containing Seurat object(s)')

# Parameters for UMAP dimensional reduction
############################################################################
## Number of dimensions to use as input into UMAP
# par_RunUMAP_dims= 25
group4_options.add_argument('--par_RunUMAP_dims_step4', type=int, help='Number of dimensions to use as input into UMAP')

## Number of neighbouring points to use in local approximation of manifold structure
# par_RunUMAP_n.neighbors= 45
group4_options.add_argument('--par_RunUMAP_n_neighbors_step4', type=int, help='Number of neighbouring points to use in local approximation of manifold structure')

# Parameters for doublet detection and removal (optional)
############################################################################
## If you want to remove predicted doublets from downstream analyses set the following to "yes"
## If you want to keep predicted doublets for further analysis set the following to "no"
# par_dropDN= "yes"
group4_options.add_argument('--par_dropDN', action='store_true', help='If you want to remove predicted doublets from downstream analyses')

## Number of principal components to use as input doublet analysis.
## This can be determined by the bend in the by elbow plot produced in Step 3
# par_PCs= 25
group4_options.add_argument('--par_PCs', type=int, help='Number of principal components to use as input doublet analysis')

## The number of artificial doublets to generate. DoubletFinder is largely invariant to this parameter. We suggest keeping 0.25
# par_pN= 0.25
group4_options.add_argument('--par_pN', type=float, help='The number of artificial doublets to generate')

## Logical representing whether SCTransform was used during original Seurat object pre-processing
# par_sct= FALSE
group4_options.add_argument('--par_sct', action='store_true', help='Logical representing whether SCTransform was used during original Seurat object pre-processing')

##rate_nExp: the doublet rate according to the number of cells
# par_rate_nExp=0.076
group4_options.add_argument('--par_rate_nExp', type=float, help='The doublet rate according to the number of cells')

## Expected doublet rate for each sample. First list sample IDs, then list the corresponding expected doublet rate for each sample depending on the number of recovered or loaded cells. Sample names should be the same ones used in the library.csv file used for Step 1.
# par_sample_names= c("Control1","Parkinson1")
group4_options.add_argument('--par_sample_names_step4', type=str, help='List of sample names')
# par_expected_doublet_rate= c(0.05,0.05)
group4_options.add_argument('--par_expected_doublet_rate', type=str, help='List of expected doublet rates')

#####################
# Parser for step 5 #
#####################

group5_options = parser.add_argument_group('Step 5')

# If you want to save an RNA expression matrix and metadata data frame set the following to "yes"
############################################################################
# par_save_RNA= "yes"
group5_options.add_argument('--par_save_RNA_step5', action='store_true', help='If you want to save an RNA expression matrix')
# par_save_metadata= "yes"
group5_options.add_argument('--par_save_metadata_step5', action='store_true', help='If you want to save a metadata dataframe')

# If you already have a processed Seurat RDS object(s), and did not perform Step 4 of scRNAbox use this to add the path to the directory containing you Seurat object(s).
# Make sure that no other files/objects are present in the directory besides Seurat RDS objects.
# Uncomment the line to activate the parameter
############################################################################
#par_seurat_object= "/path/to/directory/containing/seurat/objects"
group5_options.add_argument('--par_seurat_object_step5', type=is_dir, help='Directory containing Seurat object(s)')

# If you only have one Seurat object and want to skip integration set the following to "yes"
############################################################################
# par_one_seurat= "no"
group5_options.add_argument('--par_one_seurat', action='store_true', help='If you only have one Seurat object and want to skip integration')

# If you have multiple Seurat objects, choose whether you want to integrate or merge the objects
############################################################################
## Integrate Seurat objects
# par_integrate_seurat= "yes"
group5_options.add_argument('--par_integrate_seurat', action='store_true', help='If you want to integrate Seurat objects')

## Merge Seurat objects
# par_merge_seurat= "no"
group5_options.add_argument('--par_merge_seurat', action='store_true', help='If you want to merge Seurat objects')

# Parameters for normalization and scaling
# Even if you opt to skip integration, adjust the following parameters
############################################################################
## Assay to perform normalization and scaling (prior to integration). For most use cases this will be RNA
# par_DefaultAssay= "RNA"
group5_options.add_argument('--par_DefaultAssay', type=str, help='Assay to perform normalization and scaling')

## Normalization method
# par_normalization.method= "LogNormalize"
group5_options.add_argument('--par_normalization_method_step5', type=str, help='Normalization method')

## Scale factor
# par_scale.factor= 10000
group5_options.add_argument('--par_scale_factor_step5', type=int, help='Scale factor')

# Parameters for integration
############################################################################
## Method for detecting top variable features. vst, mean.var.plot (mvp), dispersion (disp)
# par_selection.method= "vst"
group5_options.add_argument('--par_selection_method_step5', type=str, help='Method for detecting top variable features')

## Number of features to select as top variable features for integration
# par_nfeatures= 2500
group5_options.add_argument('--par_nfeatures_step5', type=int, help='Number of features to select as top variable features for integration')

## Which dimensions to use from the CCA to specify the neighbour search space
# par_FindIntegrationAnchors_dim= 25
group5_options.add_argument('--par_FindIntegrationAnchors_dim', type=int, help='Which dimensions to use from the CCA to specify the neighbour search space')

# Parameters for linear dimensional reduction
# even if you opt to skip integration, adjust the following parameters
############################################################################
## Total Number of PCs to compute and store for RunPCA
# par_RunPCA_npcs= 30
group5_options.add_argument('--par_RunPCA_npcs', type=int, help='Total Number of PCs to compute and store for RunPCA')

## Which dimensions to use as input features for RunUMAP
# par_RunUMAP_dims= 25
group5_options.add_argument('--par_RunUMAP_dims_step5', type=int, help='Which dimensions to use as input features for RunUMAP')

## The number of neighbouring points used in local approximations of manifold structure.
# par_RunUMAP_n.neighbors= 45
group5_options.add_argument('--par_RunUMAP_n_neighbors_step5', type=int, help='The number of neighbouring points used in local approximations of manifold structure')

## Whether or not to perform JackStraw computation. This computation takes a long time.
# par_compute_jackstraw= "no"
group5_options.add_argument('--par_compute_jackstraw', action='store_true', help='Whether or not to perform JackStraw computation')

#####################
# Parser for step 6 #
#####################

group6_options = parser.add_argument_group('Step 6')

# If you want to save an RNA expression matrix and metadata dataframe set the following to "yes"
############################################################################
# par_save_RNA= "yes"
group6_options.add_argument('--par_save_RNA_step6', action='store_true', help='If you want to save an RNA expression matrix')
# par_save_metadata= "yes"
group6_options.add_argument('--par_save_metadata_step6', action='store_true', help='If you want to save a metadata dataframe')

# If you already have a processed Seurat RDS object, and did not perform Step 5 of scRNAbox
# use this parameter to add the path to the directory containing your Seurat object.
# Uncomment the line to activate the parameter. Note you can only have one Seurat object at this point.
############################################################################
#par_seurat_object= "/path/to/seurat.rds"
group6_options.add_argument('--par_seurat_object_step6', type=is_file, help='Path to Seurat object')

# If you skipped integration in step 5, set the following to "yes".
# If you performed integration, keep the default as no.
# This will keep the default assay as "RNA"
# If you wish to cluster on with "RNA" and not the integrated data you can do so by setting this to "yes".
############################################################################
# par_skip_integration= "no"
group6_options.add_argument('--par_skip_integration', action='store_true', help='If you skipped integration in step 5')

# Clustering parameters
############################################################################
## Number of PC to use as input to find neighbours. Can be informed by the elbow and jackstraw plots produced in Step 5.
# par_FindNeighbors_dims= 25
group6_options.add_argument('--par_FindNeighbors_dims', type=int, help='Number of PC to use as input to find neighbours')

## Number of PCs to use as input for UMAP. Can be informed by the elbow and jackstraw plots produced in Step 5.
# par_RunUMAP_dims= 25
group6_options.add_argument('--par_RunUMAP_dims_step6', type=int, help='Number of PCs to use as input for UMAP')

## Defines k for the k-nearest neighbour algorithm. The number of neighbours to include when constructing the SNN graph.
# par_FindNeighbors_k.param= 45
group6_options.add_argument('--par_FindNeighbors_k_param', type=int, help='Defines k for the k-nearest neighbour algorithm')

## Sets the cutoff for acceptable Jaccard index when computing the neighbourhood overlap for the SNN construction
# par_FindNeighbors_prune.SNN= 1/15
group6_options.add_argument('--par_FindNeighbors_prune_SNN', type=float, help='Sets the cutoff for acceptable Jaccard index when computing the neighbourhood overlap for the SNN construction')

## Value of the clustering resolution parameter. You may provide multiple resolution values. Use a value above 1.0 if you want to obtain a larger number of smaller communities.
# par_FindClusters_resolution= c(0, 0.05, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0)
group6_options.add_argument('--par_FindClusters_resolution', type=str, help='Value of the clustering resolution parameter')

# Adjusted Rand Index for clustering resolutions
############################################################################
## If you want to compute ARI, set the following to "yes". This takes a long time.
# par_compute_ARI= "yes"
group6_options.add_argument('--par_compute_ARI', action='store_true', help='If you want to compute ARI')

## Number of repetitions of generated clusters to calculate ARI.
# par_RI_reps= 25
group6_options.add_argument('--par_RI_reps', type=int, help='Number of repetitions of generated clusters to calculate ARI')

#####################
# Parser for step 7 #
#####################

group7_options = parser.add_argument_group('Step 7')

# If you want to save an RNA expression matrix and metadata data frame set the following to "yes"
############################################################################
# par_save_RNA= "yes"
group7_options.add_argument('--par_save_RNA_step7', action='store_true', help='If you want to save an RNA expression matrix')
# par_save_metadata= "yes"
group7_options.add_argument('--par_save_metadata_step7', action='store_true', help='If you want to save a metadata dataframe')

# If you already have a processed Seurat RDS object and did not perform Step 6 of scRNAbox use this parameter to add the path to the directory containing your Seurat object.
# Uncomment the line to activate the parameter
# Your Seurat object must already have clusters
############################################################################
#par_seurat_object= "/path/to/seurat.rds"
group7_options.add_argument('--par_seurat_object_step7', type=is_dir, help='Path to Seurat object')

# General parameters for cluster annotation
############################################################################
## The cluster resolution that you want to use. If you skipped integration, use par_level_cluster="RNA_snn_res.0.75", for example, if you want to proceed with a clustering resolution of 0.75
## This parameter can also be set to your annotated cell type names if you want to check expression levels or find markers within the annotated groups.
# par_level_cluster= "integrated_snn_res.0.75"
group7_options.add_argument('--par_level_cluster', type=str, help='The cluster resolution that you want to use')

# Tool 1: Marker GSEA
############################################################################
## Identify cluster specific markers
# par_run_find_marker= "yes"
group7_options.add_argument('--par_run_find_marker', action='store_true', help='If you want to identify cluster specific markers')

## Run EnrichR GSEA on cluster-specific markers. This step should follow the identification of cluster-specific markers.  Additionally, this step can only be run if your HPC allows internet access.
## Note that code is provided to run this locally if your HPC cannot access the internet.
# par_run_enrichR= "no"
group7_options.add_argument('--par_run_enrichR', action='store_true', help='If you want to run EnrichR GSEA on cluster-specific markers')

## Number of top markers based on avg_log2FC
## This is the number of markers to include on a heatmap for visualization.
# par_top_sel= 5
group7_options.add_argument('--par_top_sel', type=int, help='Number of top markers based on avg_log2FC')

## Character vector of EnrichR databases to search for enrichment
## This is only needed if you have internet connection and are running GSEA
# par_db= c("Descartes_Cell_Types_and_Tissue_2021","CellMarker_Augmented_2021","Azimuth_Cell_Types_2021")
group7_options.add_argument('--par_db', type=str, help='Character vector of EnrichR databases to search for enrichment')

# Tool 2: Visualize the expression of known marker genes
######################################################################
## Perform module score computation
# par_run_module_score= "yes"
group7_options.add_argument('--par_run_module_score', action='store_true', help='If you want to perform module score computation')

## Visualize markers (dot plot, violin plot, feature plot)
# par_run_visualize_markers= "yes"
group7_options.add_argument('--par_run_visualize_markers', action='store_true', help='If you want to visualize markers')

## Define the path to a csv file containing the genes sets for module score
# par_module_score= "/path/to/gene_sets.csv"
group7_options.add_argument('--par_module_score', type=is_file, help='csv file containing the genes sets for module score')

## List of markers that you want to visualize (dot plot, violin plot, feature plot)
## Be sure to use the official gene names
# par_select_features_list= c("gene1", "gene2", "gene3")
group7_options.add_argument('--par_select_features_list', type=str, help='List of markers that you want to visualize')

## If you want to define multiple lists of markers to visualize, you can do so with a csv file. The header should contain the list names and all features belonging to the same list should be in the same column. Uncomment the below parameter and enter the location of the csv file. This can be the same csv file used for module score.
#par_select_features_csv= "/path/to/visualize_features.csv"
group7_options.add_argument('--par_select_features_csv', type=is_file, help='csv file containing the markers that you want to visualize')

# Tool 3: Reference-based annotation parameters
######################################################################
## Seurat RDS object to use as the reference
# par_reference= "/path/to/reference_seurat_object.rds"
group7_options.add_argument('--par_reference', type=is_file, help='Seurat RDS object to use as the reference')

## Define an arbitrary name for the reference object. This will be used to name the metadata slot.
# par_reference_name= "reference"
group7_options.add_argument('--par_reference_name', type=str, help='Define an arbitrary name for the reference object')

## Name of a metadata column in the reference Seurat object that contains cell type annotations
# par_level_celltype= "Cell_Type"
group7_options.add_argument('--par_level_celltype', type=str, help='Name of a metadata column in the reference Seurat object that contains cell type annotations')

## How many dimensions to use to find transfer anchors between query and reference dataset
# par_FindTransferAnchors_dim= 50
group7_options.add_argument('--par_FindTransferAnchors_dim', type=int, help='How many dimensions to use to find transfer anchors between query and reference dataset')

## This will increase your RAM usage so set this number mindfully
# par_futureglobalsmaxSize= 60000 * 1024^2
group7_options.add_argument('--par_futureglobalsmaxSize', type=int, help='This will increase your RAM usage so set this number mindfully')

# Annotate parameters
# Annotations from each iteration will be added to the Step 7 Seurat object
######################################################################
## The clustering resolution to annotate
# par_annotate_resolution= "integrated_snn_res.0.75"
group7_options.add_argument('--par_annotate_resolution', type=str, help='The clustering resolution to annotate')

## the name of the metadata slot under which the cluster labels will be stored.
# par_name_metadata= "Celltypes1"
group7_options.add_argument('--par_name_metadata', type=str, help='The name of the metadata slot under which the cluster labels will be stored')

## A list of cluster labels. Make sure you have as the same number of labels as clusters at the defined clustering resolution. Please do not use "_" when naming cell types.
# par_annotate_labels= c("Annot1", "Annot2", "Annot3")
group7_options.add_argument('--par_annotate_labels', type=str, help='A list of cluster labels')

#####################
# Parser for step 8 #
#####################

group8_options = parser.add_argument_group('Step 8')

# If you want to save an RNA expression matrix and metadata data frame set the following to "yes"
############################################################################
# par_save_RNA= "yes"
group8_options.add_argument('--par_save_RNA_step8', action='store_true', help='If you want to save an RNA expression matrix')
# par_save_metadata= "yes"
group8_options.add_argument('--par_save_metadata_step8', action='store_true', help='If you want to save a metadata dataframe')

# If you already have a processed Seurat RDS object, and did not perform Step 7 of scRNAbox use this parameter to add the path to the Seurat object.
# Uncomment the line to activate the parameter
############################################################################
#par_seurat_object= "/path/to/seurat.rds"
group8_options.add_argument('--par_seurat_object_step8', type=is_file, help='Path to Seurat object')

# Add metadata parameters
############################################################################
## Define the column from the Seurat object that you want to use to add the new metadata
# par_merge_meta= "Sample_ID"
group8_options.add_argument('--par_merge_meta', type=str, help='Column from the Seurat object that you want to use to add the new metadata')

## Enter the path to a csv file describing new metadata that should be added to the Seurat object to facilitate DEG analysis.
## The rows should contain the data to add in the order of the levels of "Sample_ID" or the metadata slot you will use to define your samples. The column names should be the desired name of the metadata slot to add.
# par_metadata= "path/to/metadata.csv"
group8_options.add_argument('--par_metadata', type=is_file, help='Path to a csv file describing new metadata that should be added to the Seurat object')

# Run DGE parameters
# Choose which differential gene expression (DGE) methods you want to use for this submission
# Be sure adjust the appropriate txt design files.
############################################################################
## Perform cell-base DGE with all cells
# par_run_cell_based_all_cells= "yes"
group8_options.add_argument('--par_run_cell_based_all_cells', action='store_true', help='If you want to perform cell-based DGE with all cells')

## Perform cell-based on each cell type group
# par_run_cell_based_celltype_groups= "yes"
group8_options.add_argument('--par_run_cell_based_celltype_groups', action='store_true', help='If you want to perform cell-based DGE on each cell type group')

## Perform Sample-based DGE with all cells (pseudobulk)
# par_run_sample_based_all_cells= "yes"
group8_options.add_argument('--par_run_sample_based_all_cells', action='store_true', help='If you want to perform sample-based DGE with all cells')

## Perform Sample-based DGE on each cell type group (pseudobulk)
# par_run_sample_based_celltype_groups= "yes"
group8_options.add_argument('--par_run_sample_based_celltype_groups', action='store_true', help='If you want to perform sample-based DGE on each cell type group')

# Cell-replicate DGE parameters
############################################################################
## Which statistical method to use when computing DGE using individual cells as replicates
# par_statistical_method= "MAST"
group8_options.add_argument('--par_statistical_method', type=str, help='Which statistical method to use when computing DGE using individual cells as replicates')

# contrast_cell_based_all_cells
group8_options.add_argument('--step8_contrast_cell_based_all_cells', type=is_file, help='To perform cell-based DGE using all cells')

# contrast_cell_based_celltype_groups
group8_options.add_argument('--step8_contrast_cell_based_celltype_groups', type=is_file, help='To perform cell-based DGE using cell type groups')

# contrast_sample_based_all_cells
group8_options.add_argument('--step8_contrast_sample_based_all_cells', type=is_file, help='To perform sample-based DGE using all cells')

# contrast_sample_based_celltype_groups
group8_options.add_argument('--step8_contrast_sample_based_celltype_groups', type=is_file, help='To perform sample-based DGE using all cells')


#################################
# Validate the input parameters #
#################################

args = parser.parse_args()

# check if no arguments are provided
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(0)

if args.dir == None or args.steps == None:
    print("Please provide the directory and the steps to execute")
    sys.exit(1)

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# Generate the step*.txt files based on the user input and the default parameters #
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

# Create working directory
# Check if the directory already exists
if os.path.exists(args.dir):
    print("Directory already exists, use this directory to continue the analysis.")
else:
    os.makedirs(args.dir)
    # Create job_info directory
    os.makedirs(args.dir + "/job_info")
    # Create configs directory
    os.makedirs(args.dir + "/job_info/configs")
    # Create parameters directory
    os.makedirs(args.dir + "/job_info/parameters")
    # Create .tmp directory
    os.makedirs(args.dir + "/job_info/.tmp")
    # Create logs directory
    os.makedirs(args.dir + "/job_info/logs")
    # touch summary_report.txt
    open(args.dir + "/job_info/summary_report.txt", "w").close()

#################################
# Generate scrnabox_config.ini  #
#################################

config_file_path = args.dir + "/job_info/configs/scrnabox_config.ini"
print ("Generating scrnabox_config.ini file: " + config_file_path)
if os.path.exists(config_file_path):
    # Copy file with timestamp with standar mv command
    subprocess.run(["mv", config_file_path, config_file_path + "." + str(int(time.time()))])

config_file = open(config_file_path, "w")

if args.R_path:
    config_file.write("R_LIB_PATH=" + args.R_path + "\n")
else:
    # Check for R_LIB_PATH in the environment variables
    if "R_LIB_PATH" in os.environ:
        config_file.write("R_LIB_PATH=" + os.environ["R_LIB_PATH"] + "\n")
    else:
        print("Please provide the path to the R library directory using the --R_path argument or set the R_LIB_PATH environment variable.")
        sys.exit(1)

if args.method:
    config_file.write("SCRNA_METHOD=" + args.method + "\n")

if args.job_mode:
    config_file.write("JOB_MODE=" + args.job_mode + "\n")
else:
    config_file.write("JOB_MODE=local\n")

# close the file
config_file.close()

#################################
# Generate the step1_*.txt file #
#################################

step1_file_path = args.dir + "/job_info/parameters/step1_par.txt"
step1_file = open(step1_file_path, "w")

if args.par_automated_library_prep:
    step1_file.write("par_automated_library_prep=\"yes\"\n")
else:
    step1_file.write("par_automated_library_prep=\"no\"\n")

if args.par_fastq_directory:
    step1_file.write("par_fastq_directory='" + args.par_fastq_directory + "'\n")

if args.par_sample_names_step1:
    step1_file.write("par_sample_names= c(" + args.par_sample_names_step1 + ")\n")

if args.par_rename_samples:
    step1_file.write("par_rename_samples=\"yes\"\n")
else:
    step1_file.write("par_rename_samples=\"no\"\n")

if args.par_new_sample_names:
    step1_file.write("par_new_sample_names= c(" + args.par_new_sample_names + ")\n")

if args.par_paired_end_seq:
    step1_file.write("par_paired_end_seq=TRUE\n")
else:
    step1_file.write("par_paired_end_seq=FALSE\n")

if args.par_ref_dir_grch:
    step1_file.write("par_ref_dir_grch='" + args.par_ref_dir_grch + "'\n")

if args.par_r1_length:
    step1_file.write("par_r1_length=" + str(args.par_r1_length) + "\n")

if args.par_r2_length:
    step1_file.write("par_r2_length=" + str(args.par_r2_length) + "\n")

if args.par_mempercode:
    step1_file.write("par_mempercode=" + str(args.par_mempercode) + "\n")
else:
    step1_file.write("par_mempercode=30\n")

if args.par_include_introns:
    step1_file.write("par_include_introns=TRUE\n")
else:
    step1_file.write("par_include_introns=FALSE\n")

if args.par_no_target_umi_filter:
    step1_file.write("par_no_target_umi_filter=\"no\"\n")

if args.par_expect_cells:
    step1_file.write("par_expect_cells=" + str(args.par_expect_cells) + "\n")

if args.par_force_cells:
    step1_file.write("par_force_cells=" + str(args.par_force_cells) + "\n")

if args.par_no_bam:
    step1_file.write("par_no_bam=\"no\"\n")

if args.par_RNA_run_names:
    step1_file.write("par_RNA_run_names= " + args.par_RNA_run_names + "\n")

if args.par_HTO_run_names:
    step1_file.write("par_HTO_run_names= " + args.par_HTO_run_names + "\n")

if args.par_seq_run_names:
    step1_file.write("par_seq_run_names= " + args.par_seq_run_names + "\n")

if args.id:
    step1_file.write("id= " + args.id + "\n")

if args.name:
    step1_file.write("name= " + args.name + "\n")

if args.read:
    step1_file.write("read= " + args.read + "\n")

if args.pattern:
    step1_file.write("pattern= " + args.pattern + "\n")

if args.sequence:
    step1_file.write("sequence= " + args.sequence + "\n")

step1_file.close()

#################################
# Generate the step2_*.txt file #
#################################

step2_file_path = args.dir + "/job_info/parameters/step2_par.txt"
step2_file = open(step2_file_path, "w")

if args.par_save_RNA_step2:
    step2_file.write("par_save_RNA=\"yes\"\n")

if args.par_save_metadata_step2:
    step2_file.write("par_save_metadata=\"yes\"\n")

if args.par_ambient_RNA:
    step2_file.write("par_ambient_RNA=\"yes\"\n")

if args.par_min_cells_L:
    step2_file.write("par_min_cells_L=" + str(args.par_min_cells_L) + "\n")

if args.par_normalization_method_step2:
    step2_file.write("par_normalization_method= " + args.par_normalization_method_step2 + "\n")

if args.par_scale_factor_step2:
    step2_file.write("par_scale_factor=" + str(args.par_scale_factor_step2) + "\n")

if args.par_selection_method_step2:
    step2_file.write("par_selection_method= " + args.par_selection_method_step2 + "\n")

if args.par_nfeatures_step2:
    step2_file.write("par_nfeatures=" + str(args.par_nfeatures_step2) + "\n")

step2_file.close()

#################################
# Generate the step_3.txt file #
#################################

step3_file_path = args.dir + "/job_info/parameters/step3_par.txt"
step3_file = open(step3_file_path, "w")

if args.par_save_RNA_step3:
    step3_file.write("spar_save_RNA=\"yes\"\n")

if args.par_save_metadata_step3:
    step3_file.write("par_save_metadata=\"yes\"\n")

if args.par_seurat_object_step3:
    step3_file.write("par_seurat_object= " + args.par_seurat_object_step3 + "\n")

if args.par_nFeature_RNA_L:
    step3_file.write("par_nFeature_RNA_L=" + str(args.par_nFeature_RNA_L) + "\n")

if args.par_nFeature_RNA_U:
    step3_file.write("par_nFeature_RNA_U=" + str(args.par_nFeature_RNA_U) + "\n")

if args.par_nCount_RNA_L:
    step3_file.write("par_nCount_RNA_L=" + str(args.par_nCount_RNA_L) + "\n")

if args.par_nCount_RNA_U:
    step3_file.write("par_nCount_RNA_U=" + str(args.par_nCount_RNA_U) + "\n")

if args.par_mitochondria_percent_L:
    step3_file.write("par_mitochondria_percent_L=" + str(args.par_mitochondria_percent_L) + "\n")

if args.par_mitochondria_percent_U:
    step3_file.write("par_mitochondria_percent_U=" + str(args.par_mitochondria_percent_U) + "\n")

if args.par_ribosomal_percent_L:
    step3_file.write("par_ribosomal_percent_L=" + str(args.par_ribosomal_percent_L) + "\n")

if args.par_ribosomal_percent_U:
    step3_file.write("par_ribosomal_percent_U=" + str(args.par_ribosomal_percent_U) + "\n")

if args.par_remove_mitochondrial_genes:
    step3_file.write("par_remove_mitochondrial_genes=\"yes\"\n")

if args.par_remove_ribosomal_genes:
    step3_file.write("par_remove_ribosomal_genes=\"yes\"\n")

if args.par_remove_genes:
    step3_file.write("par_remove_genes= c(" + args.par_remove_genes + ")\n")

if args.par_regress_cell_cycle_genes:
    step3_file.write("par_regress_cell_cycle_genes=\"yes\"\n")

if args.par_regress_custom_genes:
    step3_file.write("par_regress_custom_genes=\"yes\"\n")

if args.par_regress_genes:
    step3_file.write("par_regress_genes= c(" + args.par_regress_genes + ")\n")

if args.par_normalization_method_step3:
    step3_file.write("par_normalization_method= " + args.par_normalization_method_step3 + "\n")

if args.par_scale_factor_step3:
    step3_file.write("par_scale_factor=" + str(args.par_scale_factor_step3) + "\n")

if args.par_selection_method_step3:
    step3_file.write("par_selection_method= " + args.par_selection_method_step3 + "\n")

if args.par_nfeatures_step3:
    step3_file.write("par_nfeatures=" + str(args.par_nfeatures_step3) + "\n")

if args.par_top:
    step3_file.write("par_top=" + str(args.par_top) + "\n")

if args.par_npcs_pca:
    step3_file.write("par_npcs_pca=" + str(args.par_npcs_pca) + "\n")

step3_file.close()

#################################
# Generate the step_4.txt file #
#################################

step4_file_path = args.dir + "/job_info/parameters/step4_par.txt"
step4_file = open(step4_file_path, "w")

if args.par_save_RNA_step4:
    step4_file.write("par_save_RNA=\"yes\"\n")

if args.par_save_metadata_step4:
    step4_file.write("par_save_metadata=\"yes\"\n")

if args.par_seurat_object_step4:
    step4_file.write("par_seurat_object= " + args.par_seurat_object_step4 + "\n")

if args.par_RunUMAP_dims_step4:
    step4_file.write("par_RunUMAP_dims=" + str(args.par_RunUMAP_dims_step4) + "\n")

if args.par_RunUMAP_n_neighbors_step4:
    step4_file.write("par_RunUMAP_n_neighbors=" + str(args.par_RunUMAP_n_neighbors_step4) + "\n")

if args.par_dropDN:
    step4_file.write("par_dropDN=\"yes\"\n")

if args.par_PCs:
    step4_file.write("par_PCs=" + str(args.par_PCs) + "\n")

if args.par_pN:
    step4_file.write("par_pN=" + str(args.par_pN) + "\n")

if args.par_sct:
    step4_file.write("par_sct=TRUE\n")

if args.par_rate_nExp:
    step4_file.write("par_rate_nExp=" + str(args.par_rate_nExp) + "\n")

if args.par_sample_names_step4:
    step4_file.write("par_sample_names= c(" + args.par_sample_names_step4 + ")\n")

if args.par_expected_doublet_rate:
    step4_file.write("par_expected_doublet_rate= c(" + args.par_expected_doublet_rate + ")\n")

step4_file.close()

#################################
# Generate the step_5.txt file #
#################################

step5_file_path = args.dir + "/job_info/parameters/step5_par.txt"
step5_file = open(step5_file_path, "w")

if args.par_save_RNA_step5:
    step5_file.write("par_save_RNA=\"yes\"\n")

if args.par_save_metadata_step5:
    step5_file.write("par_save_metadata=\"yes\"\n")

if args.par_seurat_object_step5:
    step5_file.write("step5_par_seurat_object= " + args.par_seurat_object_step5 + "\n")

if args.par_one_seurat:
    step5_file.write("par_one_seurat=\"yes\"\n")

if args.par_integrate_seurat:
    step5_file.write("par_integrate_seurat=\"yes\"\n")

if args.par_merge_seurat:
    step5_file.write("par_merge_seurat=\"yes\"\n")

if args.par_DefaultAssay:
    step5_file.write("par_DefaultAssay= " + args.par_DefaultAssay + "\n")

if args.par_normalization_method_step5:
    step5_file.write("par_normalization_method= " + args.par_normalization_method_step5 + "\n")

if args.par_scale_factor_step5:
    step5_file.write("par_scale_factor=" + str(args.par_scale_factor_step5) + "\n")

if args.par_selection_method_step5:
    step5_file.write("par_selection_method= " + args.par_selection_method_step5 + "\n")

if args.par_nfeatures_step5:
    step5_file.write("par_nfeatures=" + str(args.par_nfeatures_step5) + "\n")

if args.par_FindIntegrationAnchors_dim:
    step5_file.write("par_FindIntegrationAnchors_dim=" + str(args.par_FindIntegrationAnchors_dim) + "\n")

if args.par_RunPCA_npcs:
    step5_file.write("par_RunPCA_npcs=" + str(args.par_RunPCA_npcs) + "\n")

if args.par_RunUMAP_dims_step5:
    step5_file.write("par_RunUMAP_dims=" + str(args.par_RunUMAP_dims_step5) + "\n")

if args.par_RunUMAP_n_neighbors_step5:
    step5_file.write("par_RunUMAP_n_neighbors=" + str(args.par_RunUMAP_n_neighbors_step5) + "\n")

if args.par_compute_jackstraw:
    step5_file.write("par_compute_jackstraw=\"yes\"\n")

step5_file.close()

#################################
# Generate the step_6.txt file #
#################################

step6_file_path = args.dir + "/job_info/parameters/step6_par.txt"
step6_file = open(step6_file_path, "w")

if args.par_save_RNA_step6:
    step6_file.write("step6_par_save_RNA=\"yes\"\n")

if args.par_save_metadata_step6:
    step6_file.write("step6_par_save_metadata=\"yes\"\n")

if args.par_seurat_object_step6:
    step6_file.write("step6_par_seurat_object= " + args.par_seurat_object_step6 + "\n")

if args.par_skip_integration:
    step6_file.write("par_skip_integration=\"yes\"\n")

if args.par_FindNeighbors_dims:
    step6_file.write("par_FindNeighbors_dims=" + str(args.par_FindNeighbors_dims) + "\n")

if args.par_RunUMAP_dims_step6:
    step6_file.write("par_RunUMAP_dims=" + str(args.par_RunUMAP_dims_step6) + "\n")

if args.par_FindNeighbors_k_param:
    step6_file.write("par_FindNeighbors_k_param=" + str(args.par_FindNeighbors_k_param) + "\n")

if args.par_FindNeighbors_prune_SNN:
    step6_file.write("par_FindNeighbors_prune_SNN=" + str(args.par_FindNeighbors_prune_SNN) + "\n")

if args.par_FindClusters_resolution:
    step6_file.write("par_FindClusters_resolution= c(" + args.par_FindClusters_resolution + ")\n")

if args.par_compute_ARI:
    step6_file.write("par_compute_ARI=\"yes\"\n")

if args.par_RI_reps:
    step6_file.write("par_RI_reps=" + str(args.par_RI_reps) + "\n")

step6_file.close()

#################################
# Generate the step_7.txt file #
#################################

step7_file_path = args.dir + "/job_info/parameters/step7_par.txt"
step7_file = open(step7_file_path, "w")

if args.par_save_RNA_step7:
    step7_file.write("step7_par_save_RNA=\"yes\"\n")

if args.par_save_metadata_step7:
    step7_file.write("step7_par_save_metadata=\"yes\"\n")

if args.par_seurat_object_step7:
    step7_file.write("step7_par_seurat_object= " + args.par_seurat_object_step7 + "\n")

if args.par_level_cluster:
    step7_file.write("par_level_cluster= " + args.par_level_cluster + "\n")

if args.par_run_find_marker:
    step7_file.write("par_run_find_marker=\"yes\"\n")

if args.par_run_enrichR:
    step7_file.write("par_run_enrichR=\"yes\"\n")

if args.par_top_sel:
    step7_file.write("par_top_sel=" + str(args.par_top_sel) + "\n")

if args.par_db:
    step7_file.write("par_db= c(" + args.par_db + ")\n")

if args.par_run_module_score:
    step7_file.write("par_run_module_score=\"yes\"\n")

if args.par_run_visualize_markers:
    step7_file.write("par_run_visualize_markers=\"yes\"\n")

if args.par_module_score:
    step7_file.write("par_module_score= " + args.par_module_score + "\n")

if args.par_select_features_list:
    step7_file.write("par_select_features_list= c(" + args.par_select_features_list + ")\n")

if args.par_select_features_csv:
    step7_file.write("par_select_features_csv= " + args.par_select_features_csv + "\n")

if args.par_reference:
    step7_file.write("par_reference= " + args.par_reference + "\n")

if args.par_reference_name:
    step7_file.write("par_reference_name= " + args.par_reference_name + "\n")

if args.par_level_celltype:
    step7_file.write("par_level_celltype= " + args.par_level_celltype + "\n")

if args.par_FindTransferAnchors_dim:
    step7_file.write("par_FindTransferAnchors_dim=" + str(args.par_FindTransferAnchors_dim) + "\n")

if args.par_futureglobalsmaxSize:
    step7_file.write("par_futureglobalsmaxSize=" + str(args.par_futureglobalsmaxSize) + "\n")

if args.par_annotate_resolution:
    step7_file.write("par_annotate_resolution= " + args.par_annotate_resolution + "\n")

if args.par_name_metadata:
    step7_file.write("par_name_metadata= " + args.par_name_metadata + "\n")

if args.par_annotate_labels:
    step7_file.write("par_annotate_labels= c(" + args.par_annotate_labels + ")\n")

step7_file.close()

#################################
# Generate the step_8.txt file #
#################################

step8_file_path = args.dir + "/job_info/parameters/step8_par.txt"
step8_file = open(step8_file_path, "w")

if args.par_save_RNA_step8:
    step8_file.write("step8_par_save_RNA=\"yes\"\n")

if args.par_save_metadata_step8:
    step8_file.write("step8_par_save_metadata=\"yes\"\n")

if args.par_seurat_object_step8:
    step8_file.write("step8_par_seurat_object= " + args.par_seurat_object_step8 + "\n")

if args.par_merge_meta:
    step8_file.write("par_merge_meta= " + args.par_merge_meta + "\n")

if args.par_metadata:
    step8_file.write("par_metadata= " + args.par_metadata + "\n")

if args.par_run_cell_based_all_cells:
    step8_file.write("par_run_cell_based_all_cells=\"yes\"\n")

if args.par_run_cell_based_celltype_groups:
    step8_file.write("par_run_cell_based_celltype_groups=\"yes\"\n")

if args.par_run_sample_based_all_cells:
    step8_file.write("par_run_sample_based_all_cells=\"yes\"\n")

if args.par_run_sample_based_celltype_groups:
    step8_file.write("par_run_sample_based_celltype_groups=\"yes\"\n")

if args.par_statistical_method:
    step8_file.write("par_statistical_method= " + args.par_statistical_method + "\n")

step8_file.close()

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# Generate the command line         #
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

# Execute launch_scrnabox.sh script
command = [ "launch_scrnabox.sh", "-d", args.dir, "--steps", args.steps ]

# Add msd if provided
if args.msd:
    command.append("--msd")

# Add markergsea if provided
if args.markergsea:
    command.append("--markergsea")

# Add knownmarkers if provided
if args.knownmarkers:
    command.append("--knownmarkers")

# Add referenceannotation if provided
if args.referenceannotation:
    command.append("--referenceannotation")

# Add annotate if provided
if args.annotate:
    command.append("--annotate")

# Add addmeta if provided
if args.addmeta:
    command.append("--addmeta")

# Add rundge if provided
if args.rundge:
    command.append("--rundge")

# Add seulist if provided
if args.seulist:
    command.append("--seulist")

# Add rcheck if provided
if args.rcheck:
    command.append("--rcheck")

print ("Executing command:\n" + " ".join(command) + "\n")
result = subprocess.run(command)

if result.returncode != 0:
    print("Error while executing launch_scrnabox.sh script")
    sys.exit(1)