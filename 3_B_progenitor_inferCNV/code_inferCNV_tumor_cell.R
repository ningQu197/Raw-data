# https://www.jianshu.com/p/6b44e511f641
# I

BiocManager::install("infercnv")
library("infercnv")

# 
########################################### 
# WASH7P  chr1    14363   29806
# LINC00115       chr1    761586  762902
# NOC2L   chr1    879584  894689
# MIR200A chr1    1103243 1103332

# 

# 
library(infercnv)
library(Seurat)
library(tidyverse)
B_progenitor = readRDS("./B_progenitor_cells.rds")
reference_data = readRDS("../1_preprocessing/mydata_cluster.rds")
reference_data = subset(reference_data, Type=="PBMMC")
reference_data = subset(reference_data, cell_type=="Cytotoxic NK/T cells")  #
data = as.matrix(mydata@assays$RNA@counts)
annotation = subset(mydata@meta.data, select=c("cell_type"))
genes = read.table("./infercnv_genes.txt", header=T, row.names=1, sep="\t")
#################################################################################################
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=data,
	annotations_file=annotation, delim="\t",
	gene_order_file=genes,
	ref_group_names=c("Cytotoxic NK/T cells")
)
#
infercnv_obj = infercnv::run(infercnv_obj,
	cutoff=0.1, 
	out_dir="./inferCNV_output",
	cluster_by_groups=TRUE, 
	analysis_mode="subclusters",
	HMM_type="i3",
	denoise=TRUE,
	HMM_report_by="subcluster",
	HMM=TRUE,
	output_format="pdf",
	no_prelim_plot=T
)
# analysis_mode = c("samples", "subclusters", "cells"), 默认是samples
# 

# default：i6
# (i6 or i3): 
# i6: infercnv (0, 0.5, 1, 1.5, 2, >2) 
# i3: infercnv (del, neutral, amp) 
# 
# State 1: 0x: complete loss
# State 2: 0.5x: loss of one copy
# State 3: 1x: neutral
# State 4: 1.5x: addition of one copy
# State 5: 2x: addition of two copies
# State 6: 3x: essentially a placeholder for >2x copies but modeled as 3x

# 
# State 1: del
# State 2: neutral
# State 3: amp

# HMM_report_by
#####################################################################################################
# 

# (2)infercnv.references.txt : the 'normal' cell matrix data values
# (3)infercnv.observations.txt : the tumor cell matrix data values

# (4)run.final.infercnv_obj
# @expr.data：
# @reference_grouped_cell_indices：
# @observation_grouped_cell_indices：
#########################################################################################
# (5)HMM_CNV_predictions.*.cell_groupings
# 
# 
# cell_group_name                         cell
# malignant_MGH36.malignant_MGH36_s1      MGH36_P3_E06
# malignant_MGH36.malignant_MGH36_s1      MGH36_P10_E07
#########################################################################################
# (6)HMM_CNV_predictions.*.pred_cnv_genes.dat  
# cell_group_name gene_region_name        state   gene    chr     start   end
# malignant_MGH36.malignant_MGH36_s1      chr1-region_2   2       ACOT7   chr1    6281253 6296032
# malignant_MGH36.malignant_MGH36_s1      chr1-region_2   2       NOL9    chr1    6324329 6454451

# HMM_CNV_predictions.HMMi6.hmm_mode-cells.Pnorm_0.5.pred_cnv_regions.dat  
# 
# cell_group_name	cnv_name	state	chr	start	end
# AAAGAACGTTCGGTAT-1	chr1-region_2	2	chr1	29192657	39034636
# AAAGAACGTTCGGTAT-1	chr1-region_4	4	chr1	183023420	214552449
# AAAGAACGTTCGGTAT-1	chr2-region_7	4	chr2	97646062	127887185







