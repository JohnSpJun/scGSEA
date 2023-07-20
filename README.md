# scGSEA

<!-- A brief text description of the module, usually one sentence in length. -->
**Description**: scGSEA is an extension of ssGSEA designed to improve the assessment of pathway activity in single-cell data by addressing sparsity and ensuring stable enrichment scoring. This module is intended to be used subsequent to [Seurat.Clustering](https://github.com/genepattern/Seurat.Clustering) module and the user can supply the Seurat RDS file.

<!-- This field is for the author or creator of the module. If the algorithm of the module is from a published paper, this is usually the first or corresponding author from the paper. If the module algorithm is unpublished, this is usually the developer of the module itself. This field can simply be a name of a person or group. -->
**Authors**: John Jun; UCSD - Mesirov Lab, UCSD

<!--This field is used for responding to help requests for the module, and should be an email address or a link to a website with contact information or a help forum. -->
**Contact**: [Forum Link](https://groups.google.com/forum/?utm_medium=email&utm_source=footer#!forum/genepattern-help).  

## Summary

scGSEA is an extension of ssGSEA tailored for single-cell data analysis. It addresses the challenges of sparsity and unreliable enrichment scoring by employing specialized normalization methods and scoring metrics. By utilizing scGSEA, scientists can explore and interpret pathway activity and functional alterations within heterogeneous populations of cells, thereby advancing our understanding of complex biological systems.

## Parameters
<!-- short description of the module parameters and their default values, as well as whether they are required -->

| Name | Description <!--short description--> | Default Value |
---------|--------------|----------------
| input_file * |  File to be read in RDS format |
| chip_file  | Chip file used for conversion to gene symbols |
| gene_set_database_file *  | Gene set data in GMT format |
| output_file_name * | The basename to use for output file | scGSEA_scores

\*  required

## Input Files
<!-- longer descriptions of the module input files. Include information about format and/or preprocessing...etc -->

1. `input_file`
    This is the Seurat RDS file from the Seurat.Clustering module.
2. `chip_file`  
    This parameter’s drop-down allows you to select CHIP files from the [Molecular Signatures Database (MSigDB)](https://www.gsea-msigdb.org/gsea/msigdb/index.jsp) on the GSEA website. This drop-down provides access to only the most current version of MSigDB. You can also upload your own gene set file(s) in [CHIP](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#CHIP:_Chip_file_format_.28.2A.chip.29) format.
4. `gene_set_database_file`
    * This parameter’s drop-down allows you to select gene sets from the [Molecular Signatures Database (MSigDB)](https://www.gsea-msigdb.org/gsea/msigdb/index.jsp) on the GSEA website. This drop-down provides access to only the most current version of MSigDB. You can also upload your own gene set file(s) in [GMT](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29) format.
    * If you want to use files from an earlier version of MSigDB you will need to download them from the archived releases on the [website](https://www.gsea-msigdb.org/gsea/downloads.jsp).
5. `output_file_name`  
    The prefix used for the name of the output GCT and CSV file. If unspecified, output prefix will be set to `scGSEA_scores`. The output CSV and GCT files will contain the projection of input dataset onto a space of gene set enrichments scores.
6. `cluster_data_label`  
    The name of the metadata label within the input Seurat object. This label will be used to access the annotations utilized for aggregating cells. The default value for this parameter is `seurat_clusters`, which is the metadata label for cluster annotations generated upon running Seurat.Clustering module. Use the default value when using the RDS file generated from the [Seurat.Clustering](https://github.com/genepattern/Seurat.Clustering) module.
    
## Output Files
<!-- list and describe any files output by the module -->

1. `<output_file_name>.csv`   
    This is a gene set by cell cluster data consisted of scGSEA scores. 
2. `<output_file_name>.gct`   
    This is a gene set by cell cluster data consisted of scGSEA scores. The HeatmapViewer module can accept this file as input for generating heatmap visualizations.
3. `cluster_expression.csv`   
    This is a gene by cell cluster data consisted of normalized gene expression level. 
4. `stdout.txt`  
    This is standard output from the script.

For more details, please refer to the [full documentation](https://github.com/genepattern/scGSEA/blob/develop/docs/documentation.md).
