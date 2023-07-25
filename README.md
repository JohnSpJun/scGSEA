# scGSEA

**Description**: scGSEA is an extension of ssGSEA tailored for single-cell data analysis. It addresses the challenge of sparsity by employing a normalization method and scoring metric chosen to minimize any variability. By utilizing scGSEA, scientists can explore and interpret pathway activity and functional alterations within heterogeneous populations of cells.

**Authors**: John Jun; UCSD - Mesirov Lab, UCSD

**Contact**: [Forum Link](https://groups.google.com/forum/?utm_medium=email&utm_source=footer#!forum/genepattern-help).  

## Parameters

<table>
    <thead>
        <tr>
            <th>Parameter Group</th>
            <th>Name</th>
            <th>Description</th>
            <th>Default Value</th>
        </tr>
    </thead>
    <tbody>
        <tr>
            <td colspan="1" rowspan="4" align="center">Input Files</td>
            <td>input file *</td> 
            <td>File containing raw counts or mRNA abundance estimates</td>
            <td></td>
        </tr>
        <tr>
            <td>gene set database file *</td>
            <td>Gene sets in GMT format</td>
            <td></td>
        </tr>
        <tr>
            <td>chip file</td>
            <td>Chip file used for conversion to gene symbols</td>
            <td></td>
        </tr>
        <tr>
            <td>output file name *</td>
            <td>Basename to use for output file</td>
            <td align = "center"><i>scGSEA_scores</i></td>
        </tr>
        <tr>
            <td colspan="1" rowspan="2" align="center">Cell Grouping Data * <br><i font size="3">Only use one of the two</i></td>
            <td>metacell data label</td> 
            <td>Metadata label for cell grouping (metacell) information; clustering data</td>
            <td align = "center"><i>seurat_clusters</i></td>
        </tr>
        <tr>
            <td>metacell data file</td> 
            <td>Metadata file for cell grouping (metacell) information; clustering data</td>
            <td></td>
        </tr>
        <tr>
            <td align="center"> Multi-threading</td>
            <td>n_thread</td>
            <td>Number of CPUs to utilize for parallel computing</td>
            <td align="center">3</td>
    </tbody>
</table>

\* Required

## Input Files

1. `input file`  
    This is a file containing unnormalized gene expression data in raw read counts or estimated RNA abundance. The scGSEA module supports multiple input file formats including Seurat RDS, H5seurat, H5ad formats as well as 10x Market Exchange (MEX) and HDF5 (h5) formats. For a Seurat object, the $RNA@counts slot will be used. For an AnnData object, the raw.X slot will be used. 
   * If you come across the following message in the `stderr.txt` file, please verify that the input file contains unnormalized raw counts data.
    &nbsp;<pre><code>The raw counts matrix was not composed of integer values. This may represent an issue with the processing pipeline. Please be advised...</code></pre>
   * If you have used `kallisto` or `salmon.alevin` for alignment, please disregard the message about the raw counts data not being in integer format; the aforementioned tools generate estimated RNA abundances, which may consist of non-integer count values.
   * For 10x MEX file format, please compress the folder containing the three files (barcodes.tsv, matrix.mtx, features.tsv) and supply the `.zip` file.
2. `gene set database file`
    * This parameter’s drop-down allows you to select gene sets from the [Molecular Signatures Database (MSigDB)](https://www.gsea-msigdb.org/gsea/msigdb/index.jsp) on the GSEA website. This drop-down provides access to only the most current (2023) version of MSigDB. You can also upload your own gene set file(s) in [GMT](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29) format.
    * If you want to use files from an earlier version of MSigDB you will need to download them from the archived releases on the [GSEA website](https://www.gsea-msigdb.org/gsea/downloads.jsp).
3. `chip file`  
    This parameter’s drop-down allows you to select CHIP files from the [Molecular Signatures Database (MSigDB)](https://www.gsea-msigdb.org/gsea/msigdb/index.jsp) on the GSEA website. This drop-down provides access to only the most current version (2023) of MSigDB. [How do I choose a chip file?](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/CHIP_File_Selection_Help)

4. `output file name`  
    The prefix used for the name of the output GCT and CSV file. The default output prefix is <i>scGSEA_scores</i>. The output CSV and GCT files will contain a gene set x metacell matrix of enrichments scores.

### Cell Grouping Data
5. `metacell data label`  
    The name of the metadata label for cell grouping information within the input Seurat/AnnData object. This label will be used to access the cell grouping information utilized for aggregating cells to create metacells. The default value for this parameter is <i>seurat_clusters</i>, which is the metadata label for the slot that stores cell-to-cluster mapping generated by the Seurat's [FindClusters](https://satijalab.org/seurat/reference/findclusters) method. Otherwise, provide the appropriate metadata label for the slot that stores cell grouping information.
6. `metacell data file`  
    If your input file is `10x HDF5` or `10x MEX` format, a separate cell grouping data file (tab-delimited .txt file) must be supplied here. The first column, "Name", would have cell names and the second column, "Metacell", would have metacell (cell group) names. The grouping information in this file is used to aggregate cells prior to computing scGSEA scores. Therefore, if you have `10X HDF5` or `10x MEX` formatted files and do not have a metacell data file, please perform clustering using a clustering method of your choice.
    
## Output Files
<!-- list and describe any files output by the module -->

1. `<output_file_name>.csv`   
    This is a gene set x metacell matrix consisted of scGSEA scores. 
2. `<output_file_name>.gct`   
    This is a gene set x metacell matrix consisted of scGSEA scores.

For more details, please refer to the [full documentation](https://github.com/genepattern/scGSEA/blob/develop/docs/documentation.md).
