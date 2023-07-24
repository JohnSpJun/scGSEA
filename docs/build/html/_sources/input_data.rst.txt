.. scGSEA documentation master file, created by
   sphinx-quickstart on Thu Jul 20 09:37:25 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Input Data
==================================

The scGSEA function requires the following data for input.

Input Files
----------------------------------
Input file\ :sup:`*`
    * File containing unnormalized gene expression data in raw read counts or estimated RNA abundance. Supported input files formats include Seurat RDS, H5seurat, H5ad formats as well as 10x Market Exchange (MEX) and HDF5 (h5) formats. 
    * For a Seurat object, the ``$RNA@counts`` slot will be used. For an AnnData object, the ``raw.X`` slot will be used.
Gene set database file\ :sup:`*`
    Gene sets database file(s) from `GSEA website <https://www.gsea-msigdb.org/gsea/msigdb>`_. Upload your own gene set data using the `GMT file format <https://www.genepattern.org/file-formats-guide#GMT>`_ if the gene set(s) of interest is not listed as a choice from MsigDB.
Chip file
    * Individual gene annotations `chip file <https://www.genepattern.org/file-formats-guide#CHIP>`_ containing Gene ID to Gene Symbol mappings. Both Human and mouse chip files are available for download on `GSEA-MSigDB website <https://www.gsea-msigdb.org/gsea/downloads.jsp>`_. 
    * Leave blank if the gene names are already in gene symbols.
Output file name\ :sup:`*`
    Base name of the output files to be saved.

Cell Grouping Data
----------------------------------
Cell grouping data (e.g. clustering data, celltype annotation data, etc.) to be used to aggregate cells to create metacells. **Please use only one of the two parameters.**

Metacell data label
    * Name of the label for cell grouping information. Use with .rds, .h5seurat, and .h5ad file formats. The default label is *seurat_clusters* that is generated from Seurat's `FindClusters <https://satijalab.org/seurat/reference/findclusters>`_ function.
Metacell data file
    * Tab-delimited .txt file containing cell grouping information. Use with zipped 10x MEX file or 10x HDF5 (h5) file formats. 
    * The first column, **Name**, should consist of cell names (IDs) and the second column, **Metacell**, should consist of cell grouping information i.e. cluster, celltype annotation, etc.

.. list-table:: `example_metacell_data_file.txt`
    :align: center
    :widths: auto
    :header-rows: 1
    
    * - Name
      - Metacell
    * - AAACCTGAGCACCGTC-1
      - Smooth muscle
    * - AAACCTGAGGGCTCTC-1
      - Smooth muscle
    * - AAACCTGAGTTTAGGA-1
      - Mast cell
    * - AAACCTGCAGCTGTTA-1
      -	PDGFRA
    * - AAACCTGGTAACGCGA-1
      - Muscle macrophage
    * - AAACCTGGTAATCGTC-1
      - Fibroblast
    * - AAACCTGTCATCGGAT-1
      - T Cell
    * - AAACCTGTCTCTGCTG-1
      - Fibroblast
    * - AAACGGGCAAATTGCC-1
      - Smooth muscle
    * - AAACGGGCAGACAAGC-1
      - Glia
    * - AAACGGGCAGACGCAA-1
      - PDGFRA
    * - AAACGGGGTGTAATGA-1
      - PDGFRA
    * - AAACGGGTCAACACGT-1
      - Fibroblast
    * - AAACGGGTCCGTTGCT-1
      - PDGFRA
    * - AAACGGGTCTGCGACG-1
      - Fibroblast
  

`*` Required