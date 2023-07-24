.. scGSEA documentation master file, created by
   sphinx-quickstart on Thu Jul 20 09:37:25 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Introduction
==================================

scGSEA (single-cell Gene Set Enrichment Analysis)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
scGSEA (single-cell Gene Set Enrichment Analysis) is an extension of `ssGSEA <https://cloud.genepattern.org/gp/pages/index.jsf?lsid=urn:lsid:broad.mit.edu:cancer.software.genepattern.module.analysis:00270:10.1.0>`_ \ :sup:`[2]` (single-sample Gene Set Enrichment Analysis) specifically designed for analyzing single-cell data. ssGSEA is a computational method used to assess the activation or repression of an a priori set of genes associated with a particular biological process or pathway in an individual sample (Barbie et al., Nature 2009). Genes are ranked by their absolute expression in the sample and assessed for their over-representation, i.e., enrichment, at the top or bottom of the list.

While ssGSEA was designed for use with bulk gene expression data, scGSEA addresses the challenge posed by the sparsity which may cause some variability in enrichment scoring due to the large number of genes which are not expressed or lowly expressed in single-cell datasets. To overcome these limitations, scGSEA introduces modifications to the normalization and scoring metrics used in ssGSEA to ensure stability in the computed gene set enrichment scores independent of the order that genes with the same mRNA abundance (or number of reads) appear in the ranked list, which gives rise to the enrichment score. Instead of the normalization options offered by ssGSEA, i.e.,  log-transformation or rank normalization,  scGSEA utilizes centered-log ratio normalization across cells. scGSEA employs the Kolmogorov-Smirnov (KS) statistic for the scoring metric. This approach has been found empirically to better reduce variability in the enrichment scores resulting from reordering of “tied” genes compared to other normalizations and scoring metrics. 

References
^^^^^^^^^^^
    1. Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L., Gillette, M. A., et al. (2005). Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. Proceedings of the National Academy of Sciences of the United States of America, 102(43), 15545-15550. http://doi.org/10.1073/pnas.0506580102
    2. Barbie, D. A., Tamayo, P., Boehm, J. S., et al. (2009). Systematic RNA interference reveals that oncogenic KRAS-driven cancers require TBK1. Nature. 2009;462:108-112. http://doi.org/10.1038/nature08460