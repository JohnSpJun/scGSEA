#!/usr/bin/env python3

from scGSEA.scgsea_helper import *
import argparse
import pandas as pd
import logging
import os
import scanpy as sc
import sys
import zipfile
import humanfriendly
from timeit import default_timer as timer

Log_Format = "[%(asctime)s] %(message)s"
logging.basicConfig(level=logging.INFO, 
                    stream = sys.stdout, 
                    filemode = "w",
                    format = Log_Format,
                    datefmt='%Y-%m-%d %H:%M:%S')

logger = logging.getLogger()

# scGSEA main function
# Performs Centered-Log Ratio normalization, aggregates cells based on the metacell data, 
# and computes the enrichment scores for each metacell using the KS statistics.
def scgsea(
    input_file,
    gene_set_file,
    chip_file,
    gene_set_database_file,
    output_file_name = "scGSEA_scores",
    cluster_data_label = "seurat_clusters",
    cluster_data_file = None,
    n_job = 3
):
    logger.info('==========================================================')
    logger.info("About to read the object named")

    split_name = input_file.split(".")
    file_extension = split_name[-1]

    # User supplies AnnData object with annotation information
    if check_file_extension(input_file, "h5ad"):
        scanpy_obj = sc.read_h5ad(input_file)
        use_metacell_data_file = False

    # User supplies raw counts data and annotation file
    elif check_file_extensions(input_file, ["h5", "hdf5"]):
        scanpy_obj = sc.read_10x_h5(input_file, gex_only = True)
        use_metacell_data_file = True

    elif check_file_extension(input_file, "zip"):
        with zipfile.ZipFile(input_file, 'r') as zip_ref:
            zip_ref.extractall("./data")

        data_paths = find_files(".")
        logger.info("Files inside the MEX directory")
        logger.info(data_paths)

        common_prefix = os.path.commonprefix(data_paths)
        MEX_prefix = common_prefix.split("/")[-1]
        if MEX_prefix != "":
            MEX_dir = common_prefix[:common_prefix.find(MEX_prefix)]
        else:
            MEX_dir = common_prefix
        logging.info(MEX_dir)
        logger.info("Searhing for MEX_dir {} using prefix = {}".format(MEX_dir, MEX_prefix))

        scanpy_obj = sc.read_10x_mtx(MEX_dir, gex_only = True, prefix = MEX_prefix)
        use_metacell_data_file = True
    else:
        raise ValueError("Input file format unrecognized! Please refer to the documentation for further assistance.")

    logger.info(input_file)
    logger.info("has been read to memory!")

    # Normalize the raw data
    logger.info("Normalizing the data...")
    logger.info('==========================================================')
    scanpy_obj = CLR_normalize(scanpy_obj)
    logger.info("Normalizing complete! ")

    annotation_slot = cluster_data_label
    if use_metacell_data_file:
        logger.info("Processing user-supplied cell annotation file")
        cell_annotation_df = pd.read_csv(cluster_data_label, sep = "\t", index_col = 0)
        annotation_slot = cell_annotation_df.columns[0]
        scanpy_obj.obs = scanpy_obj.obs.merge(cell_annotation_df, left_index = True, right_index = True)

    # Group by metacell 
    logger.info("Aggregating cells by clusters...")
    logger.info('==========================================================')
    metacell_exp = AverageExpression(scanpy_obj, annotation_slot)
    logger.info("Aggregating cells by clusters complete!")

    if chip_file:
        logger.info("Loading CHIP file to convert to Gene Symbol")
        logger.info('==========================================================')
        chip = read_chip(chip_file)
        metacell_exp = convert_to_gene_symbol(chip, metacell_exp)
        logger.info("Loaded CHIP file!\n")

    else:
        logger.info("The first 10 gene identifiers found in the expression matrix are: ")
        logger.info(metacell_exp.index[:10])
        logger.info("Chip file was not selected and the original gene identifiers are used without conversion. \
            If the gene identifiers are not in gene symbols, the outputted .gct file may be empty.")
    
    # Save the cluster expression matrix
    metacell_exp.to_csv("metacell_exp_converted.csv", sep=",", mode = "w")

    # Load the gene set database files
    logger.info(f"Loading gene set database files")
    logger.info('==========================================================')

    gs, gs_desc = read_gmts(gene_set_database_file)
    logger.info("Loaded gene set file!\n")

    logger.info("Number of gene sets loaded for scGSEA: {}".format(gs.shape[0]))

    logger.info("Running scGSEA...")
    logger.info('==========================================================')
    start = timer()
    scGSEA_scores = run_ssgsea_parallel(
        metacell_exp,
        gs,
        n_job = n_job,
        file_path = None
    )
    end = timer()

    scGSEA_scores.to_csv(output_file_name + ".csv", sep=",", mode = 'w')

    write_gct(scGSEA_scores, output_file_name, gs_desc)

    logger.info("We are done!")
    logger.info(f"scGSEA Runtime using {n_job} CPUs: ", humanfriendly.format_timespan(end - start))