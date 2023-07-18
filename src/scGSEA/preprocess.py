#!/usr/bin/env python3

from scgsea_helper import *
import argparse
import pandas as pd
import logging
import os
import scanpy as sc
import sys
import zipfile
import humanfriendly
from timeit import default_timer as timer

parser = argparse.ArgumentParser()

# ~~~~Module Required Arguments~~~~~ #
parser.add_argument("--input_file",
                    type=str,
                    help="RDS file to load",
                    default='False')

parser.add_argument("--cluster_data_label",
                    type=str,
                    help="Metadata label to use for aggregating cells",
                    default='seurat_clusters')

parser.add_argument("--output_file_name",
                    type=str,
                    help="filename to use for output files",
                    default = None)

# ~~~~Optional Arguments~~~~ #
parser.add_argument("--cluster_data_file",
                    type=str,
                    help="File containing cell annotation/grouping data",
                    default = None)

args = parser.parse_args()

Log_Format = "[%(asctime)s] %(message)s"
logging.basicConfig(level=logging.INFO, 
                    stream = sys.stdout, 
                    filemode = "w",
                    format = Log_Format,
                    datefmt='%Y-%m-%d %H:%M:%S')

logger = logging.getLogger()

################################################################################
#Begin Running the functions
################################################################################

logger.info('==========================================================')
logger.info("About to read the object named")

split_name = args.input_file.split(".")
file_extension = split_name[-1]    

# User supplies h5ad object with annotation information
if check_file_extension(args.input_file, "h5ad"):
    scanpy_obj = sc.read_h5ad(args.input_file)
    userInput = False

# User supplies raw counts data and annotation file
elif check_file_extensions(args.input_file, ["h5", "hdf5"]):
    scanpy_obj = sc.read_10x_h5(args.input_file, gex_only = True)
    userInput = True

elif check_file_extension(args.input_file, "zip"):
    with zipfile.ZipFile(args.input_file, 'r') as zip_ref:
        zip_ref.extractall("./data")

    data_paths = find_files(".")
    logger.info("Files inside the MEX directory")
    logger.info(data_paths)

    common_prefix = os.path.commonprefix(data_paths)

    logging.info(common_prefix)
    MEX_prefix = common_prefix.split("/")[-1]
    logging.info(MEX_prefix)
    if MEX_prefix != "":
        MEX_dir = common_prefix[:common_prefix.find(MEX_prefix)]
    else:
        MEX_dir = common_prefix
    logging.info(MEX_dir)
    logger.info("Searhing for MEX_dir {} using prefix = {}".format(MEX_dir, MEX_prefix))

    scanpy_obj = sc.read_10x_mtx(MEX_dir, gex_only = True, prefix = MEX_prefix)
    userInput = True
else:
    raise ValueError("Input file format unrecognized! Please refer to the documentation for further assistance.")

logger.info(args.input_file)
logger.info("has been read to memory!")

# Normalize the raw data
logger.info("Normalizing the data...")
logger.info('==========================================================')
scanpy_obj = CLR_normalize(scanpy_obj)
logger.info("Normalizing complete! ")

annotation_slot = args.cluster_data_label
if userInput:
    logger.info("Processing user-supplied cell annotation file")
    cell_annotation_df = pd.read_csv(args.cluster_data_file, sep = "\t", index_col = 0)
    annotation_slot = cell_annotation_df.columns[0]
    scanpy_obj.obs = scanpy_obj.obs.merge(cell_annotation_df, left_index = True, right_index = True)

# Group by metacell 
logger.info("Aggregating cells by clusters...")
logger.info('==========================================================')
aggregated_obj = AverageExpression(scanpy_obj, annotation_slot)
logger.info("Aggregating cells by clusters complete!")

logger.info("Saving metacell average expression profile...")
logger.info('==========================================================')
aggregated_obj.to_csv("cluster_expression_normalized.csv")
logger.info("cluster_expression_normalized.csv file saved!")
logger.info("Aggregating cells by clusters complete!")