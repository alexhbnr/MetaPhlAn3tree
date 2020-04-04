################################################################################
# Write config file from Argparse object
#
# Alex Huebner, 02/04/2020
################################################################################

import json
import os
import sys


def argparse_to_json(args):
    """Write selected Argparse arguments to config file for Snakemake."""
    # Subset options from Argparse
    options = {k: v for k, v in args.items()
               if k in ['tmpdir']}
    # Generate other parameters
    options['urllist'] = f"{options['tmpdir']}/repgenomes_urls.txt"
    if os.path.basename(os.environ['CONDA_PREFIX']) == "miniconda":
        options['condaenv'] = os.environ['CONDA_PREFIX']
    else:
        options['condaenv'] = os.environ['CONDA_PREFIX'] + "/../../"
    # Write JSON to file
    with open(f"{options['tmpdir']}/snakemake_config.json", "wt") as outfile:
        json.dump(options, outfile)


def argparse_to_clusterconfig(args, jsontemplate):
    """Write selected Argparse arguments to cluster config file for Snakemake."""
    # Subset options from Argparse
    options = {k: v for k, v in args.items()
               if k in ['tmpdir']}
    with open(jsontemplate, "rt") as jsonfile:
        snakemake_config = json.load(jsonfile)
    # Add the file paths for the log files
    logfileprefix = options['tmpdir'] + '/cluster_logs/scheduler.%j'
    snakemake_config['__default__']['out'] =  logfileprefix + '.out'
    snakemake_config['__default__']['err'] =  logfileprefix + '.err'
    with open(options['tmpdir'] + '/snakemake_cluster.json', "wt") as outfile:
        json.dump(snakemake_config, outfile)
