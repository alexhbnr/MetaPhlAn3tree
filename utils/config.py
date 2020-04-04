################################################################################
# Write config file from Argparse object
#
# Alex Huebner, 02/04/2020
################################################################################

import json
import sys


def argparse_to_json(args):
    """Write select Argparse arguments to config file for Snakemake."""
    # Subset options from Argparse
    options = {k: v for k, v in args.items()
               if k in ['tmpdir']}
    # Generate other parameters
    options['urllist'] = f"{options['tmpdir']}/repgenomes_urls.txt"
    # Write JSON to file
    with open(f"{options['tmpdir']}/snakemake_config.json", "wt") as outfile:
        json.dump(options, outfile)



