#!/usr/bin/env python
###############################################################################
# 
#
# Alex Huebner, 02/04/2020
###############################################################################

import argparse
from glob import glob
import os
import re
import subprocess
import sys

from src import species


# Args = {'databasedir': '/projects1/users/huebner/miniconda3/envs/metaphlan2/bin/metaphlan_databases',
        # 'metaphlanversion': 'latest',
        # 'genbankurl': 'https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt',
        # 'tmpdir': '/projects1/users/huebner/tmp/metaphlantree'}


def main():
    """Generate phylogenetic tree based on species present in MetaPhlAn3.

       Provides the functionality to use MASH and PhyloPhlAn3 to generate a
       phylogenetic tree for MetaPhlAn2 version >= 2.9 that can estimate
       distances based on underlying phylogenetic information such as UniFrac
       or PhILR.
    """
    # Extract list of species from MetaPhlAn database pickle
    if os.path.isdir(Args['databasedir']):
        db_versions = [re.search(r'mpa_v([0-9]+)_.+.pkl',
                                 os.path.basename(db)).group(1)
                       for db in glob(f"{Args['databasedir']}/*.pkl")]
        db_versions.sort()
        if Args['metaphlanversion'] == 'latest':
            db_version = db_versions[-1]
        elif Args['metaphlanversion'].replace("v", "") in db_versions:
            db_version = Args['metaphlanversion'].replace("v", "")
        else:
            print(f"The database version {Args['metaphlanversion']} is not "
                  f"present in the database directory {Args['databasedir']}. "
                  "The following database versions are currently available: "
                  f"{', '.join(['v' + db for db in db_versions])}. Either "
                  "pick from the available or download the database using "
                  "MetaPhlAn.",
                  file=sys.stderr)
            sys.exit(1)
        species_cont = species.Species(db_version, Args['databasedir'])
        print("Extract the strain information from the MetaPhlAn database",
              file=sys.stderr)
        species_cont.extract_strains()
    else:
        print(f'The directory {Args["databasedir"]} does not exist. Specify a '
               'valid directory that contains the MetaPhlAn databases.',
               file=sys.stderr)
        sys.exit(1)

    # Download the overview of RefSeq genomes and join information with genomes
    print("Download the GCA assembly summary from NCBI",
          file=sys.stderr)
    if not os.path.isdir(Args['tmpdir']):
        os.makedirs(Args['tmpdir'])
    subprocess.run(f'cd  {Args["tmpdir"]} && wget -N -nH '
                    '--user-agent=Mozilla/5.0 --relative -r --no-parent '
                    '--reject "index.html*" --cut-dirs=2 -e robots=off '
                    f'{Args["genbankurl"]}', shell=True)
    print("Join the GCA assembly summary information with the MetaPhlAn "
          "database information", file=sys.stderr)
    species_cont.join_genbank(Args['tmpdir'] + "/" +
                              os.path.basename(Args['genbankurl'])) 
    print(f"Fetch information for missing genomes from NCBI Assembly directly.",
          file=sys.stderr)
    species_cont.get_missing_information()
    print("Determine the representative genomes of "
          f"{species_cont.genomes.shape[0]} genomes present in the database",
          file=sys.stderr)
    species_cont.determine_representative_genomes()
    print(f"Identified {len(species_cont.representative_genomes)} genomes. "
           "Download genomes from NCBI.", file=sys.stderr)
    species_cont.write_url_list(Args['tmpdir'] + "/repgenomes_urls.txt")

# Argument parser
Parser = argparse.ArgumentParser(description='Generate phylogenetic tree based '
                                 'on species present in MetaPhlAn3')
Parser.add_argument('-o', '--output', help='name of the tree file')
Parser.add_argument('--metaphlanversion', default="latest",
                    help='specify version of MetaPhlAn for which tree should '
                    'be constructed, e.g. v293 ["latest"]')
Parser.add_argument('--genbankurl', default='https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt',
                    help='URL to summary of the assembly status of genomes '
                    'present in NCBI GenBank '
                    '[https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt]')
Parser.add_argument('--tmpdir', required=True,
                    help='path to the folder for storing temporary output')
Parser.add_argument('--databasedir', required=True,
                    help='path to folder "metaphlan_databases"')
Args = vars(Parser.parse_args())

if __name__ == '__main__':
    main()
