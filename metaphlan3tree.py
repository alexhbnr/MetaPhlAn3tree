#!/usr/bin/env python
###############################################################################

__author__ = ('Alexander Huebner (alexhbnr@gmail.com)')
__version__ = '0.1'
__date__ = '28 April 2020'

import argparse
from glob import glob
import os
from pathlib import Path
import re
import subprocess
import sys

from ete3 import Tree
import pandas as pd

from src import database, species
from utils import config


def main():
    """Generate phylogenetic tree based on species present in MetaPhlAn3.

       Provides the functionality to use MASH and PhyloPhlAn3 to generate a
       phylogenetic tree for MetaPhlAn2 version >= 2.9 that can estimate
       distances based on underlying phylogenetic information such as UniFrac
       or PhILR.
    """
    # Configure snakemake execution
    print("Configure Snakemake for execution.")
    config.argparse_to_json(Args)
    if Args['snakemakedir'] is None:
        Args['snakemakedir'] = os.path.dirname(os.path.realpath(__file__)) + "/snakemake"
    if Args['clusterconfig_template'] is None:
        Args['clusterconfig_template'] = Args['snakemakedir'] + "/cluster_config_template.json"
    config.argparse_to_clusterconfig(Args,
                                     Args['clusterconfig_template'])
    os.makedirs(Args['tmpdir'] + "/logs", exist_ok=True)
    if not Args['local']:
        os.makedirs(Args['tmpdir'] + "/cluster_logs", exist_ok=True)
        Args['cluster_cmd'] = f"--cluster-config {Args['tmpdir']}/snakemake_cluster.json --cluster '{Args['cluster_cmd']}'"
        Args['cluster_cmd'] += f" --resources cores={Args['max_resources']}"
        Args['cluster_cmd'] += f" {Args['snakemake_args']}"
        print(Args['cluster_cmd'])
    else:
        Args['cluster_cmd'] = ''

    # Extract list of species from MetaPhlAn database pickle
    if (not os.path.isfile(Args['tmpdir'] + "/repgenomes_urls.txt") or
            Args['force']):
        print("1. Extract all species from MetaPhlAn database and prepare URLS "
              "for download from NCBI.", file=sys.stderr)
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
            print("\tExtract the strain information from the MetaPhlAn database",
                  file=sys.stderr)
            species_cont.extract_strains()
        else:
            print(f'The directory {Args["databasedir"]} does not exist. Specify a '
                  'valid directory that contains the MetaPhlAn databases.',
                  file=sys.stderr)
            sys.exit(1)

        # Download the overview of RefSeq genomes and join information with genomes
        print("\tDownload the GCA assembly summary from NCBI\n",
              file=sys.stderr)
        if not os.path.isdir(Args['tmpdir']):
            os.makedirs(Args['tmpdir'])
        subprocess.run(f'cd  {Args["tmpdir"]} && wget -N -nH '
                       '--user-agent=Mozilla/5.0 --relative -r --no-parent '
                       '--reject "index.html*" --cut-dirs=2 -e robots=off '
                       f'{Args["genbankurl"]}', shell=True)
        print("\tJoin the GCA assembly summary information with the MetaPhlAn "
              "database information", file=sys.stderr)
        species_cont.join_genbank(Args['tmpdir'] + "/" +
                                  os.path.basename(Args['genbankurl']))
        print(f"\tFetch information for missing genomes from NCBI Assembly directly.",
              file=sys.stderr)
        species_cont.get_missing_information()
        print("\tDetermine the representative genomes of "
              f"{species_cont.genomes.shape[0]} genomes present in the database",
              file=sys.stderr)
        if Args['taxnames'] is None:
            taxnames = []
        else:
            taxnames = [line.rstrip().replace("s__", "")
                        for line in open(Args['taxnames'], 'rt')]
        species_cont.subset_taxa(taxnames)
        species_cont.determine_representative_genomes()
        print(f"\tIdentified {len(species_cont.representative_genomes)} genomes.\n"
              "\tPrepare URL list for download of genomes from NCBI.", file=sys.stderr)
        species_cont.write_url_list(Args['tmpdir'] + "/repgenomes_urls.txt")
        species_cont.genomes.to_csv(Args['tmpdir'] + "/genomes.tsv", sep="\t",
                                    index=False)

    if (not os.path.isfile(Args['tmpdir'] + "/done/download_representative_genomes") or
            Args['force']):
        print("\nDownload the representative genomes from NCBI\n", file=sys.stderr)
        subprocess.run(f"snakemake -s {Args['snakemakedir']}/download_genomes.Snakefile "
                       f"--configfile {Args['tmpdir']}/snakemake_config.json "
                       f"{Args['cluster_cmd']} "
                       "--restart-times 5 -k "
                       f"-j {Args['nproc']}", shell=True,
                       stderr=open(Args['tmpdir'] +
                                   "/logs/snakemake-download_representative_genomes.log",
                                   "at"))

    if (not os.path.isfile(Args['tmpdir'] + "/done/install_marker_database") or
            Args['force']):
        print("Prepare marker gene database of Segata et al. (2013) for tree "
              "building", file=sys.stderr)
        database.write_superconfig_aa(Args['tmpdir'])
        phylophlan_config = database.load_check_config(Args['tmpdir'],
                                                       Args['nproc'])
        database.install_marker_database(Args['tmpdir'], phylophlan_config)
        Path(Args['tmpdir'] + '/done/install_marker_database').touch(exist_ok=True)

    if (not os.path.isfile(Args['tmpdir'] + "/done/fake_proteomes") or
            Args['force']):
        print("Uncompress FastA files downloaded from NCBI, align against the "
              "protein marker database of Segata et al. (2013) using DIAMOND "
              "blastx, identify and extract marker genes, and translate into "
              "amino acid sequences", file=sys.stderr)
        subprocess.run(f"snakemake -s {Args['snakemakedir']}/fake_proteomes.Snakefile "
                       f"--configfile {Args['tmpdir']}/snakemake_config.json "
                       f"{Args['cluster_cmd']} "
                       "--restart-times 5 "
                       f"-j {Args['nproc']}", shell=True,
                       stderr=open(Args['tmpdir'] +
                                   "/logs/snakemake-fake_proteomes.log", "at"))

    if (not os.path.isfile(Args['tmpdir'] + "/done/protein_markers") or
            Args['force']):
        print("Clean the fake proteomes, align against the protein marker "
              "database of Segata et al. (2013) using DIAMOND blastp, and "
              "identify and extract marker gene sequences.", file=sys.stderr)
        subprocess.run(f"snakemake -s {Args['snakemakedir']}/protein_markers.Snakefile "
                       f"--configfile {Args['tmpdir']}/snakemake_config.json "
                       f"{Args['cluster_cmd']} "
                       "--restart-times 5 "
                       f"-j {Args['nproc']}", shell=True,
                       stderr=open(Args['tmpdir'] +
                                   "/logs/snakemake-protein_markers.log", "at"))

    if (not os.path.isfile(Args['tmpdir'] + "/done/alignment") or
            Args['force']):
        print("Identify the marker genes that are present at least in four "
              "genomes and make alignments, trim non-variant sites and remove "
              "samples consisting of >= 90% gaps", file=sys.stderr)
        subprocess.run(f"snakemake -s {Args['snakemakedir']}/alignment.Snakefile "
                       f"--configfile {Args['tmpdir']}/snakemake_config.json "
                       f"{Args['cluster_cmd']} "
                       "--restart-times 5 "
                       f"-j {Args['nproc']}", shell=True,
                       stderr=open(Args['tmpdir'] +
                                   "/logs/snakemake-alignment.log", "at"))

    if Args['skip_redefining']:
        if (not os.path.isfile(Args['tmpdir'] + "/done/simple_tree") or
                Args['force']):
            print("Build gene trees using FastTree, and generate single tree "
                  "using ASTRAL.", file=sys.stderr)
            subprocess.run(f"snakemake -s {Args['snakemakedir']}/simple_tree.Snakefile "
                           f"--configfile {Args['tmpdir']}/snakemake_config.json "
                           f"{Args['cluster_cmd']} "
                           "--restart-times 5 "
                           f"-j {Args['nproc']}", shell=True,
                           stderr=open(Args['tmpdir'] +
                                       "/logs/snakemake-tree.log", "at"))
    else:
        if (not os.path.isfile(Args['tmpdir'] + "/done/tree") or
                Args['force']):
            print("Build gene trees using FastTree, resolve polytomies using "
                  "DendroPy, re-define trees using RAxML, and generate single "
                  "tree using ASTRAL.", file=sys.stderr)
            subprocess.run(f"snakemake -s {Args['snakemakedir']}/tree.Snakefile "
                           f"--configfile {Args['tmpdir']}/snakemake_config.json "
                           f"{Args['cluster_cmd']} "
                           "--restart-times 5 "
                           f"-j {Args['nproc']}", shell=True,
                           stderr=open(Args['tmpdir'] +
                                       "/logs/snakemake-tree.log", "at"))

    if (not os.path.isfile(Args['output']) or Args['force']):
        print("Annotate the tree with taxonomic information and write to output "
              "file.", file=sys.stderr)
        tree_annot_df = pd.read_csv(Args['tmpdir'] + "/genomes.tsv",
                                    sep="\t")[['kingdom', 'phylum', 'class',
                                               'order', 'family', 'genus',
                                               'species', 'GCAid']]

        # Prepare identifiers for annotation
        def concat(r):
            prefices = ['kingdom', 'phylum', 'class', 'order',
                        'family', 'genus', 'species', 't']
            return "|".join([f'{p[0]}__{u}' for p, u in zip(prefices, r)])
        tree_annot_df['label'] = tree_annot_df.apply(lambda r: concat(r), axis=1)
        tree_annot_df = tree_annot_df.set_index(['GCAid'])

        if Args['skip_redefining']:
            treefn = Args['tmpdir'] + "/MetaPhlAn3tree.FastTree_Astral.tre"
        else:
            treefn = Args['tmpdir'] + "/MetaPhlAn3tree.RAxML_Astral.tre"
        tre = Tree(open(treefn, "rt").readline())
        for l in tre.iter_leaves():
            l.name = tree_annot_df.loc[l.name.replace(".faa", ""), 'label']
        tre.write(outfile=Args['output'])
    else:
        print(f"The tree output file {Args['output']} already exists and the "
              "option '--force' has not been enabled to re-run it. Activate "
              "this option to re-run all steps or remove the output file prior "
              "to running the script.", file=sys.stderr)
        sys.exit(1)

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
Parser.add_argument('--taxnames',
                    help='subset tree building to species names listed in file;'
                    ' one species name per line')
Parser.add_argument('--skip_redefining', action='store_true',
                    help='skip removing polytomies and re-defining the gene '
                    'trees using RAxML')
Parser.add_argument('--snakemakedir',
                    help='folder with SnakeMake workflows [./snakemake]')
Parser.add_argument('--local', action='store_true',
                    help='do not submit PhyloPhlAn steps to cluster but run locally')
Parser.add_argument('--clusterconfig_template',
                    help='path to the cluster-config template JSON file '
                    '[./snakemake/cluster_config_template.json]')
Parser.add_argument('--cluster_cmd', default='sbatch --mem {cluster.mem} '
                                             '-p {cluster.partition} '
                                             '-t {cluster.time} '
                                             '-o {cluster.out} '
                                             '-e {cluster.err} -n {threads}',
                    help="command provided to option '--cluster' of snakemake")
Parser.add_argument('--max_resources', default=160, type=int,
                    help='maximum number of CPUs to be used on cluster [160]')
Parser.add_argument('--snakemake_args', default='',
                    help='additional arguments to pass to Snakemake')
Parser.add_argument('--nproc', default=8,
                    help='number of parallel jobs to submit to cluster; if '
                         '--local, number of maximum processors to run '
                         'Snakemake with [8]')
Parser.add_argument('--force', action='store_true',
                    help='ignore checkpoints and re-run all steps')
Args = vars(Parser.parse_args())

if __name__ == '__main__':
    main()
