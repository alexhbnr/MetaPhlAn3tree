###############################################################################
# Download the PhyloPhlAn database used by Segata et al. (2013)
#
# Alex Huebner, 04/04/2020
###############################################################################

import os
import subprocess

import pandas as pd
import phylophlan.phylophlan as phylophlan


def write_superconfig_aa(tmpdir):
    """Generate supertree_aa config."""
    subprocess.run("phylophlan_write_config_file "
                   f"-o {tmpdir}/supertree_aa.cfg "
                   "-d a "
                   "--db_aa diamond "
                   "--map_dna diamond "
                   "--map_aa diamond "
                   "--msa mafft "
                   "--trim trimal "
                   "--gene_tree1 fasttree "
                   "--gene_tree2 raxml "
                   "--tree1 astral "
                   "--overwrite "
                   "--verbose", shell=True)


def load_check_config(tmpdir, nproc):
    """Load and check supertree_aa config."""
    configs = phylophlan.read_configs(tmpdir + "/supertree_aa.cfg",
                                      verbose=True)
    phylophlan.check_configs(configs, verbose=True)
    phylophlan.check_dependencies(configs, nproc, verbose=True)
    return configs


def install_marker_database(tmpdir, config):
    """Install marker database from Segata et al. (2013)."""
    if not os.path.isdir(tmpdir + "/database"):
        os.makedirs(tmpdir + "/database")
    database_outputfn = "".join([tmpdir, "/database/",
                                 os.path.basename(phylophlan.DATABASE_DOWNLOAD_URL) \
                                 .replace('?dl=1', '')])
    phylophlan.download(phylophlan.DATABASE_DOWNLOAD_URL, database_outputfn,
                        overwrite=False, verbose=True)
    database_urls = pd.read_csv(database_outputfn, sep="\t",
                                index_col=['#database_name']) \
        .loc['phylophlan']
    phylophlan.download_and_unpack_db('phylophlan',
                                      database_urls['database_url'],
                                      database_urls['database_md5'],
                                      tmpdir + "/database",
                                      update=False, verbose=True)
    phylophlan.check_database('phylophlan', tmpdir + "/database", verbose=True)
    db_type, db_dna, db_aa = phylophlan.init_database('phylophlan', tmpdir + "/database",
                                           'a', config, 'db_dna', 'db_aa',
                                           verbose=True)


