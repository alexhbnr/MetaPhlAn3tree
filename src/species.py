###############################################################################
# Parse the list of species present in the MetaPhlAn database and select
# representative genomes for tree construction.
#
# Alex Huebner, 02/04/2020
###############################################################################

import bz2
import os
import pickle
import re
import sys

from Bio import Entrez
import pandas as pd
from pandas.api.types import CategoricalDtype


# Standard regex searches
extract_gcaid = re.compile(r'.+\|t__((GCA|GCF|PRJNA)_*[0-9]+)$')
extract_taxonomy = re.compile(r'(.+)\|t__(GCA|GCF|PRJNA)_*[0-9]+$')
split_rank = re.compile(r'([a-z])__(.+)')

# Taxonomy abbreviations
tax_abbr = {'k': 'kingdom',
            'p': 'phylum',
            'c': 'class',
            'o': 'order',
            'f': 'family',
            'g': 'genus',
            's': 'species'}

# Entrez
Entrez.email = "A.N.Other@example.com"


class Species:
    """Container for information on representative genomes."""

    def __init__(self, db_version, db_dir):
        self.db_version = db_version
        self.db_dir = db_dir
        self.pkl_fn = f'{self.db_dir}/mpa_v{self.db_version}_CHOCOPhlAn_201901.pkl'

    def extract_strains(self):
        """Extract genome IDs of strains present in database file."""
        def parse_strain_info(strainname, strainnumbers):
            # Parse strain information by extracting GCA id
            gca_id = extract_gcaid.search(strainname).group(1)
            taxid = int(strainnumbers[0].split("|")[-1])
            taxonomy = pd.DataFrame([split_rank.search(t).groups()
                                     for t in extract_taxonomy.search(strainname).
                                     group(1).split("|")], columns=["abbr", "unit"])
            taxonomy['rank'] = taxonomy['abbr'].map(tax_abbr)
            return taxonomy.drop(['abbr'], axis=1).set_index('rank').T \
                .assign(GCAid=gca_id) \
                .assign(taxid=taxid) \
                .assign(label=strainname) \
                .assign(genomelength=strainnumbers[1])[['GCAid', 'label',
                                                        'taxid', 'genomelength'] +
                                                       list(tax_abbr.values())] \
                .reset_index(drop=True)

        with bz2.open(self.pkl_fn, "rb") as bzfile:
            strains = pickle.load(bzfile)['taxonomy']
        genome_df = [parse_strain_info(name, numbers)
                     for name, numbers in strains.items()
                     if not name.startswith("k__Vir") ]
        self.genomes = pd.concat(genome_df)

    def join_genbank(self, genbankfn):
        """Add information from GenBank about assembly status."""
        genbank = pd.read_csv(genbankfn, sep="\t", skiprows=1,
                              usecols=['# assembly_accession', 'refseq_category',
                                       'taxid', 'species_taxid', 'assembly_level',
                                       'release_type', 'genome_rep', 'ftp_path']) \
            .rename({'taxid': 'strain_taxid'}, axis=1)
        genbank['GCAid'] = genbank['# assembly_accession'].str.replace(r"\.[0-9]+$", "")
        self.genomes = self.genomes.merge(genbank, how="left", on=['GCAid'])

    def get_missing_information(self):
        """Fill in missing information from NCBI Assembly."""
        # Identify missing entries from GenBank table
        missing_gcas = self.genomes.loc[(self.genomes['assembly_level'].isnull() &
                                         ~self.genomes['GCAid'].str.startswith("PRJ")),
                                        'GCAid'].tolist()
        if len(missing_gcas) > 0:
            print(f"\t\t{len(missing_gcas)} genomes have missing information. Search "
                  "on NCBI Assembly for them.", file=sys.stderr)
            # Retrieve for each missing entry from NCBI Assembly
            gca_entries = []
            for gcaid in missing_gcas:
                eterm_handle = Entrez.esearch(db="assembly", term=gcaid, report="full")
                eterm_record = Entrez.read(eterm_handle)
                esummary_handle = Entrez.esummary(db="assembly", id=eterm_record['IdList'][0], report='full')
                esummary_record = Entrez.read(esummary_handle, validate=False)
                gca_entries.append([esummary_record['DocumentSummarySet'] \
                                    ['DocumentSummary'][0][f]
                                    for f in ['AssemblyAccession', 'RefSeq_category',
                                              'Taxid', 'SpeciesTaxid', 'AssemblyStatus',
                                              'ReleaseType',
                                              'PartialGenomeRepresentation',
                                              'FtpPath_GenBank']])
            # Generate Pandas DataFrame for merging
            gca_dfs = pd.DataFrame(gca_entries, columns=['# assembly_accession',
                                                         'refseq_category',
                                                         'strain_taxid', 'species_taxid',
                                                         'assembly_level', 'release_type',
                                                         'genome_rep', 'ftp_path']) \
                    .assign(GCAid=missing_gcas)
            gca_dfs['genome_rep'] = gca_dfs['genome_rep'].replace({'false': 'Full',
                                                                   'true': 'Partial'})
            # Replace missing values
            missing_genomes = self.genomes.loc[self.genomes['assembly_level'].isnull(),
                                               self.genomes.columns.tolist()[:11]].copy()
            missing_genomes = missing_genomes.merge(gca_dfs, how="left", on=['GCAid'])
            self.genomes = pd.concat([self.genomes.loc[~self.genomes['assembly_level'].isnull()],
                                      missing_genomes])
        else:
            print("\t\tNo genomes have missing information. Continue.", file=sys.stderr)

    def subset_taxa(self, taxalist):
        """Subset representative genomes to selected ones for downloading."""
        if len(taxalist) > 0:
            self.genomes_set = self.genomes.loc[self.genomes['label'].str \
                                                .extract(r'.+\|(s__.+)\|t__GC[AF]_[0-9]+',
                                                         expand=False).isin(taxalist)]
            skipped_taxa = [t for t in taxalist
                            if t.replace("s__", "") not in self.genomes_set['species'].tolist()]
            print(f"\tFrom {len(taxalist)} species names in the specified list "
                  f"of taxa, {self.genomes_set['species'].unique().shape[0]} were "
                  f"found in MetaPhlAn database. The following species names "
                  f"are not present in the database: {', '.join(skipped_taxa)}.",
                  file=sys.stderr)
            if self.genomes_set.shape[0] > 0:
                print(f"\tContinue with remaining {len(taxalist) - len(skipped_taxa)}"
                      " species.", file=sys.stderr)
            else:
                print("\tNo species left. Select species present in the "
                      "database.", file=sys.stderr)
                sys.exit(1)
        else:
            self.genomes_set = self.genomes

    def determine_representative_genomes(self):
        """Determine the representative genome by evaluating assembly info."""
        assembly_stats = self.genomes_set[['GCAid', 'taxid', 'assembly_level',
                                           'release_type', 'genome_rep']].copy() \
            .reset_index(drop=True)

        # Convert assembly columns into ranked categorial data types
        asslevel_type = CategoricalDtype(categories=['Complete Genome',
                                                     'Chromosome', 'Scaffold',
                                                     'Contig'], ordered=True)
        assembly_stats.loc[:, ('assembly_level')] = assembly_stats.loc[:, ('assembly_level')] \
            .astype(asslevel_type) \
            .cat.codes
        release_type = CategoricalDtype(categories=['Major', 'Minor'],
                                        ordered=True)
        assembly_stats.loc[:, ('release_type')] = assembly_stats.loc[:, ('release_type')] \
            .astype(release_type) \
            .cat.codes
        genrep_type = CategoricalDtype(categories=['Full', 'Partial'],
                                       ordered=True)
        assembly_stats.loc[:, ('genome_rep')] = assembly_stats.loc[:, ('genome_rep')] \
            .astype(genrep_type) \
            .cat.codes
        # Calculate overall rank
        assembly_stats['rank'] = assembly_stats.loc[:, ~assembly_stats.columns
                                                    .isin(['taxid', 'GCAid'])].sum(axis=1)
        # Subset genome with lowest rank per taxid and sample randomly multiple
        # with equal rank
        self.representative_genomes = assembly_stats.iloc[assembly_stats
                                                          .groupby(['taxid'])
                                                          ['rank'].idxmin()] \
                                                    .groupby(['taxid']) \
                                                    .apply(lambda d: d.sample(1))['GCAid'] \
                                                    .unique().tolist()

    def write_url_list(self, outfn):
        """Write list of links to genomes for representative genomes to file."""
        def generate_url(ftppath):
            return f'{ftppath}/{os.path.basename(ftppath)}_genomic.fna.gz'

        link_df = self.genomes.loc[self.genomes['GCAid'].isin(self.representative_genomes),
                                   ('GCAid', 'ftp_path')].drop_duplicates().copy()

        empty_ftp_gcaids = link_df.loc[link_df['ftp_path'] == ""]['GCAid'].tolist()
        if len(empty_ftp_gcaids) > 0:
            print(f"\t\t{len(empty_ftp_gcaids)} genomes with a mixing FTP link. "
                  "Fetch from NCBI Assembly", file=sys.stderr)
            ftp_entries = []
            for gcaid in empty_ftp_gcaids:
                eterm_handle = Entrez.esearch(db="assembly", term=gcaid, report="full")
                eterm_record = Entrez.read(eterm_handle)
                esummary_handle = Entrez.esummary(db="assembly", id=eterm_record['IdList'][0], report='full')
                esummary_record = Entrez.read(esummary_handle, validate=False)
                ftp_entries.append(esummary_record['DocumentSummarySet'] \
                                   ['DocumentSummary'][0]['FtpPath_RefSeq'])
            link_df.loc[link_df['ftp_path'] == "", 'ftp_path'] = ftp_entries
            print(f"\t\tSkipping {link_df.loc[link_df['ftp_path'] == ''].shape[0]} "
                  f"genomes ({', '.join(link_df.loc[link_df['ftp_path'] == '', 'GCAid'].tolist())})"
                  " from analysis because of no FTP link information.", file=sys.stderr)

        link_df = link_df.loc[link_df['ftp_path'] != ""]
        link_df['url'] = link_df['ftp_path'].map(generate_url)
        link_df[['GCAid', 'url']].to_csv(outfn, sep="\t", index=False)
