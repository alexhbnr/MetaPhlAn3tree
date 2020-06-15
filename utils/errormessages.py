import sys

CHECKPOINT_FILES = {'download_representative_genomes': ('download_genomes.Snakefile',
                                                        'download_representative_genomes',
                                                        'download_representative_genomes.log'),
                    'fake_proteomes': ('fake_proteomes.Snakefile',
                                       'fake_proteomes',
                                       'fake_proteomes.log'),
                    'protein_markers': ('protein_markers.Snakefile',
                                        'protein_markers',
                                        'protein_markers.log'),
                    'alignment': ('alignment.Snakefile',
                                  'alignment',
                                  'alignment.log'),
                    'tree': ('tree.Snakefile',
                             'tree',
                             'tree.log')}


def print_errormessage(analysis, tmpdir):
    """Generate error message for missing checkpoint files."""
    checkpoint = CHECKPOINT_FILES[analysis]
    print(f"The execution of the Snakefile '{checkpoint[0]}' was not successful"
          f". The expected checkpoint file '{tmpdir + '/done/' + checkpoint[1]}' "
          "was not created. Please check the log file of the pipeline, "
          f"'{tmpdir + '/logs/snakemake-' + checkpoint[2]}' for the problem.",
          file=sys.stderr)
