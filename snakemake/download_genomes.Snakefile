################################################################################
# Download representatitve genomes from GenBank
################################################################################

import os

import pandas as pd

workdir: config['tmpdir']


GENOMES = pd.read_csv(config['urllist'], sep="\t", index_col=['GCAid'])

rule all:
    input:
        "done/download_representative_genomes"

rule create_checkpoint:
    input: 
        expand("fasta/{gcaid}.fna.gz", gcaid=GENOMES.index.values)
    output: touch("done/download_representative_genomes")

rule download_genome:
    output:
        "fasta/{gcaid}.fna.gz"
    message: "Download genome {wildcards.gcaid}"
    params: 
        url = lambda wildcards: GENOMES.loc[wildcards.gcaid, 'url'],
        tmpdir = config['tmpdir']
    shell:
        """
        cd fasta && \
        wget -N -nH --user-agent=Mozilla/5.0 --relative -r --no-parent \
            --reject "index.html*" --cut-dirs=2 --cut-dirs 7 \
            {params.url}
        mv $(basename {params.url}) {wildcards.gcaid}.fna.gz
        """
