################################################################################
# Uncompress the representatitve genomes downloaded from GenBank
#
# Alex Huebner, 04/04/2020
################################################################################

workdir: config['tmpdir']


GENOMES, = glob_wildcards("fasta/{gcaid}.fna.gz")

rule all:
    input:
        "done/uncompress_fnas"

rule create_checkpoint:
    input: 
        expand("tmp/uncompressed/{gcaid}.fna", gcaid=GENOMES)
    output: touch("done/uncompress_fnas")

rule decompress_genomes:
    output:
        "tmp/uncompressed/{gcaid}.fna"
    message: "Uncompress genome {wildcards.gcaid}"
    params: 
        fna = "fasta/{gcaid}.fna.gz",
        condaenv = config['condaenv'] + '/etc/profile.d/conda.sh'
    threads: 4
    shell:
        """
        set +u
        source {params.condaenv}
        conda activate phylophlan
        set -u
        pigz -dck -p {threads} {params.fna} > {output} 
        """
