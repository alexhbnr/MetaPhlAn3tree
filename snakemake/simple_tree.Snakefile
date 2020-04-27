################################################################################
# Build gene trees using FastTree, and generate single tree using ASTRAL.
#
# Alex Huebner, 27/04/2020
################################################################################

import os
import sys

import dendropy
import phylophlan.phylophlan as phylophlan
from phylophlan.phylophlan import info, error

workdir: config['tmpdir']

# Load subsitution model
sub_mod = phylophlan.load_substitution_model(os.path.dirname(phylophlan.__file__) +
                                             "/phylophlan_substitution_models/phylophlan.tsv")

# Marker alignments
MARKERS, = glob_wildcards("tmp/fragmentary/{marker}.aln")
MARKERS = [marker for marker in MARKERS if os.path.getsize(f"tmp/fragmentary/{marker}.aln") > 0]

rule all:
    input:
        "done/simple_tree"

rule create_checkpoint:
    input: 
        "MetaPhlAn3tree.FastTree_Astral.tre"
    output:
        touch("done/simple_tree")

rule gene_tree1:
    output:
        "tmp/gene_tree1/{marker}.tre"
    message: "Build initial gene tree for marker {wildcards.marker} using FastTree"
    params:
        condaenv = config['condaenv'] + '/etc/profile.d/conda.sh',
        aln = "tmp/fragmentary/{marker}.aln",
    resources:
        cores = 8
    threads: 8
    shell:
        """
        set +u
        source {params.condaenv}
        conda activate phylophlan
        set -u
        export OMP_NUM_THREADS={threads}
        FastTreeMP -mlacc 2 \
                   -slownni \
                   -spr 4 \
                   -fastest \
                   -mlnni 4 \
                   -no2nd \
                   -lg \
                   -out {output} \
                   < {params.aln} 
        """

rule merging_gene_trees:
    input:
        expand("tmp/gene_tree1/{marker}.tre", marker=MARKERS)
    output:
        "tmp/gene_trees_FastTree.tre"
    message: "Concatenate all FastTree tree files into a single one"
    params:
        dir = "tmp/gene_tree1"
    run:
        phylophlan.merging_gene_trees(params.dir, output[0], verbose=True)

rule astral:
    input:
        "tmp/gene_trees_FastTree.tre"
    output:
        "MetaPhlAn3tree.FastTree_Astral.tre"
    message: "Summarise FastTree trees using Astral"
    resources:
        cores = 1
    params:
        condaenv = config['condaenv'] + '/etc/profile.d/conda.sh'
    threads: 1
    shell:
        """
        set +u
        source {params.condaenv}
        conda activate phylophlan
        set -u
        java -Xmx400g -jar ${{CONDA_PREFIX}}/bin/astral.5.7.1.jar -i {input} -o {output}
        """
