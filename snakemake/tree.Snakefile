################################################################################
#
#
# Alex Huebner, 
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
        "done/tree"

rule create_checkpoint:
    input: 
        expand("tmp/gene_tree2/{marker}.tre", marker=MARKERS)
    output:
        touch("done/tree")

rule gene_tree1:
    output:
        "tmp/gene_tree1/{marker}.tre"
    message: "Build initial gene tree for marker {wildcards.marker} using FastTree"
    params:
        aln = "tmp/fragmentary/{marker}.aln",
    resources:
        cores = 8
    threads: 8
    shell:
        """
        set +u
        source /projects1/users/huebner/miniconda3/etc/profile.d/conda.sh
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

rule resolve_polytomies:
    input:
        "tmp/gene_tree1/{marker}.tre"
    output:
        "tmp/gene_tree1_polytomies/{marker}.tre"
    message: "Resolve polytomies of tree of marker {wildcards.marker}"
    resources:
        cores = 1
    run:
        # Function phylophlan.resolve_polytomies_rec without the
        # overhead of the multiprocessing module
        tree = dendropy.Tree.get(path=input[0], schema="newick", preserve_underscores=True)
        tree.resolve_polytomies()
        tree.write(path=output[0], schema="newick")

rule gene_tree2:
    input:
        "tmp/gene_tree1_polytomies/{marker}.tre"
    output:
        "tmp/gene_tree2/{marker}.tre"
    message: "Redefine tree of marker {wildcards.marker} using RAxML"
    resources:
        cores = 9
    params:
        aln = "tmp/fragmentary/{marker}.aln",
        workdir = f"{os.path.abspath(config['tmpdir'])}/tmp/gene_tree2"
    threads: 8
    shell:
        """
        set +u
        source /projects1/users/huebner/miniconda3/etc/profile.d/conda.sh
        conda activate phylophlan
        set -u
        raxmlHPC-PTHREADS-SSE3 -m PROTCATLG -p 1989 \
            -t {input} \
            -w {params.workdir} \
            -s {params.aln} \
            -n {wildcards.marker}.tre \
            -T {threads}
        ln -s {params.workdir}/RAxML_bestTree.{wildcards.marker}.tre {output}
        """
