################################################################################
# Build gene trees using FastTree, resolve polytomies using DendroPy, re-define
# trees using RAxML, and generate single tree using ASTRAL.
################################################################################

import os
import sys

import dendropy
from Bio import SeqIO
import phylophlan.phylophlan as phylophlan
from phylophlan.phylophlan import info, error

workdir: config['tmpdir']

# Load subsitution model
sub_mod = phylophlan.load_substitution_model(os.path.dirname(phylophlan.__file__) +
                                             "/phylophlan_substitution_models/phylophlan.tsv")

# Marker alignments
MARKERS, = glob_wildcards("tmp/fragmentary/{marker}.aln")
MARKERS = [marker for marker in MARKERS if os.path.getsize(f"tmp/fragmentary/{marker}.aln") > 0]
MARKERS = [marker for marker in MARKERS
           if len(list(SeqIO.parse(open(f"tmp/fragmentary/{marker}.aln", "rt"), "fasta"))) > 3]
print(f"Number of markers: {len(MARKERS)}")

localrules: astral

rule all:
    input:
        "done/tree"

rule create_checkpoint:
    input: 
        "MetaPhlAn3tree.RAxML_Astral.tre"
    output:
        touch("done/tree")

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
        cores = 8
    params:
        condaenv = config['condaenv'] + '/etc/profile.d/conda.sh',
        aln = "tmp/fragmentary/{marker}.aln",
        wdir = f"{os.path.abspath(config['tmpdir'])}/tmp/gene_tree2"
    threads: 8
    shell:
        """
        set +u
        source {params.condaenv}
        conda activate phylophlan
        set -u
        for f in info log result; do
            if [[ -f {params.wdir}/RAxML_${{f}}.{wildcards.marker}.tre ]]; then
                rm {params.wdir}/RAxML_${{f}}.{wildcards.marker}.tre
            fi
        done
        raxmlHPC-PTHREADS-SSE3 -m PROTCATLG -p 1989 \
            -t {input} \
            -w {params.wdir} \
            -s {params.aln} \
            -n {wildcards.marker}.tre \
            -T {threads}
        ln -s {params.wdir}/RAxML_bestTree.{wildcards.marker}.tre {output}
        """

rule merging_gene_trees:
    input:
        expand("tmp/gene_tree2/{marker}.tre", marker=MARKERS)
    output:
        "tmp/gene_trees_RAxML.tre"
    message: "Concatenate all RAxML tree files into a single one"
    params:
        dir = "tmp/gene_tree2"
    run:
        phylophlan.merging_gene_trees(params.dir, output[0], verbose=True)

rule astral:
    input:
        "tmp/gene_trees_RAxML.tre"
    output:
        "MetaPhlAn3tree.RAxML_Astral.tre"
    message: "Summarise RAxML trees using Astral"
    resources:
        cores = 1
    params:
        condaenv = config['condaenv'] + '/etc/profile.d/conda.sh',
    threads: 1
    shell:
        """
        set +u
        source {params.condaenv}
        conda activate phylophlan
        set -u
        java -Xmx400g -jar ${{CONDA_PREFIX}}/bin/astral.5.7.1.jar -i {input} -o {output}
        """
