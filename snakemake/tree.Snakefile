################################################################################
# Build gene trees using FastTree, and redefine using RAxML.
################################################################################

from glob import glob
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

localrules: link_alignments

rule all:
    input:
        "done/tree"

rule create_checkpoint:
    input: 
        "MetaPhlAn3tree.RAxML.tre"
    output:
        touch("done/tree")

rule link_alignments:
    output:
        "tmp/final_alignments/{marker}.aln"
    message: "Link alignment of marker {wildcards.marker} to folder for contatenation"
    params:
        aln = "tmp/fragmentary/{marker}.aln"
    shell:
        "ln -s ${{PWD}}/{params.aln} {output}"

rule concatenate_alignments:
    input:
        expand("tmp/final_alignments/{marker}.aln", marker=MARKERS)
    output:
        "tmp/markers.aln"
    message: "Concatenate all markers to single alignment"
    params:
        dir = "tmp/final_alignments/",
        species = [os.path.basename(i).replace(".gz", "") for i in glob("tmp/markers_aa/*.faa.gz")]
    run:
        phylophlan.concatenate(params.species, params.dir, output[0], False, True)

rule tree1:
    input:
        "tmp/markers.aln"
    output:
        "tmp/MetaPhlAn3tree.FastTree.tre"
    message: "Build initial tree using FastTree"
    params:
        condaenv = config['condaenv'] + '/etc/profile.d/conda.sh',
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
                   < {input} 
        """

rule resolve_polytomies:
    input:
        "tmp/MetaPhlAn3tree.FastTree.tre"
    output:
        "tmp/MetaPhlAn3tree.FastTree_nopolytomies.tre"
    message: "Resolve polytomies of tree"
    resources:
        cores = 1
    run:
        # Function phylophlan.resolve_polytomies_rec without the
        # overhead of the multiprocessing module
        tree = dendropy.Tree.get(path=input[0], schema="newick", preserve_underscores=True)
        tree.resolve_polytomies()
        tree.write(path=output[0], schema="newick")

rule tree2:
    input:
        aln = "tmp/markers.aln",
        tre = "tmp/MetaPhlAn3tree.FastTree_nopolytomies.tre"
    output:
        "MetaPhlAn3tree.RAxML.tre"
    message: "Redefine tree using RAxML"
    resources:
        cores = 16
    params:
        condaenv = config['condaenv'] + '/etc/profile.d/conda.sh',
        wdir = f"{os.path.abspath(config['tmpdir'])}/tmp/tree2"
    threads: 16
    shell:
        """
        set +u
        source {params.condaenv}
        conda activate phylophlan
        set -u
        mkdir -p {params.wdir}
        for f in info log result; do
            if [[ -f {params.wdir}/RAxML_${{f}}.concatenated.tre ]]; then
                rm {params.wdir}/RAxML_${{f}}.concatenated.tre
            fi
        done
        raxmlHPC-PTHREADS-SSE3 -m PROTCATLG -p 1989 \
            -t {input.tre} \
            -w {params.wdir} \
            -s {input.aln} \
            -n concatenated.tre \
            -T {threads}
        ln -s {params.wdir}/RAxML_bestTree.concatenated.tre {output}
        """
