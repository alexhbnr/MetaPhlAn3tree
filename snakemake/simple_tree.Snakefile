################################################################################
# Build gene trees using FastTree, and generate single tree using ASTRAL.
################################################################################

from glob import glob
import os
import sys

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
        "MetaPhlAn3tree.IQtree.tre"
    output:
        touch("done/simple_tree")

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
        species = [os.path.basename(i).replace("fasta", "faa") for i in glob("tmp/clean_aa/*.fasta")]
    run:
        phylophlan.concatenate(params.species, params.dir, output[0], False, True)

rule tree1:
    input:
        "tmp/markers.aln"
    output:
        "MetaPhlAn3tree.IQtree.tre"
    message: "Generate phylogenetic tree using IQtree"
    params:
        prefix = "tmp/MetaPhlAn3tree.IQtree"
    threads: 24
    shell:
        """
        iqtree -quiet -nt {threads} -m LG \
            -s {input} \
            -pre {params.prefix}
        cp {params.prefix}.treefile {output}
        """

