################################################################################
# Identify the marker genes that are present at least in four genomes and make
# alignments, trim non-variant sites and remove samples consisting of >= 90%
# gaps
################################################################################

from collections import Counter
import glob
import gzip
import os
from pathlib import Path
import tqdm

from Bio import SeqIO  # Biopython requires NumPy
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
import numpy as np
from skbio import TabularMSA, Protein

import phylophlan.phylophlan as phylophlan
from phylophlan.phylophlan import info, error

workdir: config['tmpdir']

def evaluate_no_markers(wildcards):
    """Evaluates the number of markers that were present at least four
       genomes."""
    checkpoint_output = checkpoints.inputs2markers.get(**wildcards).output[0]
    return expand("tmp/fragmentary/{marker}.aln",
                  marker=glob_wildcards(os.path.join(checkpoint_output, "{marker}.faa")).marker)

rule all:
    input:
        "done/alignment"

rule create_checkpoint:
    input:
        evaluate_no_markers
    output:
        touch("done/alignment")

checkpoint inputs2markers:
    output:
        directory("tmp/markers")
    message: "Identify marker genes present in input genomes"
    resources:
        cores = 1
    params:
        input_folder = "tmp/markers_aa",
        output_folder = "tmp/markers",
        extension = ".faa.gz",
        output_extension = ".faa",
        min_num_entries = 4
    threads: 1
    run:
        # Function phylophlan.inputs2markers with enabling the reading of
        # gzipped FastA files
        os.makedirs(params.output_folder, exist_ok=False)
        markers2inputs = {}

        for i, f in enumerate(glob.iglob(os.path.join(params.input_folder,
                                                      '*' + params.extension))):
            if i % 100 == 0 and i > 0:
                print(i, end='\r')
            inp, _ = os.path.splitext(os.path.basename(f))

            for idd, seq in SimpleFastaParser(gzip.open(f, 'rt')):
                marker = idd.split(' ')[0].split(':')[0].split('_')[-1]

                if marker in markers2inputs:
                    markers2inputs[marker].append(SeqRecord(Seq(seq), id=inp, description=''))
                else:
                    markers2inputs[marker] = [SeqRecord(Seq(seq), id=inp, description='')]

        for marker, sequences in markers2inputs.items():
            if len(sequences) >= params.min_num_entries:
                with open(os.path.join(params.output_folder, marker + params.output_extension), 'w') as f:
                    SeqIO.write(sequences, f, 'fasta')
            elif verbose:
                info('"{}" discarded, not enough inputs ({}/{})\n'.format(marker, len(sequences), params.min_num_entries))

rule msas:
    input:
        "tmp/markers/{marker}.faa"
    output:
        "tmp/msas/{marker}.aln"
    message: "Align genomes for proteome marker {wildcards.marker}"
    resources:
        cores = 1
    threads: 4
    params:
        condaenv = config['condaenv'] + '/etc/profile.d/conda.sh'
    shell:
        """
        set +u
        source {params.condaenv}
        conda activate phylophlan
        set -u
        mafft --quiet \
              --anysymbol \
              --thread {threads} \
              --auto \
              {input} > {output}
        """

rule trim_non_variant:
    input:
        "tmp/msas/{marker}.aln"
    output:
        temp("tmp/trim_non_variant/{marker}.aln")
    message: "Trim non-variant sites from marker {wildcards.marker}"
    resources:
        cores = 1
    params:
        non_variant_threshold = 0.99
    run:
        # Function phylophlan.trim_not_variant_rec without the overhead
        # of the multiprocessing module
        inp, out, thr, verbose = [input[0], output[0],
                                  params.non_variant_threshold, True]
        info('Trimming not variant "{}"\n'.format(inp))
        aln = TabularMSA.read(inp, constructor=Protein, lowercase=False)

        def evaluate_nonvariant(site):
            sitefreq = site.frequencies()
            if "-" in sitefreq:
                del sitefreq['-']
            nongap_samples = sum(sitefreq.values())
            remove_site = [(count / nongap_samples) >= params.non_variant_threshold
                           for aa, count in sitefreq.items()]
            return any(remove_site)

        nonvariant_sites = [evaluate_nonvariant(aln.iloc[:, i])
                            for i in tqdm.tqdm(range(aln.shape[1]))]
        if np.sum(nonvariant_sites) < aln.shape[1]:
            aln.loc[:, nonvariant_sites].write(out, "fasta")
        elif verbose:
            info('"{}" discarded because no columns retained while removing not variant sites (thr: {})\n'.format(inp, thr))

rule remove_fragmentary_entries:
    input:
        "tmp/trim_non_variant/{marker}.aln"
    output:
        "tmp/fragmentary/{marker}.aln"
    message: "Remove fragmentary sites from marker {wildcards.marker}"
    resources:
        cores = 1
    params:
        fragmentary_threshold = 0.9,
        min_num_entries = 4
    run:
        # Function phylophlan.remove_fragmentary_entries_rec without the
        # overhead of the multiprocessing module
        inp, out, frag_thr, min_num_entries, verbose = [input[0], output[0],
                                                        params.fragmentary_threshold,
                                                        params.min_num_entries,
                                                        True]
        info('Fragmentary "{}"\n'.format(inp))
        inp_aln = TabularMSA.read(inp, constructor=Protein, lowercase=False)
        gap_frequencies = inp_aln.gap_frequencies(axis="position") / inp_aln.shape[1]

        if np.sum(gap_frequencies < frag_thr) > 0:
            inp_aln[gap_frequencies < frag_thr, :].write(out, "fasta")
        else:
            Path(out).touch()
