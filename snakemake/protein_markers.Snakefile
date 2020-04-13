################################################################################
# Clean the fake proteomes, align against the protein marker database of Segata
# et al. (2013) using DIAMOND blastp, and identify and extract marker gene
# sequences.
#
# Alex Huebner, 07/04/2020
################################################################################

import gzip
from Bio import SeqIO  # Biopython requires NumPy
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import phylophlan.phylophlan as phylophlan
from phylophlan.phylophlan import info, error

workdir: config['tmpdir']

GENOMES, = glob_wildcards("tmp/fake_proteomes/{gcaid}.faa")

rule all:
    input:
        "done/protein_markers"

rule create_checkpoint:
    input: 
        expand("tmp/markers_aa/{gcaid}.faa.gz", gcaid=GENOMES)
    output:
        touch("done/protein_markers")

rule clean_input_proteomes:
    output:
        "tmp/clean_aa/{gcaid}.fasta"
    message: "Clean proteome of genome {wildcards.gcaid}"
    resources:
        cores = 1
    params:
        faa = "tmp/fake_proteomes/{gcaid}.faa"
    threads: 1
    run:
        # Function phylophlan.fake_proteome_rec without the overhead
        # of the multiprocessing module
        inp, out = (params.faa, output[0])
        inp_clean, _ = os.path.splitext(os.path.basename(inp))
        info('Cleaning "{}"\n'.format(inp))

        # http://biopython.org/DIST/docs/api/Bio.Alphabet.IUPAC.ExtendedIUPACProtein-class.html
        # B = "Asx"; Aspartic acid (R) or Asparagine (N)
        # X = "Xxx"; Unknown or 'other' amino acid
        # Z = "Glx"; Glutamic acid (E) or Glutamine (Q)
        # J = "Xle"; Leucine (L) or Isoleucine (I), used in mass-spec (NMR)
        # U = "Sec"; Selenocysteine
        # O = "Pyl"; Pyrrolysine
        output = (SeqRecord(Seq(e[1].replace('B', 'X').replace('Z', 'X').replace('J', 'X').replace('U', 'X').replace('O', 'X')),
                            id='{}_{}'.format(inp_clean, counter), description='')
                    for counter, e in enumerate(SimpleFastaParser(open(inp))) if e[1])

        with open(out, 'w') as f:
            SeqIO.write(output, f, "fasta")

rule diamond_blastp:
    input:
        "tmp/clean_aa/{gcaid}.fasta"
    output:
        #temp("tmp/map_aa/{gcaid}.b6o.bkp")
        "tmp/map_aa/{gcaid}.b6o.bkp"
    message: "Align genome {wildcards.gcaid} against protein database using BLASTP"
    resources:
        cores = 8
    params:
        condaenv = config['condaenv'] + '/etc/profile.d/conda.sh',
        db = config['tmpdir'] + "/database/phylophlan/phylophlan"
    threads: 8
    shell:
        """
        set +u
        source {params.condaenv}
        conda activate phylophlan
        set -u
        diamond blastp --quiet \
                       --threads {threads} \
                       --outfmt 6 \
                       --more-sensitive \
                       --id 50 \
                       --max-hsps 35 \
                       -k 0 \
                       --query {input} \
                       --db {params.db} \
                       --out {output}
        """

rule gene_marker_selection:
    input:
        "tmp/map_aa/{gcaid}.b6o.bkp"
    output:
        "tmp/map_aa/{gcaid}.b6o"
    message: "Select gene markers for genome {wildcards.gcaid}"
    resources:
        cores = 1
    params:
        min_num_proteins = 1
    run:
        matches = phylophlan.best_hit(input[0], False)        
        if len(matches) >= params.min_num_proteins:
            with open(output[0], 'w') as f:
                f.write('{}\n'.format('\n'.join(['\t'.join(m) for m in matches])))

rule marker_gene_extraction:
    input:
        blast = "tmp/map_aa/{gcaid}.b6o",
        faa = "tmp/clean_aa/{gcaid}.fasta"
    output:
        "tmp/markers_aa/{gcaid}.faa.gz"
    message: "Extract marker genes for genome {wildcards.gcaid}"
    resources:
        cores = 1
    run:
        # Function phylophlan.gene_markers_extraction_rec without the overhead
        # of the multiprocessing module
        out_file, src_file, b6o_file, min_num_markers, frameshifts = ([output[0], input[1],
                                            input[0], 1, False])
        out_file_seq = []
        contig2marker2b6o = {}
        info('Extracting "{}"\n'.format(b6o_file))

        for l in open(b6o_file):
            row = l.strip().split('\t')
            contig = row[0]
            marker = row[1]
            start = int(row[2])
            end = int(row[3])
            rev = bool(int(row[4]))

            if (contig in contig2marker2b6o) and (marker in contig2marker2b6o[contig]):
                error('contig: {} and marker: {} already present into contig2marker2b6o'.format(contig, marker))

            if contig not in contig2marker2b6o:
                contig2marker2b6o[contig] = {}

            contig2marker2b6o[contig][marker] = (end, start, rev) if end < start else (start, end, rev)

        for record in SimpleFastaParser(open(src_file)):
            fid = record[0].split(' ')[0]

            if fid not in contig2marker2b6o:
                continue

            for marker in contig2marker2b6o[fid]:
                s, e, rev = contig2marker2b6o[fid][marker]
                idd = '{}_{}:'.format(fid, marker)
                seq = Seq(record[1][s - 1:e]) if (s - 1 >= 0) else None

                if not seq:  # skip empty sequences
                    continue

                if rev:
                    idd += 'c'
                    seq = seq.reverse_complement()

                out_file_seq.append(SeqRecord(seq, id='{}{}-{}'.format(idd, s, e), description=''))

                if frameshifts:
                    if not rev:
                        if record[1][s:e]:  # skip empty sequences
                            out_file_seq.append(SeqRecord(Seq(record[1][s:e]),
                                                            id='{}{}-{}'.format(idd, s + 1, e), description=''))

                        if record[1][s + 1:e]:  # skip empty sequences
                            out_file_seq.append(SeqRecord(Seq(record[1][s + 1:e]),
                                                            id='{}{}-{}'.format(idd, s + 2, e), description=''))
                    else:
                        if (s - 1 >= 0) and (e - 1 >= 0) and (record[1][s - 1:e - 1]):  # skip empty sequences
                            out_file_seq.append(SeqRecord(Seq(record[1][s - 1:e - 1]).reverse_complement(),
                                                            id='{}{}-{}'.format(idd, s, e - 1), description=''))

                        if (s - 1 >= 0) and (e - 2 >= 0) and (record[1][s - 1:e - 2]):  # skip empty sequences
                            out_file_seq.append(SeqRecord(Seq(record[1][s - 1:e - 2]).reverse_complement(),
                                                            id='{}{}-{}'.format(idd, s, e - 2), description=''))

        len_out_file_seq = int(len(out_file_seq) / 3) if frameshifts else len(out_file_seq)

        if out_file_seq and (len_out_file_seq >= min_num_markers):
            with gzip.open(out_file, 'wt') as f:
                SeqIO.write(out_file_seq, f, 'fasta')
