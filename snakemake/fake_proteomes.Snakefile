################################################################################
# Uncompress the representatitve genomes downloaded from GenBank
#
# Alex Huebner, 04/04/2020
################################################################################

import gzip

from Bio import SeqIO  # Biopython requires NumPy
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from phylophlan.phylophlan import info, error

workdir: config['tmpdir']


GENOMES, = glob_wildcards("fasta/{gcaid}.fna.gz")

rule all:
    input:
        "done/fake_proteomes"

rule create_checkpoint:
    input: 
        expand("tmp/fake_proteomes/{gcaid}.faa", gcaid=GENOMES)
    output:
        touch("done/fake_proteomes")

rule decompress_genomes:
    output:
        temp("tmp/uncompressed/{gcaid}.fna")
    message: "Uncompress genome {wildcards.gcaid}"
    params: 
        condaenv = config['condaenv'] + '/etc/profile.d/conda.sh',
        fna = "fasta/{gcaid}.fna.gz"
    resources:
        cores = 4
    threads: 4
    shell:
        """
        set +u
        source {params.condaenv}
        conda activate phylophlan
        set -u
        pigz -dck -p {threads} {params.fna} > {output} 
        """

rule diamond_blastx:
    input:
        "tmp/uncompressed/{gcaid}.fna"
    output:
        temp("tmp/map_dna/{gcaid}.b6o.bkp")
    message: "Align genome {wildcards.gcaid} against protein database using BLASTX"
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
        diamond blastx --quiet \
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
        "tmp/map_dna/{gcaid}.b6o.bkp"
    output:
        "tmp/map_dna/{gcaid}.b6o"
    message: "Select gene markers for genome {wildcards.gcaid}"
    resources:
        cores = 1
    params:
        min_num_proteins = 1
    run:
        matches = phylophlan.largest_cluster(input[0], True)        
        if len(matches) >= params.min_num_proteins:
            with open(output[0], 'w') as f:
                f.write('{}\n'.format('\n'.join(['\t'.join(m) for m in matches])))

rule marker_gene_extraction:
    input:
        blast = "tmp/map_dna/{gcaid}.b6o",
        fna = "tmp/uncompressed/{gcaid}.fna"
    output:
        "tmp/markers_dna/{gcaid}.fna.gz"
    message: "Extract marker genes for genome {wildcards.gcaid}"
    resources:
        cores = 1
    run:
        # Function phylophlan.gene_markers_extraction_rec without the overhead
        # of the multiprocessing module
        out_file, src_file, b6o_file, min_num_markers, frameshifts = ([output[0], input.fna,
                                            input.blast, 1, True])
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

rule fake_proteomes:
    input:
        "tmp/markers_dna/{gcaid}.fna.gz"
    output:
        "tmp/fake_proteomes/{gcaid}.faa"
    message: "Generate fake proteome for genome {wildcards.gcaid}"
    resources:
        cores = 1
    params:
        min_protein_length = 50
    threads: 1
    run:
        inp, out, min_len_protein = (input[0], output[0], params.min_protein_length)
        proteome = []
        info('Generating "{}"\n'.format(inp))
        for idd, seq in SimpleFastaParser(gzip.open(inp, "rt")):
            idd = idd.split(' ')[0]
            s, e = idd.split(':')[-1].split('-')
            if s.startswith('c'):
                s = s[1:]
            while (len(seq) % 3) != 0:
                seq = seq[:-1]
            seq_t = Seq.translate(Seq(seq), to_stop=True)
            if len(seq_t) >= min_len_protein:
                proteome.append(SeqRecord(seq_t, id=idd, description=''))
        with open(out, 'w') as f:
            SeqIO.write(proteome, f, 'fasta')
