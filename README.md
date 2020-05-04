MetaPhlAn3tree
================

MetaPhlAn3tree is a pipeline that wraps the functionality of [PhyloPhlAn3](https://github.com/biobakery/phylophlan) (Asnicar et al. In press) in order to automatise the generation of a phylogenetic tree based on a set of 400 universal markers genes (Segata et al. 2013) for a set of microbial species detected by [MetaPhlAn3](https://github.com/biobakery/MetaPhlAn/tree/3.0).

Provided with the MetaPhlAn3 databases that was used to infer the taxonomic composition, it will extract and download a representative of each species in the database or a subset of specified taxa and uses these as input to PhyloPhlAn3. The PhyloPhlAn3 pipeline is not directly run itself but its functions are wrapped into Snakemake (Köster and Rahmann 2012) pipelines in order to allow for a better scalability across clusters.

The result is a phylogenetic tree that can be used as underlying tree for calculating distances between metagenomic samples using methods that incorporate the phylogenetic distances of taxa, e.g. UniFrac (Lozupone et al. 2007) or PhILR (Silverman et al. 2017).

Literature
==========

Asnicar, Francesco, Andrew Maltez Thomas, Francesco Beghini, Claudia Mengoni, Serena Manara, Paolo Manghi, Qiyun Zhu, et al. In press. “Precise Phylogenetic Analysis of Microbial Isolates and Genomes from Metagenomes Using Phylophlan 3.0.”

Köster, Johannes, and Sven Rahmann. 2012. “Snakemake—a Scalable Bioinformatics Workflow Engine.” *Bioinformatics* 28 (19). Oxford University Press: 2520–2.

Lozupone, Catherine A, Micah Hamady, Scott T Kelley, and Rob Knight. 2007. “Quantitative and Qualitative *β* Diversity Measures Lead to Different Insights into Factors That Structure Microbial Communities.” *Appl. Environ. Microbiol.* 73 (5). Am Soc Microbiol: 1576–85.

Segata, Nicola, Daniela Börnigen, Xochitl C Morgan, and Curtis Huttenhower. 2013. “PhyloPhlAn Is a New Method for Improved Phylogenetic and Taxonomic Placement of Microbes.” *Nature Communications* 4 (1). Nature Publishing Group: 1–11.

Silverman, Justin D, Alex D Washburne, Sayan Mukherjee, and Lawrence A David. 2017. “A Phylogenetic Transform Enhances Analysis of Compositional Microbiota Data.” *Elife* 6. eLife Sciences Publications Limited: e21887.
