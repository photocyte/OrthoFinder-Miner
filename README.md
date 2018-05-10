# OrthoFinder Miner

This is a software tool designed to take the output from [OrthoFinder](https://github.com/davidemms/OrthoFinder) and different expression analysis programs, apply different criteria, and then produce an final filtered list of genes for experimental follow up (e.g. cDNA cloning and recombinant expression).

## OrthoFinder run expectations 
* For known/reference species: Input the peptide FASTA file of the [Uniprot proteomes](https://www.uniprot.org/help/proteome) into OrthoFinder. E.g. the [Arabidopsis thaliana](https://www.uniprot.org/proteomes/UP000006548) Uniprot proteome.
* For unpublished species: Input the peptide FASTA file(s) produced from an [Transdecoder](https://github.com/TransDecoder/TransDecoder/wiki)(v5.2.0+) *in silico* translated *de novo* transcriptome into OrthoFinder.

Run this OrthoFinder analysis on your own computer / cluster.

## Criteria to filter by

* 1. Sequence similarity/homology to a provided gene
* 2. Greatly higher or lower expression in a given species versus another.
* 3. Not a conserved gene (aka, not in a "1-1-1" orthogroup)
* 4. Genes with a similar expression pattern to a target gene / gene expression similarity clustering.

## Program goals
* Parsing of Orthofinder results (Orthogroups.csv)
* Parsing of RSEM expression analysis results
* Parsing of Kallisto expression analysis results
