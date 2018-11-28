# OrthoFinder Miner

This is a comparative genomic software tool designed to take the output from [OrthoFinder](https://github.com/davidemms/OrthoFinder) and different expression analysis programs (e.g. [RSEM](https://deweylab.github.io/RSEM/) or [Kallisto](https://pachterlab.github.io/kallisto/)), apply different criteria, and then produce an final filtered list of genes for experimental follow up (e.g. by cDNA cloning and recombinant expression).

## OrthoFinder run expectations 
* For known/reference species: Input the peptide FASTA file of the [Uniprot proteomes](https://www.uniprot.org/help/proteome) into OrthoFinder. E.g. the [*Arabidopsis thaliana*](https://www.uniprot.org/proteomes/UP000006548) Uniprot proteome. Prefiltering down to a single gene per gene-isoform cluster, such as with [this](https://github.com/photocyte/filter_uniprot_to_best_isoform) script might be a good idea to help with downstream intepretation.
* For unpublished species: Input the peptide FASTA file(s) produced from an [Transdecoder](https://github.com/TransDecoder/TransDecoder/wiki)(v5.2.0+) *in silico* translated *de novo* transcriptome into OrthoFinder.

Run this OrthoFinder analysis on your own computer / cluster.

For ease of readability, it is recommended that you filter the FASTA record descriptions away from the Transdecoder produced FASTSA file.  [Seqkit](https://github.com/shenwei356/seqkit) can accomplish this:

```seqkit replace -p " .+" -r "" INPUT.fa > OUTPUT.fa```

## Criteria to filter by

* 1. Sequence similarity/homology to a provided gene via being in the same Orthogroup **(Now implemented)**
* 2. Greatly higher or lower expression in a given species versus another.
* 3. Not a conserved gene (aka, not in a "1-1-1" orthogroup). Also known as a "Reciprocally direct orthogroup" or RDOG **(Now implemented)**
* 4. Genes with a similar expression pattern to a target gene / gene expression similarity clustering.

## Program goals
* Parsing of Orthofinder results & metadata (Orthogroups.csv, SpeciesIDs.txt, SequenceIDs.txt)
* Parsing of RSEM expression analysis results
* Parsing of Kallisto expression analysis results
