# OrthoFinder Miner

This is a software tool designed to take the output from [OrthoFinder](https://github.com/davidemms/OrthoFinder) and apply different criteria, to produce an eventual list of genes for experimental follow up.

## Criteria to filter by

* 1. Sequence similarity/homology to a provided gene
* 2. Greatly higher or lower expression in a given species versus another.
* 3. Not a conserved gene (aka, not in a "1-1-1" orthogroup)
* 4. Genes with a similar 

## Goals

* Parsing of RSEM expression analysis results
* Parsing of Kallisto expression analysis results
