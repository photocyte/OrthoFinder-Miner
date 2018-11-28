seqkit replace -p " .+" -r "" uniprot_arabidopsis_medicago_NC-lyases.fa.gz | grep ">" | tr -d ">" > NC-lyases.txt

OFO="/lab/solexa_weng/Seq_data/Projects/Sophia_Xu/20180503_Galega_OrthoFinder/peptide_files/Results_May04/"
SL="NC-lyases.txt"
python3 ../orthofinder-miner.py --orthofinder_output $OFO --similarity_list $SL | sort | uniq > NC-lyase_OGs.txt
