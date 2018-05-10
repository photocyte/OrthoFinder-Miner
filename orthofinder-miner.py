import Bio
import argparse
import glob

parser = argparse.ArgumentParser(description='Proof of concept script for comparative genomics with OrthoFinder')
parser.add_argument('--orthofinder_output', required=True, nargs=1,metavar="DIR",help="OrthoFinder output directory")

args = parser.parse_args()

class Orthogroup:
    def __init__(self,line,species):
        tab_split = line.split('\t')
        orthogroup = tab_split[0].strip()
        
        self.species = species ##Have to store the species in relation to the orthogroup, as the orthogroup is only informative with respect to species it was calculated against.
        ##Yes it isn't very memory efficient, but oh well
        ##Alternatively could implement an "Orthogroup holder" class which keeps track of the species

        ##Implement a list of lists to hold the per species genes.
        ##Index of first list is the species
        ##Index of second list is the gene
        self.per_species_gene_ids = []
        per_line_gene_ids = tab_split[1:].strip() ##strip removes the newline at the end of the line
        
        for k in range(0,len(per_line_gene_ids)):
            genes = per_line_gene_ids[k].split(",")
            for z in range(0,len(genes)):
                genes[z] = genes[z].strip() ##Remove the prefix / suffix whitespace
            s = species[k]
            if len(genes) == 1 and genes[0] == "":
                genes = [] ##Change it to an empty list
            self.per_species_gene_ids.append(genes)
               
        self.name = orthogroup
        assert len(self.species) == len(self.per_species_gene_ids)

    def get_gene_ids_via_species_index(self,i):
        return self.per_species_gene_ids[i]

    def get_gene_ids_via_species_name(self,name):
        i = self.species.index(name)
        return self.get_gene_ids_via_species_index(i)

    def get_all_gene_ids(self):
        the_genes = []
        for i in range(0,len(self.per_species_gene_ids)):
            the_genes += self.get_gene_ids_via_species_index(i)
        return the_genes        

    def is_reciprocally_direct_orthogroup(self):
         for i in range(0,len(self.per_species_gene_ids)):
            if len(self.get_gene_ids_via_species_index(i)) != 1:
                return False
         return True
            
    def __str__(self):
        return ",".join(self.get_all_gene_ids())

    def __print__(self):
        return {self.name:self.get_all_gene_ids}

##Load species.
##Can also get this from the header of the Orthogroups.csv file
##Note that the species ID within OrthoFinder has the .fasta or .fa removed
path = args.orthofinder_output[0]+"WorkingDirectory/SpeciesIDs.txt"
species = []
with open(path) as f:
    species = f.read().splitlines()

for i in range(0,len(species)):
    species[i] = species[i].split(" ")[1]
print("Loaded",len(species),"species from OrthoFinder results")

##Load sequence IDs
sequence_name_to_species_dict = dict()
species_to_sequences = dict()
path = args.orthofinder_output[0]+"WorkingDirectory/SequenceIDs.txt"
sequence_ids = []
with open(path) as f:
    sequence_ids = f.read().splitlines()
for i in range(0,len(sequence_ids)):
    splitline = sequence_ids[i].split("_")
    species_index = int(splitline[0])
    s = species[species_index]
    sequence_index = int(splitline[1].split(":")[0])
    sequence_name = ":".join("_".join(splitline[1:]).split(":")[1:]).strip()
    sequence_name_to_species_dict[sequence_name] = s
    if s not in species_to_sequences.keys():
        species_to_sequences[s] = [sequence_name]
    else:
        species_to_sequences[s].append(sequence_name)
print("done")

##Finding the Orthologues folder
Orthologues_dir_path = glob.glob(args.orthofinder_output[0]+"WorkingDirectory/Orthologues_*")
if len(Orthologues_dir_path) == 0:
    print("No Orthologues folder found. Are you sure you ran OrthoFinder completely?")
elif len(Orthologues_dir_path) > 1:
    print("Multiple Orthologues folder found. Only one can be used. Exiting.")
    exit()

####Load the Orthogroups.csv file into a datastructure
####Structure of the Orthogroups.csv file:
####Different species separated by tabs, genes within a given species separated by commas
print("Starting to load the Orthogroups.csv file...")
orthogroups = dict()
path = args.orthofinder_output[0]+"WorkingDirectory/Orthogroups.csv"
handle = open(path,"rU")
i=0
for line in handle.readlines():
    if i == 0: ##Header line
        i+=1
        continue
    og = Orthogroup(line,species)
    orthogroups[og.name] = og
    ##print(og)
    i+=1

print("Loaded",len(orthogroups),"orthogroups from the Orthogroups.csv file.")
print("OG0014256",orthogroups["OG0014256"].is_reciprocally_direct_orthogroup())
print("OG0000001",orthogroups["OG0000001"].is_reciprocally_direct_orthogroup())

##Finding 1-1-1 orthogroups
