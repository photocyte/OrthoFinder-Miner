import Bio
import argparse

parser = argparse.ArgumentParser(description='Proof of concept script for comparative genomics with OrthoFinder')
parser.add_argument('--orthofinder_output', required=True, nargs=1,metavar="DIR",help="OrthoFinder output directory")

args = parser.parse_args()

##Load species.
path = args.orthofinder_output[0]+"WorkingDirectory/SpeciesIDs.txt"
species = []
with open(path) as f:
    species = f.read().splitlines()

for i in range(0,len(species)):
    species[i] = species[i].split(" ")[1]
print(species)

##Load sequence IDs
path = args.orthofinder_output[0]+"WorkingDirectory/SequenceIDs.txt"
sequence_ids = []
with open(path) as f:
    sequence_ids = f.read().splitlines()
for i in range(0,len(sequence_ids)):
    splitline = sequence_ids[i].split("_")
    species_index = splitline[0]
    sequence_index = splitline[1].split(":")[0]
    sequence_name = ":".join("_".join(splitline[1:]).split(":")[1:]).strip()
    print(species_index,sequence_index,sequence_name)
exit()

####Load the Orthogroups file into a datastructure
print("Starting to load the Orthogroups.csv file...")
orthogroup_dict = dict()
handle = open(args.orthofinder_output[0],"rU")
for line in handle.readlines():
    splitline = line.split('\t')
    orthogroup = splitline[0].strip()
    if orthogroup == "OG0014256" or orthogroup == "OG0000001":
        print splitline

    protein_ids = splitline[1:]
    while '' in protein_ids:
        protein_ids.remove('')
    orthogroup_dict[orthogroup] = protein_ids


print("Finished loading the Orthogroups.csv file.")
    #print orthogroup_dict
