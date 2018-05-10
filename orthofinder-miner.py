import Bio
import argparse

parser = argparse.ArgumentParser(description='Proof of concept script for comparative genomics with OrthoFinder')
parser.add_argument('--orthofinder_output', required=True, nargs=1,metavar="DIR",help="OrthoFinder output directory")

args = parser.parse_args()

##Load species.
species = list(open(args.orthofinder_output+"WorkingDirectory/SpeciesIDs.txt","rU").readlines())
print(species)

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
