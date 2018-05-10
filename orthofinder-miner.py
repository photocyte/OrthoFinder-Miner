import Bio
import argparse

parser = argparse.ArgumentParser(description='Proof of concept script for comparative genomics with OrthoFinder')
parser.add_argument('--orthofinder_output', required=True, nargs=1,metavar="(FILE)",help="OrthoFinder orthogroups.csv file")

args = parser.parse_args()


