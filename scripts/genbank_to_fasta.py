"""
This script converts the reference sequences which are in genbank 
format into fasta files that can be used by nextclade.

"""
import glob
from Bio import SeqIO

directory = "ingest/data/references/"

for genbankfile in glob.glob(f"{directory}/*.gb"):
    
    outputfile = genbankfile.replace(".gb", ".fasta")
    
    count = SeqIO.convert(genbankfile, "genbank",outputfile, "fasta")