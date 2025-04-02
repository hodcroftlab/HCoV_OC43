#To run this script use: 
# python extract_gene_from_whole_genome.py --genbank_file ../ingest/data/references/nl63_full_reference.gb --output_directory ../data/references --product_names "spike protein" "membrane protein" "nucleocapsid protein"
#use  nargs='+', if want multiple protein names
import os
import argparse
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--genbank_file', required=True, help='Input reference file in genbank format (.gb)')
    parser.add_argument('--output_fasta', required=True, help='Output directory for FASTA files')
    parser.add_argument('--product_name',  required=True, help='Product names of proteins to extract - ensure its the same as CDS names in reference.gb file')
    
    return parser.parse_args()


def extract_protein(genbank_file, output_fasta, product_name):
    os.makedirs(os.path.dirname(output_fasta), exist_ok=True)  
    
    if product_name.lower() == "whole_genome":
        for record in SeqIO.parse(genbank_file, "genbank"):
            with open(output_fasta, "w") as output_handle:
                output_handle.write(f">{record.id}\n{record.seq}\n")
            print(f"Whole genome sequence saved to {output_fasta}")
        return

    found = False
    for record in SeqIO.parse(genbank_file, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                product_list = feature.qualifiers.get("product", [])
                if any(product_name.lower() in p.lower() for p in product_list):
                    nucleotide_sequence = feature.location.extract(record.seq)

                    with open(output_fasta, "w") as output_handle:
                        output_handle.write(f">{product_name}\n{nucleotide_sequence}\n")

                    print(f"Protein {product_name} saved to {output_fasta}")
                    found = True
                    return  
    
    if not found:
        print(f"Protein {product_name} not found in the GenBank file.")

if __name__ == "__main__":
    args = parse_args()
    extract_protein(args.genbank_file, args.output_fasta, args.product_name)


# Input GenBank file
#genbank_file = "../ingest/data/references/nl63_full_reference.gb"

# Output FASTA file

#output_directory = "../data/references" #make sure you've created the directory first

# Enter the name of the proteins you want to extract - ensure it's the same as what's in the CDS of the reference.gb file 
# product_names = [
#     "spike protein",
#     "envelope protein",
#     "membrane protein",
#     "nucleocapsid protein"
# ] 

# def extract_proteins(genbank_file, output_fasta, product_names):
    
# # Parse the GenBank file
#     for record in SeqIO.parse(genbank_file, "genbank"):
#         for product_name in product_names:  
#             found = False
#             for feature in record.features:
#                 #if feature.type == "CDS" and product_name in feature.qualifiers.get("product", []): #use exact same name as in CDS
#                 if feature.type == "CDS":
#                     product_list = feature.qualifiers.get("product", [])
#                     if any(product_name.lower() in p.lower() for p in product_list):    #using a partial match instead to be a bit more flexible
#                         # Extract the nucleotide sequence using the feature's location 
#                         nucleotide_sequence = feature.location.extract(record.seq)

#                         # Create output file and replace spaces with underscores
#                         output_file = os.path.join(output_fasta, f"{product_name.replace(' ', '_')}.fasta")
                        

#                         # Write to FASTA format
#                         with open(output_fasta, "w") as output_handle:
#                             output_handle.write(f">{product_name}\n{nucleotide_sequence}\n")
#                         print(f"Protein {product_name} saved to {output_fasta}")
#                         found = True
#                         break
#             if not found:
#                 print(f"Protein {product_name} not found in the GenBank file.")
            
# if __name__ == "__main__":
#     args = parse_args()
    
# os.makedirs(args.output_directory, exist_ok=True)
                
# extract_proteins(args.genbank_file, args.output_directory, args.product_names)