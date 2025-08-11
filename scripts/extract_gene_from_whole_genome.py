#To run this script use: 
# python extract_gene_from_whole_genome.py --genbank_file ../ingest/data/references/nl63_full_reference.gb --output_directory ../data/references --product_names "spike protein" "membrane protein" "nucleocapsid protein"
#use  nargs='+', if want multiple protein names
import os
import argparse
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--genbank_file', required=True, help='Input reference file in genbank format (.gb)')
    parser.add_argument('--output_fasta', required=True, help='Output directory for FASTA files')
    parser.add_argument('--product_name',  required=True, help='Product names of proteins to extract - ensure its the same as CDS names in reference.gb file')
    parser.add_argument('--output_genbank', required=True, help='Output genbank file')
    
    return parser.parse_args()


def extract_protein(genbank_file, output_fasta, product_name, output_genbank):
    os.makedirs(os.path.dirname(output_fasta), exist_ok=True)  
    
    if product_name.lower() == "whole-genome":
        for record in SeqIO.parse(genbank_file, "genbank"):
            with open(output_fasta, "w") as output_handle:
                output_handle.write(f">{record.id}\n{record.seq}\n")
            print(f"Whole genome sequence saved to {output_fasta}")
            
            with open(output_genbank, "w") as output_handle:
                SeqIO.write(record, output_handle, "genbank")
            print(f"Whole genome GenBank file saved to {output_genbank}")
        return

    found = False
    for record in SeqIO.parse(genbank_file, "genbank"):
        
        locus = record.annotations.get('locus', 'Unknown_Locus')
        date = record.annotations.get('date', 'Unknown_Date')
        accession = record.id  
        version = record.annotations.get('version', 'Unknown_Version')
        source = record.annotations.get('source', 'Unknown_Source')
        organism = record.annotations.get('organism', 'Unknown_Organism')
                    
        for feature in record.features:
            if feature.type == "CDS":
                product_list = feature.qualifiers.get("product", [])
                if any(product_name.lower() in p.lower() for p in product_list):
                    nucleotide_sequence = feature.location.extract(record.seq)

                    with open(output_fasta, "w") as output_handle:
                        output_handle.write(f">{product_name}\n{nucleotide_sequence}\n")
                    print(f"Protein {product_name} saved to {output_fasta}")
                    
                   
                    extracted_record = record[feature.location.start:feature.location.end]
                    extracted_record.id = f"{product_name}_extracted_{accession}"
                    extracted_record.description = f"{product_name} gene sequence ({organism}, {source}, {date}, Accession: {accession}, Version: {version})"

                    extracted_record.annotations['locus'] = locus
                    extracted_record.annotations['date'] = date
                    extracted_record.annotations['accession'] = accession
                    extracted_record.annotations['version'] = version
                    extracted_record.annotations['source'] = source
                    extracted_record.annotations['organism'] = organism
                    
                    # Find the original source feature
                    original_source = next((f for f in record.features if f.type == "source"), None)

                    # Use its qualifiers if available, or set default
                    if original_source:
                        new_qualifiers = original_source.qualifiers.copy()
                    else:
                        new_qualifiers = {
                            "organism": [organism],
                            "mol_type": ["genomic RNA"],
                            "country": ["Unknown"],
                            "isolate": ["Unknown"],
                            "strain": ["Unknown"],
                            "source": [source]
                        }

                    # Create a new source feature for the extracted region
                    source_feature = SeqFeature(
                        FeatureLocation(0, len(extracted_record.seq)), 
                        type="source",
                        qualifiers=new_qualifiers
                    )


                    extracted_record.features.append(source_feature)

                    with open(output_genbank, "w") as output_handle:
                        SeqIO.write(extracted_record, output_handle, "genbank")
                    print(f"Extracted gene GenBank file saved to {output_genbank}")
                    
                    found = True
                    return  
    
    if not found:
        print(f"Protein {product_name} not found in the GenBank file.")

if __name__ == "__main__":
    args = parse_args()
    extract_protein(args.genbank_file, args.output_fasta, args.product_name, args.output_genbank)


