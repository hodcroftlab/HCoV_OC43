# Snakemake execution:

# To run a default protein xy run:
# snakemake  auspice/HCoV_OC43_protein_xy.json --cores 9

# To run a default whole genome run (>6400bp):
# snakemake auspice/HCoV_OC43_genome.json --cores 9

# Excluding envelope for OC43 because not well annotated 
###############
wildcard_constraints:
    seg="spike|nucleocapsid|membrane|whole_genome"  

# Define segments to analyze
segments = ["spike", "nucleocapsid",  "membrane", "whole_genome"] # This is only for the expand in rule all

# Expand augur JSON paths
rule all:
    input:
        #augur_jsons = expand("auspice/HCoV_OC43_{segs}.json", segs=segments) ## TODO: replace <your_virus> with actual virus name (Ctrl+H)
        augur_jsons = expand("auspice/HCoV_OC43_{segs}-accession.json", segs=segments) ## TODO: replace <your_virus> with actual virus name (Ctrl+H)

##############################
# Rule to handle input and config files
###############################

rule files:
    input:
        sequence_length =   "{seg}",
        dropped_strains =   "config/dropped_strains.txt",
        reference =         "ingest/data/references/oc43_full_reference.gb"    ,
        lat_longs =         "config/lat_longs.tsv",
        auspice_config =    "{seg}/config/auspice_config.json",
        colors =            "config/colors.tsv",
        clades =            "{seg}/config/clades_genome.tsv",
        regions=            "config/geo_regions.tsv",
        metadata=           "data/metadata.tsv",
        extended_metafile=  "data/meta_manual.tsv",  ###TODO: Add an empty tsv file to this path or metadata for one of your sequences

files = rules.files.input

##############################
# Download from NBCI Virus with ingest snakefile
###############################

# workdir: "ingest"
# include: "ingest/Snakefile"
# workdir: ".."


rule fetch:
    input:
        dir = "ingest"
    output:
        sequences="data/sequences.fasta",
        metadata=files.metadata
    params:
        seq="ingest/data/sequences.fasta",
        meta="ingest/data/metadata.tsv"
    shell:
        """
        cd {input.dir} 
        snakemake --cores 9 all
        cd ../
        cp -u {params.seq} {output.sequences}
        cp -u {params.meta} {output.metadata}
        """


###################################################

##############################
# Update strain names
# If strain name == accession -> fetching real strain names from genbank
# Depending on how many sequences you have, it will run for a long time! >30min. Comment out to skip!
###############################

rule update_strain_names:
    message:
        """
        Updating strain information in metadata.
        """
    input:
        file_in =  files.metadata
    params:
        backup = "data/strain_names_previous_run.tsv" 
    output:
        file_out = "data/updated_strain_names.tsv"
    shell:
        """
        time bash scripts/update_strain.sh {input.file_in} {params.backup} {output.file_out}
        """



##############################
# BLAST
# blast fasta files for HCoV-OC463 proteins
# Extract protein-specific or whole-genome sequences
###############################

rule extract:
    input: 
        genbank_file = files.reference
    output: 
        extracted_fasta = "{seg}/results/extracted.fasta",    
        extracted_genbank = "{seg}/results/extracted.gbk" 
    params:
        product_name = "{seg}"
    shell:
        """
        python scripts/extract_gene_from_whole_genome.py \
        --genbank_file {input.genbank_file} \
        --output_fasta {output.extracted_fasta} \
        --product_name {params.product_name} \
        --output_genbank {output.extracted_genbank}

        """

rule blast:
    input: 
        blast_db_file = rules.extract.output.extracted_fasta,  # Provide a BLAST reference
        seqs_to_blast = rules.fetch.output.sequences
    output:
        blast_out = "{seg}/results/blast_out.csv"
    params:
        blast_db = "{seg}/results/blast_database"
    shell:
        """
        sed -i 's/-//g' {input.seqs_to_blast}
        makeblastdb -in {input.blast_db_file} -out {params.blast_db} -dbtype nucl
        blastn -task blastn -query {input.seqs_to_blast} -db {params.blast_db} \
            -outfmt '10 qseqid sseqid pident length mismatch gapopen qstart qend sstart \
            send evalue bitscore qcovs' \
            -out {output.blast_out} -evalue 0.0005
        """


rule blast_sort:
    input:
        blast_result = rules.blast.output.blast_out,  # BLAST output for specific proteins
        input_seqs = rules.fetch.output.sequences
    output:
        sequences = "{seg}/results/sequences.fasta",
        blast_length= "{seg}/results/blast_{seg}_length.tsv"
    params:
        range = "{seg}",  # Determines which protein (or whole genome) is processed
        min_length = lambda wildcards: {"spike": 2472, "nucleocapsid": 810, "membrane": 450, "whole_genome": 20000}[wildcards.seg],  # Min length 
        max_length = lambda wildcards: {"spike": 4120, "nucleocapsid": 1450, "membrane": 750, "whole_genome": 30738}[wildcards.seg]  # Max length added 100 to actual length
    shell:
        """
        python scripts/blast_sort.py --blast {input.blast_result} \
            --seqs {input.input_seqs} \
            --out_seqs {output.sequences} \
            --out_length {output.blast_length} \
            --range {params.range} \
            --min_length {params.min_length} \
            --max_length {params.max_length}

        
        """
# snakemake -c 9 "temp/nucleocapsid/blast_out.csv" -f
# snakemake -c 9 nucleocapsid/results/sequences.fasta

##############################
# AUGUR CURATE AND MERGE
# Change the format of the dates in the metadata
# Attention: ```augur curate``` only accepts iso 8 formats; please make sure that you save e.g. Excel files in the correct format
# Merge with other metadata files you might have

###############################

# rule curate_meta_dates:
#     message:
#         """
#         Cleaning up metadata with augur curate and merge your metadata with the one from ingest
#         """
#     input:
#         meta = files.metadata,
#         blast_length = rules.blast_sort.output.blast_length, 
#     params:
#         strain_id_field = "accession",
#         sorted_meta_accessions = "temp/{seg}/meta_accessions.txt",
#         sorted_blast_length_accessions = "temp/{seg}/blast_accessions.txt",
#         common_accessions = "temp/{seg}/common_accessions.txt",
    
#     output:
#         final_metadata = "{seg}/results/metadata.tsv",
#         filtered_meta = "{seg}/results/filtered_metadata.tsv",
#         filtered_blast_length = "{seg}/results/filtered_blast_length.tsv"
#     shell:
#         """
    
#         #Preventing Augur merge sqlite right/full join error 
#         cut -f1 {input.meta} | sort > {params.sorted_meta_accessions}
#         cut -f1 {input.blast_length} | sort > {params.sorted_blast_length_accessions}
#         comm -12 {params.sorted_meta_accessions} {params.sorted_blast_length_accessions} > {params.common_accessions}
#         (grep -Ff {params.common_accessions} {input.meta}) > {output.filtered_meta}
#         ( grep -Ff {params.common_accessions} {input.blast_length}) > {output.filtered_blast_length}
#         sort -k1,1 {output.filtered_meta} > {output.filtered_meta}
#         sort -k1,1 {output.filtered_blast_length} > {output.filtered_blast_length}



#         # Merge curated metadata
#         # augur merge --metadata meta={output.filtered_meta} extended_meta={output.filtered_blast_length} \
#         #     --metadata-id-columns {params.strain_id_field} \
#         #     --metadata-delimiters '\t' \
#         #     --output-metadata {output.final_metadata}

#         paste {output.filtered_meta} {output.filtered_blast_length} > {output.final_metadata}
#         """
###########################################################

rule index_sequences:
    message:
        """
        Creating an index of sequence composition for filtering.
        """
    input:
        sequences = rules.blast_sort.output.sequences
    output:
        sequence_index = "{seg}/results/sequence_index.tsv"
    #conda: "ingest/workflow/envs/nextstrain.yaml"
    shell:
        """
        augur index \
            --sequences {input.sequences} \
            --output {output.sequence_index}
        """

rule filter:
    message:
    #- {params.sequences_per_group} sequence(s) per {params.group_by!s}
        """
        Filtering to
          
          - from {params.min_date} onwards
          - excluding strains in {input.exclude}
        """
    input:
        sequences = rules.blast_sort.output.sequences,
        sequence_index = rules.index_sequences.output.sequence_index,
        # metadata = rules.curate_meta_dates.output.final_metadata,
        metadata = files.metadata,
        exclude = files.dropped_strains
    output:
        sequences = "{seg}/results/filtered.fasta",
        metadata = "{seg}/results/filtered_metadata.tsv"
    params:
        #group_by = "country year month", --group-by {params.group_by} \
        #sequences_per_group = 20,  --sequences-per-group {params.sequences_per_group} \
        min_date = 1960, #Set to 1960, as this was when HCoV first discovered, but NL63 2004
        strain_id_field= "accession",
        #min_length = 20000 #--min-length {params.min_length}
        
    #conda: "ingest/workflow/envs/nextstrain.yaml"
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --sequence-index {input.sequence_index} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id_field} \
            --exclude {input.exclude} \
            --exclude-ambiguous-dates-by year \
            --output {output.sequences} \
            --output-metadata {output.metadata} \
            --min-date {params.min_date} \
            
        """

###########################
#Nextclade needs the reference in fasta format, whereas Augur prefers the .gb format
###########################

rule align:
    message:
            """
            Aligning sequences to {input.reference} using Nextalign.
            """
    input:
        sequences = rules.filter.output.sequences,
        reference = rules.extract.output.extracted_fasta
    output:
        alignment = "{seg}/results/aligned.fasta"

    params:
            nuc_mismatch_all = 10,
            nuc_seed_length = 30
    shell:
        """
        nextclade run \
        {input.sequences}  \
        --input-ref {input.reference}\
        --allowed-mismatches {params.nuc_mismatch_all} \
        --min-length {params.nuc_seed_length} \
        --include-reference false \
        --retry-reverse-complement true \
        --output-fasta {output.alignment} 
        """

##############################
# Building a tree
###############################

rule tree:
    message:
        """
        Creating a maximum likelihood tree
        """
    input:
        alignment = rules.align.output.alignment

    output:
        tree = "{seg}/results/tree_raw.nwk"

    threads: 9
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --nthreads {threads}\
            --output {output.tree}
        """

# ##############################
# # Refine to a timeline
# ###############################

rule refine:
    message:
        """
        Refining tree by rerooting and resolving polytomies
          - estimate timetree
          - use {params.coalescent} coalescent timescale
          - estimate {params.date_inference} node dates
          - filter tips more than {params.clock_filter_iqd} IQDs from clock expectation
        """
    input:
        tree = rules.tree.output.tree,
        alignment = rules.align.output.alignment,
        metadata = rules.filter.output.metadata
    output:
        tree = "{seg}/results/tree.nwk",
        node_data = "{seg}/results/branch_lengths.json"
    params:
        coalescent = "opt",
        date_inference = "marginal",
        clock_filter_iqd = 3, # set to 6 if you want more control over outliers
        strain_id_field ="accession",
        # clock_rate = 0.004, # remove for estimation by augur; check literature
        # clock_std_dev = 0.0015

    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id_field} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --timetree \
            --coalescent {params.coalescent} \
            --date-confidence \
            --date-inference {params.date_inference} \
            --clock-filter-iqd {params.clock_filter_iqd}
        """


            # --clock-rate {params.clock_rate}\
            # --clock-std-dev {params.clock_std_dev} \
# ##############################
# # Ancestral sequences and amino acids
# ###############################

rule ancestral:
    message: "Reconstructing ancestral sequences and mutations"
    input:
        tree = rules.refine.output.tree,
        alignment = rules.align.output.alignment

    output:
        node_data = "{seg}/results/nt_muts.json"

    params:
        inference = "joint"
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --keep-ambiguous\
            --inference {params.inference}
        """
 
rule translate:
    message: "Translating amino acid sequences"
    input:
        tree = rules.refine.output.tree,
        node_data = rules.ancestral.output.node_data,
        reference = rules.extract.output.extracted_genbank
    output:
        node_data = "{seg}/results/aa_muts.json"

    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output-node-data {output.node_data}
        """

rule traits:
    message: "Inferring ancestral traits for {params.traits!s}"
    input:
        tree = rules.refine.output.tree,
        metadata = rules.filter.output.metadata
    output:
        node_data = "{seg}/results/traits.json"
        
    params:
        traits = "country",
        strain_id_field= "accession"
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id_field} \
            --output-node-data {output.node_data} \
            --columns {params.traits} \
            --confidence
        """

# ##############################
# # Assign clades or subgenotypes based on list provided
# ###############################
# rule clades: 
#     message: "Assigning clades according to nucleotide mutations"
#     input:
#         tree=rules.refine.output.tree,
#         aa_muts = rules.translate.output.node_data,
#         nuc_muts = rules.ancestral.output.node_data,
#         clades = files.clades # TODO: assign mutations to specific clades
#     output:
#         clade_data = "{seg}/results/clades.json"

#     shell:
#         """
#         augur clades --tree {input.tree} \
#             --mutations {input.nuc_muts} {input.aa_muts} \
#             --clades {input.clades} \
#             --output-node-data {output.clade_data}
#         """

# #########################
# #  EXPORT
# #########################
rule export:
    message: "Creating auspice JSONs"
    input:
        tree = rules.refine.output.tree,
        # metadata = rules.curate_meta_dates.output.final_metadata,
        metadata = files.metadata,
        branch_lengths = rules.refine.output.node_data,
        traits = rules.traits.output.node_data,
        nt_muts = rules.ancestral.output.node_data,
        aa_muts = rules.translate.output.node_data,
        # clades = rules.clades.output.clade_data,
        colors = files.colors,
        lat_longs = files.lat_longs,
        auspice_config = files.auspice_config
    params:
        strain_id_field= "accession"

    output:
        auspice_json = "auspice/HCoV_OC43_{seg}-accession.json"
        
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id_field} \
            --node-data {input.branch_lengths} {input.traits} {input.nt_muts} \
                {input.aa_muts} \
            --colors {input.colors} \
            --lat-longs {input.lat_longs} \
            --auspice-config {input.auspice_config} \
            --output {output.auspice_json}
        """
        
# ###############################
# # Change from accession to strain name view in tree
# ################################

# rule rename_json:
#     input:
#         auspice_json= rules.export.output.auspice_json,
#         metadata =files.metadata,
#     output:
#         auspice_json="auspice/HCoV_OC43_{seg}.json"
#     params:
#         strain_id_field="accession",
#         display_strain_field= "strain"
#     shell:
#         """
#         python3 scripts/set_final_strain_name.py --metadata {input.metadata} \
#                 --metadata-id-columns {params.strain_id_field} \
#                 --input-auspice-json {input.auspice_json} \
#                 --display-strain-name {params.display_strain_field} \
#                 --output {output.auspice_json}

#         mkdir -p auspice/accession/ && mv {input.auspice_json} auspice/accession/
#         """

rule clean:
    message: "Removing directories: {params}"
    params:
        "*/results/*",
        "auspice/*",
        "temp/*", 
        "data/metadata.tsv",
        "data/sequences.fasta",
        "ingest/data/metadata.tsv",
        "ingest/sequences.fasta"

    shell:
        "rm -rfv {params}"