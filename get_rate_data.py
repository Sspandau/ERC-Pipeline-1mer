#! /usr/bin/env python3
​
import os.path as osp
import os
import re
import sys
​
from statsmodels.stats.multitest import multipletests
​
from pipeline import register_previous_run, generate_network, _30mya_cutoff, _20mya_cutoff, find_tree, get_rates
from utilities import safe_mkdir, rho_sorted_neighbors, safe_phylo_read, mammal_taxa_info
​
​
def main():
    args = sys.argv[1:]
​
    assert len(args) >= 4, "Arguments: python3 get_rate_data.py full/20my/30my pipeline_directory output protein_id_dir [protein_id_dir 2]"
​
    tree = safe_phylo_read(osp.join("data", "finished_mam_timetree.nwk"))
    taxon_set = args[0].lower()
​
    assert taxon_set in ("full", "20my", "30my"), "You must specify the taxon set to get rates for."
​
    taxa_list = None
    if taxon_set == '20my':
        taxa_list = _20mya_cutoff
    elif taxon_set == '30my':
        taxa_list = _30mya_cutoff
​
    pipeline_dir = args[1]
    assert osp.exists(pipeline_dir), f"Pipeline directory {pipeline_dir} doesn't exist!"
    register_previous_run(pipeline_dir)
​
    output = args[2]
    assert output.endswith(".csv"), f"The output ({output}) must be a .csv file!"
​
    protein1_dir = args[3]
    protein_1_files = {}
​
    for filename in os.listdir(protein1_dir):
        positions = filename.strip('_')
        aa_pos = int(positions[4])
        protein_1_files[aa_pos] = filename #create list of filenames for each aa position
        protein_1_files = sorted(protein_1_files)
        #lsorted = sorted(protein_1_files, key=lambda x: int(os.path.splitext(x)[0]))
​
    if len(args) > 4:
        protein2_dir = args[4]
    else:
        protein2_dir = None
    
    if protein2_dir:
        protein_2_files = {}
        for filename in os.listdir(protein2_dir):
            positions = filename.strip('_')
            aa_pos = int(positions[4])
            protein_2_files[aa_pos] = filename
    else:
        protein_2_files = None
        
    with open(output, 'w') as f:
        f.write("protein1,protein2,taxon,taxon order,rate1,rate2,time\n")
        
    my_file.write("protein,taxon,taxon order,time,rate\n")
        
    for pos in protein_1_files: #for loop for first protein's aa positions
        print(protein, protein_1_files[pos])
        protein = find_tree(protein_1_files[protein])
        if protein_2_files: #for loop for second protein's aa positions
            for protein2 in protein_2_files:
                protein2 = find_tree(protein_2_files[protein2])
                name2taxa = mammal_taxa_info(name_as_key=True)
                print(f"Getting rate data for {args[3]} with {args[4]}...")
                taxa, rates = get_rates(tree, True, taxa_list, protein, protein2)
                for (taxon, rate1, rate2, time) in zip(taxa, rates[1], rates[2], rates[0]):
                    order = name2taxa[taxon.strip()].order
                    if protein == protein_1_list[0] and protein2 == protein_2_list[0]:
                        f.write(f"{args[3]},{args[4]},{taxon},{order},{rate1},{rate2},{time}\n") #write names of taxons, orders in first few columns
                    else:
                        f.write(f"{rate1},{rate2}\n") #append rates after adding first aa positions' rates
        else: #no for loop for second protein if no second protein
            name2taxa = mammal_taxa_info(name_as_key=True)
            print(f"Getting rate data for {args[3]}...")
            taxa, rates = get_rates(tree, True, taxa_list, protein)
            #with open(output, 'w') as f:
            my_file = open(output, 'a')
            # lines = my_file.readlines()
​
            # Create a list to store the data
            data = []
​
            # Loop through the data and store it in the list
            for (taxon, rate1, time) in zip(taxa, rates[1], rates[0]):
                order = name2taxa[taxon.strip()].order
                if protein == protein_1_files[0]:
                    data.append([args[3], taxon, order, time, rate1])
                else:
                    data.append([rate1])
​
            # Write the header to the file
​
            # Write the data to the file, column by column
            for i in range(len(data[0])):
                column = [row[i] for row in data]
                my_file.write(','.join(str(x) for x in column) + '\n')
​
            # Close the file
            my_file.close()
​
if __name__ == "__main__":
    main()
