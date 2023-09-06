from Bio.Seq import Seq
import json
import argparse
import csv


parser = argparse.ArgumentParser()

parser.add_argument('--input_nt_file', type=str, help='path to alignment nt-muts json file')
parser.add_argument('--input_escape_summary_file', type=str, help='path to DMS csv file')
parser.add_argument('--output_total_escape_file', type=str, help='path to output json file with total escapes')
parser.add_argument('--output_mutation_summary_file', type=str, help='path to output json file with summary of escape mutations')
parser.add_argument('--output_mutation_list_file', type=str, help='path to output json file with mutations listed as string for tree')
parser.add_argument('--exclude_cleavage_site', type=bool, default=False, help='option to exclude mutations in the cleave site (defaults to False)')

args = parser.parse_args()

input_nt_file = args.input_nt_file
input_escape_summary_file = args.input_escape_summary_file
output_total_escape_file = args.output_total_escape_file
output_mutation_summary_file = args.output_mutation_summary_file
output_mutation_list_file = args.output_mutation_list_file
exclude_cleavage_site = args.exclude_cleavage_site


#alignment = SeqIO.parse(open(alignment_file),'fasta')

# CDS of gsgd reference sequence spans from nt 21 to nt 1727, so use 21 and 1728 as start and end of range
cds_start = 21
cds_end = 1728

# gsgd reference sequence (which all seqs in tree were aligned to) contains an extra aa in the polybasic cleavage
# site with respect to the Indiana/2022 DMS strain. thus, the DMS site numbering is n-1 compared to the alignment
# numbering after this gap. to fix, can just exclude this codon when pulling out each CDS from the alignment.
# this gap spans from nt 1053 to nt 1055, so use 1053 and 1056 as start and end of range
gap_start = 1053
gap_end = 1056


aa_seqs = {}

with open(input_nt_file) as input_nt:
    nt_file = json.load(input_nt)
    for node, value in nt_file['nodes'].items():
        cds_str = value['sequence']
        cds = Seq(cds_str)
        cds = cds[cds_start:gap_start] + cds[gap_end:cds_end]
        cds = cds.replace('-','n')
        aa = cds.translate()
        aa_seqs[node] = str(aa)


# for seq in alignment:
#     name = seq.name
#     # get the CDS, minus the polybasic cleavage site gap
#     cds = seq.seq[cds_start:gap_start] + seq.seq[gap_end:cds_end]
#     # replace '-' with 'n' for bioseq translation (codons with 'n' translated as 'X', codons with '-' raise error)
#     cds = cds.replace("-","n")
#     # translate the CDS and save the translation as a str in the dict
#     aa = cds.translate()
#     aa_seqs[name] = str(aa)

# convert the escape summary CSV into a nested dictionary
# this makes it much faster to pull the escape values than filtering the entire df each time
escape_dict = {}

with open(input_escape_summary_file, 'r') as escape_summary:
    # read escape_summary csv
    reader = csv.reader(escape_summary)
    # pull out header names
    headers = next(reader)
    # find indices of columns needed -- if header order changes, will not need to manually edit
    # however, headers 'site', 'mutant', and 'escape' must maintain naming
    [site_index, mutant_index, escape_index] = [headers.index(name) for name in ['site', 'mutant', 'escape']]
    # iterate through rows of data
    for row in reader:
        site = int(row[site_index])
        mutant = row[mutant_index]
        # if escape is blank, set it to 0; otherwise, read in as a float
        escape = float(row[escape_index]) if row[escape_index] != '' else 0
        # if escape is negative (or 0 / was blank) we do not need to include in the dictionary
        # only want to sum escape values that are positive (per Jesse)
        if escape > 0:
            # if site already in escape_dict, add to its nested dict
            if site in escape_dict:
                escape_dict[site][mutant] = escape
            # otherwise, add to escape_dict and create nested dict
            else:
                escape_dict[site] = {mutant: escape}


# exclude GsGd seq from analysis, not on tree
exclude_from_analysis = ['AF144305.1']

cleavage_site_aas = range(339, 346)

escape_values = {}
mutation_summary_dict = {}

# iterate through all non-excluded aa_seqs and determine their total_escape
for name, aa_seq in [(name, aa_seq) for name, aa_seq in aa_seqs.items() if name not in exclude_from_analysis]:
    total_escape = 0
    site = 0
    mutations = {}
    for aa in aa_seq:
        site = site + 1
        # if exclude_cleavage_site set to True, and the site is in the cleavage site, skip; otherwise, try adding that site's escape to total_escape
        if not (exclude_cleavage_site and site in cleavage_site_aas):
            # if the mutation is found in escape_dict, add it to the strain's total_escape
            try:
                total_escape += escape_dict[site][aa]
                mutations[str(site)+aa] = escape_dict[site][aa]
            # if it's not in escape_dict, pass -- value either missing or was negative
            except KeyError:
                pass
    # save total_escape in the dict after rounding
    escape_values[name] = {'escape': round(total_escape,5)}
    mutation_summary_dict[name] = mutations


mutation_list_dict = {}


for node, mut_esc_dict in mutation_summary_dict.items():
    outlist = []
    for mut, esc in sorted(mut_esc_dict.items(), key=lambda x:x[1], reverse=True):
        outlist.append(f'{mut} [{esc}]')
    outstr = ', '.join(outlist)
    mutation_list_dict[node] = {'escape mutations': outstr}





# format dict for auspice
escape_json_dict = {'nodes': escape_values}
mutation_summary_json_dict = {'nodes': mutation_summary_dict}
mutation_list_json_dict = {'nodes': mutation_list_dict}


# save dict as json file
escape_json_object = json.dumps(escape_json_dict)
with open(output_total_escape_file, 'w') as output:
    output.write(escape_json_object)

mutation_summary_json_object = json.dumps(mutation_summary_json_dict)
with open(output_mutation_summary_file, 'w') as output:
    output.write(mutation_summary_json_object)

mutation_list_json_object = json.dumps(mutation_list_json_dict)
with open(output_mutation_list_file, 'w') as output:
    output.write(mutation_list_json_object)