from Bio import SeqIO
from Bio.Seq import Seq
import csv
import json
import argparse


parser = argparse.ArgumentParser()

parser.add_argument('--reference', type=str, help='path to reference genbank file')
parser.add_argument('--nt', type=str, help='path to nt-muts json file')
parser.add_argument('--summary', type=str, help='path to DMS csv file')
parser.add_argument('--site_header', type=str, help='header of column with amino acid sites, in h5 sequential numbering')
parser.add_argument('--h3_site_header', type=str, help='header of column with amino acid sites, in h3 sequential numbering (optional)', required=False)
parser.add_argument('--dms_config', type=str)
parser.add_argument('--auspice_config', type=str, help='optional path to config json file; if provided, will automatically populate colorings with DMS keys', required=False)
parser.add_argument('--output_totals', type=str, help='path to output json file with total DMS values')
parser.add_argument('--output_summary', type=str, help='path to output json file with summary of DMS mutations')
parser.add_argument('--output_list', type=str, help='path to output json file with mutations listed as string for tree')
parser.add_argument('--output_auspice_config', type=str, required=False)

args = parser.parse_args()

reference = args.reference
nt = args.nt
summary = args.summary
site_header = args.site_header
h3_site_header = args.h3_site_header
dms_config = args.dms_config
auspice_config = args.auspice_config
output_totals = args.output_totals
output_summary = args.output_summary
output_list = args.output_list
output_auspice_config = args.output_auspice_config


min_hex = '#DC2F24' # '#4041C7'
mid_hex =  '#ABABAB'
max_hex = '#4041C7' # '#DC2F24'


# determine start and end positions of HA CDS from reference
for r in SeqIO.parse(reference, 'genbank'):
    for f in r.features:
        if f.type == 'CDS' and f.qualifiers['gene'][0] == 'HA':
            start, end = f.location.start, f.location.end


# convert dms config tsv into a dict
dms_config_dict = {}

with open(dms_config, 'r') as f:
    reader = csv.reader(f, delimiter='\t')
    for entry in reader:
        if len(entry) != 0 and not entry[0].startswith("#"):
            try:
                dms_config_dict[entry[3]] = {'header': entry[0], 'filter': entry[1], 'title': entry[2]}
            except KeyError:
                pass


# convert the escape summary csv into a nested dictionary
# this makes it much faster to pull the escape values than filtering the entire df each time
value_dict = {k: {} for k in dms_config_dict.keys()}
h5_to_h3 = {}

with open(summary, 'r') as f:
    # read summary csv
    reader = csv.reader(f)
    # pull out header names
    headers = next(reader)
    # find indices of columns needed -- if header order changes, will not need to manually edit
    # header 'mutant' must maintain naming
    # site and DMS value headers supplied as args
    site_index = headers.index(site_header)
    mutant_index = headers.index('mutant')
    if h3_site_header:
        h3_site_index = headers.index(h3_site_header)
    # for each dms header supplied in dms config, determine indices
    for k,v in dms_config_dict.items():
        value_header = v['header']
        dms_config_dict[k]['index'] = headers.index(value_header)
    # iterate through rows of dms summary with phenotype values for each mutant position
    for row in reader:
        site = int(row[site_index])
        mutant = row[mutant_index]
        if h3_site_header and site not in h5_to_h3:
            h3_site = row[h3_site_index]
            h5_to_h3[str(site)] = h3_site
        # for each phenotype supplied in dms config, determine mutant + value pairs
        for key,v in dms_config_dict.items():
            value_index = v['index']
            value_filter = v['filter']
            # if the value is blank, set to 0 — otherwise, pull out float value
            value = float(row[value_index]) if row[value_index] != '' else 0
            # if the value is not 0, determine if the value must be filtered out
            # if value is 0, does not need to be added to the dict (no phenotype for that mutation)
            if value != 0:
                # if the value matches the phenotype's filter, or the phenotype's values do not need filtered, add the value to the dict
                # otherwise, do not add to the dict
                if ((value_filter == '+') & (value > 0)) | ((value_filter == '-') & (value < 0)) | ((value_filter != '+') & (value_filter != '-')):
                    if site in value_dict[key]:
                        value_dict[key][site][mutant] = value
                    else:
                        value_dict[key][site] = {mutant: value}
                

# for each node and tip on the tree, determine the amino acid sequence
# save into a node:aa dict
aa_seqs = {}

with open(nt, 'r') as f:
    # open nt json file
    nt_file = json.load(f)
    # iterate through nodes/tips
    for node, value in nt_file['nodes'].items():
        # pull out cds
        cds_str = value['sequence'][start:end]#[start:start+1032]+value['sequence'][start+1035:end] # deletion in DMS strain at positions 1033-1035 of HA relative to GsGd (1032:1035 relative to start of HA)
        # replace gaps with n for translation
        cds = Seq(cds_str).replace('-','n')
        # translate to aa
        aa = cds.translate()
        # add to dict
        aa_seqs[node] = str(aa)


# for each node and tip, determine the total dms value given the amino acid sequence
dms_values = {}
# also save the mutations and associated values
mutation_summary_dict = {}

# iterate through amino acid sequences
for name, aa_seq in aa_seqs.items():
    # add empty dicts for each node/tip
    dms_values[name] = {}
    mutation_summary_dict[name] = {}
    # iterate through all of the supplied dms phenotypes
    for key in value_dict.keys():
        # phenotype total starts at 0
        total = 0
        # site starts at 0
        site = 0
        # dict to save mutations and values to
        mutations = {}
        # iterate through each position in the aa sequence
        for aa in aa_seq:
            # add 1 to the site (for 1-indexed positions)
            site = site + 1
            # if the mutation is found in escape_dict, add it to the strain's total
            try:
                total += value_dict[key][site][aa]
                # also save the mutation and its associated value
                mutations[str(site)+aa] = value_dict[key][site][aa]
            # if it's not in escape_dict, pass -- value either missing or was filtered out
            except KeyError:
                pass
        # save total in the dict after rounding
        dms_values[name][key+'_dms_value'] = round(total,5)
        # and save the mutation's identities and values
        mutation_summary_dict[name][key] = mutations

# find min and max values for each phenotype for coloring
dms_values_min_max = {k:{'min':None,'max':None} for k in value_dict.keys()}

for k,v in dms_values_min_max.items():
    min_value = min(dms_values[strain][k+'_dms_value'] for strain in dms_values.keys())
    max_value = max(dms_values[strain][k+'_dms_value'] for strain in dms_values.keys())
    if min_value < 0:
        dms_values_min_max[k]['min'] = min_value
    if max_value > 0:
        dms_values_min_max[k]['max'] = max_value

# for each node/tip, save mutations + values as a string for displaying on tree
mutation_list_dict = {}

# iterate through nodes/tips
for node, mut_value_dicts in mutation_summary_dict.items():
    # add empty dict for the node/tip
    mutation_list_dict[node] = {}
    # iterate through the node/tip's mutations
    for key, mut_value_dict in mut_value_dicts.items():
        # convert mutations + values into a string and save in a list
        outlist = []
        for mut, esc in sorted(mut_value_dict.items(), key=lambda x:x[1], reverse=True):
            if h3_site_header:
                h5_mut = h5_to_h3[mut[:-1]] + mut[-1]
                outlist.append(f'{mut} ({h5_mut}) [{esc}]')
            else:
                outlist.append(f'{mut} [{esc}]')
        # combine all the strings
        if len(outlist) != 0:
            outstr = 'H5 (H3) [effect]: ' + ', '.join(outlist)
        else:
            outstr = ''
        # save the final combined string
        mutation_list_dict[node][dms_config_dict[key]['title']+' mutations'] = outstr


# format dict for auspice
dms_json_dict = {'nodes': dms_values}
mutation_summary_json_dict = {'nodes': mutation_summary_dict}
mutation_list_json_dict = {'nodes': mutation_list_dict}


# save dicts as json files
dms_json_object = json.dumps(dms_json_dict)
with open(output_totals, 'w') as output:
    output.write(dms_json_object)

mutation_summary_json_object = json.dumps(mutation_summary_json_dict)
with open(output_summary, 'w') as output:
    output.write(mutation_summary_json_object)

mutation_list_json_object = json.dumps(mutation_list_json_dict)
with open(output_list, 'w') as output:
    output.write(mutation_list_json_object)

# if auspice config is supplied, add colorings for each dms phenotype
if auspice_config:
    with open(auspice_config, 'r') as f:
        auspice_config_json = json.load(f)

    for k,v in dms_config_dict.items():
        min_value = dms_values_min_max[k]['min']
        max_value = dms_values_min_max[k]['max']
        if min_value and max_value:
            abs_min = min(min_value, -max_value)
            abs_max = -abs_min
            scale = [[abs_min, min_hex], [0, mid_hex], [abs_max, max_hex]]
        elif min_value:
            scale = [[min_value, min_hex], [0, mid_hex]]
        elif max_value:
            scale = [[0, mid_hex], [max_value, max_hex]]
        else:
            scale = [[0, mid_hex]]
        auspice_config_json['colorings'].append({'key': k+'_dms_value',
                            'title': v['title'],
                            'type': 'continuous',
                            'scale': scale})
        auspice_config_json['colorings'].append({'key': v['title']+' mutations',
                            'title': v['title']+' mutations',
                            'type': 'categorical'})
    
    auspice_config_json_object = json.dumps(auspice_config_json)

    # if an output file is specified, save the updated auspice config there
    if output_auspice_config:
        with open(output_auspice_config, 'w') as output:
            output.write(auspice_config_json_object)
    # otherwise, overwrite initial config
    # this is pretty bad if the build is re-run without reverting the config, since additional colorings will be appended
    else:
        with open(auspice_config, 'w') as output:
            output.write(auspice_config_json_object)