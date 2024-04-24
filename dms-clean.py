import json
import csv
import argparse


parser = argparse.ArgumentParser()

parser.add_argument('--dms_config', type=str)
parser.add_argument('--auspice_json', type=str)
parser.add_argument('--auspice_json_output', type=str)

args = parser.parse_args()

dms_config = args.dms_config
auspice_json = args.auspice_json
auspice_json_output = args.auspice_json_output


keys = []

with open(dms_config) as f:
    reader = csv.reader(f, delimiter='\t')
    for entry in reader:
        if len(entry) != 0 and not entry[0].startswith("#"):
            try:
                keys.append(entry[2]+' mutations')
            except KeyError:
                pass

with open(auspice_json) as f:
    auspice_json_f = json.load(f)

for entry in auspice_json_f['meta']['colorings']:
    if entry['key'] in keys:
        print('Removing', entry['key'])
        auspice_json_f['meta']['colorings'].remove(entry)
    else:
        print('Leaving', entry['key'])

auspice_json_object = json.dumps(auspice_json_f)
with open(auspice_json_output, 'w') as output:
    output.write(auspice_json_object)