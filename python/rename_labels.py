# find all recs in db1 with a seq_id and new_label in map, rewrite the record with new label into db2

import argparse
from receptor_utils import simple_bio_seq as simple
from receptor_utils import number_ighv

parser = argparse.ArgumentParser(description='Rename labels in a consolidated db')
parser.add_argument('db1')
parser.add_argument('db2')
parser.add_argument('map')
parser.add_argument('v_gapped')
parser.add_argument('v_ungapped')
args = parser.parse_args()

db1 = simple.read_csv(args.db1)
map = simple.read_csv(args.map)
map = {m['seq_id']: m['new_label'] for m in map if m['new_label']}

v_gapped = simple.read_fasta(args.v_gapped)
v_ungapped = simple.read_fasta(args.v_ungapped)

new_recs = []

for rec in db1:
    if rec['gene_label'] in map:
        rec['gene_label'] = 'IGHV-' + map[rec['gene_label']]
        res, aa, numb_notes = number_ighv.gap_sequence(rec['sequence'], v_gapped, v_ungapped)
        rec['sequence_gapped'] = res
        if numb_notes:
            if rec['notes']:
                rec['notes'] += ', '
            rec['notes'] += numb_notes
        new_recs.append(rec)

simple.write_csv(args.db2, new_recs)


