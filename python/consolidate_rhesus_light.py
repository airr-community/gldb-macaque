# Cross-reference allocated labels with the sources to provide a consolidated reportcurated

import csv
from posixpath import split
import receptor_utils.simple_bio_seq as simple
import number_ighv


rhgldb = {}
with open('../cottrell_et_al/consolidated_lambda_LJ_20230419.csv', 'r') as fi:
    reader = csv.DictReader(fi)
    for row in reader:
        if len(row) > 1:
            row['sequence'] = row['sequence'].upper()
            if row['sequence'] not in rhgldb:
                rhgldb[row['sequence']] = []
            if 'IGLV' in row['gene_label']:
                row['gene_type'] = 'IGLV'
            elif 'IGLJ' in row['gene_label']:
                row['gene_type'] = 'IGLJ'
            else:
                print('unexpected name in rhgldb: %s - record dropped' % row['gene_type'])
                continue
            rhgldb[row['sequence']].append(row)



# these two are only used for gapping sequences. Gaps were fixed using receptor_utils
imgt_v_gapped = simple.read_fasta('../imgt/Macaca_mulatta_IGLV_gapped.fasta')
imgt_v_ungapped = simple.read_fasta('../imgt/Macaca_mulatta_IGLV.fasta')


def add_name_and_type(rec, label, type):
    if rec['type'] is None:
        rec['type'] = type
        rec['gene_label'] = type + '-' + label
    return

all_pubids = []
gene_set_labels = {}
curated_sequences = {}

headers = ['gene_label', 'type', 'functional', 'include_in_set', 'inference_type', 'species_subgroup', 'subgroup_type', 'rhgldb', 'imgt', 'pubid', 'genbank', 'alt_names', 'notes', 'sequence', 'sequence_gapped']
consolidated_rows = []
curated_rows_added = 0

with open('macaca_mulatta_IGL.csv', 'r') as db, open('rhesus_consolidated_IGL.csv', 'w', newline='') as fo:
    reader = csv.DictReader(db)
    writer = csv.DictWriter(fo, fieldnames=headers)
    writer.writeheader()
    notfound = 0
    for row in reader:
        for seq in row['sequences'].split(','):
            if seq not in rhgldb:
                notfound += 1
                continue

            rearranged = False
            unrearranged = False
            names = []
            alt_names = []
            notes = []
            pubids = []
            rec = {'gene_label': '', 'type': None, 'rhgldb': '', 'imgt': '', 'genbank': '', 'pubid': '',
                   'inference_type': '', 'sequence': seq, 'sequence_gapped': seq, 'notes': '', 'functional': 'Y',
                   'species_subgroup': '', 'subgroup_type': '', 'alt_names': ''}
            genbanks = []

            for entry in rhgldb[seq]:
                names.append(entry['rhgldb'])
                for db in ['genbank', 'rhgldb', 'imgt']:
                    if entry[db]:
                        alt_names.append(f'{db}: {entry[db]}')
                alt_names.extend(entry['alt_names'].split(';'))
                if entry['IMGT/GenBank']:
                    genbanks.extend(entry['IMGT/GenBank'].split(';'))
                pubids.extend(entry['pubid'].split(';'))
                all_pubids.extend(entry['pubid'].split(';'))
                rec['type'] = entry['gene_label'].split('-')[0]
                rec['gene_label'] = entry['gene_label']

                if '_S' in entry['gene_label']:
                    rearranged = True
                else:
                    unrearranged = True

            rec['rhgldb'] = ', '.join(names)
            rec['pubid'] = ','.join(pubids)
            rec['imgt'] = entry['imgt']

            if 'V' in rec['type']:
                res, aa, numb_notes = number_ighv.gap_sequence(rec['sequence'], imgt_v_gapped, imgt_v_ungapped)
                rec['sequence_gapped'] = res

            rec['notes'] = ', '.join(notes)
            rec['genbank'] = ','.join(genbanks)

            if rearranged and unrearranged:
                rec['inference_type'] = 'Both'
            elif rearranged:
                rec['inference_type'] = 'Rearranged'
            elif unrearranged:
                rec['inference_type'] = 'Unrearranged'

            rec['alt_names'] = ','.join(alt_names)
            rec['include_in_set'] = 'Y'

            if 'include_in_set' in rec and rec['include_in_set'] == 'Y':
                if rec['gene_label'] not in gene_set_labels:
                    gene_set_labels[rec['gene_label']] = 0
                gene_set_labels[rec['gene_label']] += 1


            if rec['sequence'] not in curated_sequences:
                curated_sequences[rec['sequence']] = 0
            curated_sequences[rec['sequence']] += 1

            consolidated_rows.append(rec)

    print(f"Number of curated records added: {len(curated_sequences)}")

    for row in consolidated_rows:
        if row['gene_label'] in gene_set_labels and gene_set_labels[row['gene_label']] > 1:
            cn = row['curation_notes'].split(', ')
            cn.append('duplicated label in set')
            row['curation_notes'] = ', '.join(cn)
            print(f"Label is present more than once in gene set: {row['gene_label']}")
        if curated_sequences[row['sequence']] > 1:
            cn = row['curation_notes'].split(', ')
            cn.append('sequence has >1 label')
            row['curation_notes'] = ', '.join(cn)
            print(f"Sequence has more than one label: {row['sequence']}")
        writer.writerow(row)

    if notfound > 0:
        print('%d sequences in macaca_mulatta_db.csv were not found in rhesus_consolidated.csv' % notfound)
