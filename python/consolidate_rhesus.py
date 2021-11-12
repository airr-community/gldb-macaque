# Cross-reference allocated labels with the sources to provide a consolidated report

import csv
import simple_bio_seq as simple
import number_ighv


# Beware that the same sequence may be duplicated in the databases

rhgldb = {}
with open('../cottrell_et_al/gld_rev4_20211111.csv', 'r') as fi:
    reader = csv.DictReader(fi)
    for row in reader:
        if len(row) > 1 and 'H' in row['gene_type'] and len(row['other_species']) < 1:
            if row['gene_sequence'] not in rhgldb:
                rhgldb[row['gene_sequence']] = []
            if 'DH' in row['gene_type']:
                row['gene_type'] = 'IGHD'
            elif 'JH' in row['gene_type']:
                row['gene_type'] = 'IGHJ'
            elif 'VH' in row['gene_type']:
                row['gene_type'] = 'IGHV'
            else:
                print('unexpected name in kimdb: %s - record dropped' % row['gene_type'])
                continue
            rhgldb[row['gene_sequence']].append(row)


kimdb = {}
genes = simple.read_fasta('../bernat_et_al/Macaca-mulatta_Ig_Heavy_all.fasta')
for name, seq in genes.items():
    row = {}
    row['Gene Name'] = name
    if 'IGHD' in name:
        row['Type'] = 'IGHD'
    elif 'IGHJ' in name:
        row['Type'] = 'IGHJ'
    elif 'IGHV' in name:
        row['Type'] = 'IGHV'
    else:
        print('unexpected name in kimdb: %s - record dropped' % name)
        continue
    row['Gene Sequence'] = seq

    if row['Gene Sequence'] not in kimdb:
        kimdb[row['Gene Sequence']] = []
    kimdb[row['Gene Sequence']].append(row)

# these two are only used for gapping sequences
imgt_v_gapped = simple.read_fasta('../imgt_feb_2021/Macaca_mulatta_IGHV_gapped.fasta')
imgt_v_ungapped = simple.read_fasta('../imgt_feb_2021/Macaca_mulatta_IGHV.fasta')

imgt = {}
genes = simple.read_fasta('../imgt/Macaca_mulatta_IGH.fasta')
for name, seq in genes.items():
    row = {}
    row['Gene Name'] = name
    if 'IGHD' in name:
        row['Type'] = 'IGHD'
    elif 'IGHJ' in name:
        row['Type'] = 'IGHJ'
    elif 'IGHV' in name:
        row['Type'] = 'IGHV'
    else:
        print('unexpected name in imgt: %s - record dropped' % name)
        continue
    row['Gene Sequence'] = seq

    if row['Gene Sequence'] not in imgt:
        imgt[row['Gene Sequence']] = []
    imgt[row['Gene Sequence']].append(row)


def add_name_and_type(rec, label, type):
    if rec['type'] is None:
        rec['type'] = type
        rec['gene_label'] = type + '-' + label

    return

all_pubids = []

headers = ['gene_label', 'type', 'functional', 'inference_type', 'species_subgroup', 'kimdb', 'rhgldb', 'imgt', 'pubid', 'genbank', 'alt_names', 'notes', 'sequence', 'sequence_gapped']

with open('macaca_mulatta_db.csv', 'r') as db, open('rhesus_consolidated.csv', 'w', newline='') as fo:
    reader = csv.DictReader(db)
    writer = csv.DictWriter(fo, fieldnames=headers)
    writer.writeheader()
    notfound = 0
    for row in reader:
        for seq in row['sequences'].split(','):
            if seq not in kimdb and seq not in rhgldb:
                notfound += 1
                continue
            rearranged = False
            unrearranged = False
            alt_names = []
            rec = {'gene_label': '', 'type': None, 'kimdb': '', 'rhgldb': '', 'imgt': '', 'genbank': '', 'pubid': '',
                   'inference_type': '', 'sequence': seq, 'sequence_gapped': seq, 'notes': '', 'functional': 'Y',
                   'species_subgroup': '', 'alt_names': ''}
            if seq in kimdb:
                names = []
                for entry in kimdb[seq]:
                    names.append(entry['Gene Name'])
                    alt_names.append('kimdb:' + entry['Gene Name'])
                    add_name_and_type(rec, row['label'], entry['Type'])
                rec['kimdb'] = ', '.join(names)
                rearranged = True
            if seq in imgt:
                names = []
                for entry in imgt[seq]:
                    names.append(entry['Gene Name'])
                    add_name_and_type(rec, row['label'], entry['Type'])
                rec['imgt'] = ', '.join(names)
            if seq in rhgldb:
                names = []
                genbanks = []
                pubids = []
                for entry in rhgldb[seq]:
                    names.append(entry['gene_name'])
                    alt_names.append('rhgldb:' + entry['gene_name'])
                    if entry['gene_accession']:
                        genbanks.extend(entry['gene_accession'].split(';'))
                    if entry['gene_accession_other']:
                        genbanks.extend(entry['gene_accession_other'].split(';'))
                    pubids.extend(entry['pmid'].split(';'))
                    all_pubids.extend(entry['pmid'].split(';'))
                    add_name_and_type(rec, row['label'], entry['gene_type'])

                if entry['gene_confidence'] == '1' or entry['gene_confidence'] == '2':
                    unrearranged = True
                else:
                    rearranged = True


                rec['rhgldb'] = ', '.join(names)
                rec['genbank'] = ','.join(genbanks)
                rec['pubid'] = ','.join(pubids)
            if rec['type'] == 'IGHV':
                res, aa, notes = number_ighv.gap_sequence(rec['sequence'], imgt_v_gapped, imgt_v_ungapped)
                rec['sequence_gapped'] = res
                if notes:
                    rec['notes'] += notes
            if rearranged and unrearranged:
                rec['inference_type'] = 'Both'
            elif unrearranged:
                rec['inference_type'] = 'Unrearranged Only'
            else:
                rec['inference_type'] = 'Rearranged Only'

            rec['alt_names'] = ','.join(alt_names)
            writer.writerow(rec)

    if notfound > 0:
        print('%d sequences in macaca_mulatta_db.csv were not found in rhesus_consolidated.csv' % notfound)

#all_pubids = list(set(all_pubids))
#for pub in all_pubids:
#    print(pub)
