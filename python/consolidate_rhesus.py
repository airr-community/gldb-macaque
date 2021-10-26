# Cross-reference allocated labels with the sources to provide a consolidated report

import csv
import simple_bio_seq as simple


# Beware that the same sequence may be duplicated in the databases

rhgldb = {}
with open('../cottrell_et_al/RhGLDB.csv', 'r') as fi:
    reader = csv.DictReader(fi)
    for row in reader:
        if 'H' in row['Type']:
            if row['Gene Sequence'] not in rhgldb:
                rhgldb[row['Gene Sequence']] = []
            rhgldb[row['Gene Sequence']].append(row)


kimdb = {}
genes = simple.read_fasta('../bernat_et_al/Macaca-mulatta_Ig_Heavy_all.fasta')
for name, seq in genes.items():
    row = {}
    row['Gene Name'] = name
    if 'IGHD' in name:
        row['Type'] = 'DH'
    elif 'IGHJ' in name:
        row['Type'] = 'JH'
    elif 'IGHV' in name:
        row['Type'] = 'VH'
    else:
        print('unexpected name in kimdb: %s - record dropped' % name)
        continue
    row['Gene Sequence'] = seq

    if row['Gene Sequence'] not in kimdb:
        kimdb[row['Gene Sequence']] = []
    kimdb[row['Gene Sequence']].append(row)


imgt = {}
genes = simple.read_fasta('../imgt/Macaca_mulatta_IGH.fasta')
for name, seq in genes.items():
    row = {}
    row['Gene Name'] = name
    if 'IGHD' in name:
        row['Type'] = 'DH'
    elif 'IGHJ' in name:
        row['Type'] = 'JH'
    elif 'IGHV' in name:
        row['Type'] = 'VH'
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
        rec['gene_label'] = 'IGH' + type[0] + '-' + label

    return

headers = ['gene_label', 'type', 'kimdb', 'rhgldb', 'imgt', 'pubid', 'genbank', 'sequence']

with open('macaca_mulatta_db.csv', 'r') as db, open('rhesus_consolidated.csv', 'w', newline='') as fo:
    reader = csv.DictReader(db)
    writer = csv.DictWriter(fo, fieldnames=headers)
    writer.writeheader()
    for row in reader:
        for seq in row['sequences'].split(','):
            rec = {'gene_label': '', 'type': None, 'kimdb': '', 'rhgldb': '', 'imgt': '', 'genbank': '', 'pubid': '', 'sequence': seq}
            if seq in kimdb:
                names = []
                for entry in kimdb[seq]:
                    names.append(entry['Gene Name'])
                    add_name_and_type(rec, row['label'], entry['Type'])
                rec['kimdb'] = ', '.join(names)
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
                    names.append(entry['Gene Name'])
                    genbanks.append(entry['Gene Accession'])
                    pubids.extend(entry['Citations'].split(';'))
                    add_name_and_type(rec, row['label'], entry['Type'])
                rec['rhgldb'] = ', '.join(names)
                rec['genbank'] = ','.join(genbanks)
                rec['pubid'] = ','.join(pubids)
            writer.writerow(rec)


