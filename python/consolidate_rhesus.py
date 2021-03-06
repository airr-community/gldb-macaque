# Cross-reference allocated labels with the sources to provide a consolidated report

import csv
import receptor_utils.simple_bio_seq as simple
import number_ighv


# Beware that the same sequence may be duplicated in the databases

rhgldb = {}
with open('../cottrell_et_al/gld_rev4_20211111.csv', 'r') as fi:
    reader = csv.DictReader(fi)
    for row in reader:
        if len(row) > 1 and 'H' in row['gene_type']:
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

# IMGT major updates: https://www.imgt.org/IMGTgenedbdoc/dataupdates.html#assemblies

# these two are only used for gapping sequences. The gap change was introduced some time after this set was collected on 18 feb 2021
imgt_v_gapped = simple.read_fasta('../imgt_feb_2021/Macaca_mulatta_IGHV_gapped.fasta')
imgt_v_ungapped = simple.read_fasta('../imgt_feb_2021/Macaca_mulatta_IGHV.fasta')

imgt_functional_only = {v: k for k,v in simple.read_fasta('../imgt/Macaca_mulatta_IGH_29_apr_2022_Functional.fasta').items()}

# The mass change of IMGT IGHV names was announced on 4 September 2020

imgt_files = {
    'post Sept 2021': '../imgt/Macaca_mulatta_IGH_29_apr_2022.fasta',
    'pre Sept 2021': '../imgt_may_2020/Macaca_mulatta_IGH_10_may_2020.fasta',
}

kimdb_accessions = {}
kimdb_accessions['inferred'] = {}
kimdb_accessions['sanger'] = {}

with open('../bernat_et_al/bernat_inferred_sequences.csv', 'r') as fi:
    reader = csv.reader(fi)
    for row in reader:
        if '#' not in row[0] and len(row) > 1:
            kimdb_accessions['inferred'][row[0].replace(' ', '')] = row[1]

with open('../bernat_et_al/bernat_sanger_sequences.csv', 'r') as fi:
    reader = csv.reader(fi)
    for row in reader:
        if '#' not in row[0] and len(row) > 2 and row[2] == 'R':
            kimdb_accessions['sanger'][row[0].replace(' ', '')] = row[1]


seq_notes = {}
seq_notes['rhgldb'] = {}
seq_notes['kimdb'] = {}

used_seq_notes = {}
used_seq_notes['rhgldb'] = []
used_seq_notes['kimdb'] = []

for row in simple.read_csv('../sequence_notes.csv'):
    if 'rhgldb' in row['name']:
        seq_notes['rhgldb'][row['name'].replace('rhgldb:', '')] = row['notes']
    elif 'kimdb' in row['name']:
        seq_notes['kimdb'][row['name'].replace('kimdb:', '')] = row['notes']
    else:
        print('Unrecognised sequence name in notes: %s' % row)



def read_imgt():
    ret = {}

    for qualifier,  filename in imgt_files.items():
        genes = simple.read_fasta(filename)
        for name, seq in genes.items():
            if  seq not in ret:
                ret[seq] = {}

                if 'IGHD' in name:
                    ret[seq]['Type'] = 'IGHD'
                elif 'IGHJ' in name:
                    ret[seq]['Type'] = 'IGHJ'
                elif 'IGHV' in name:
                    ret[seq]['Type'] = 'IGHV'
                else:
                    print('unexpected name in imgt file %s: %s - record dropped' % (filename, name))
                    continue
                ret[seq]['Gene Sequence'] = seq
                ret[seq]['Names'] = {}

            if name not in ret[seq]['Names']:
                ret[seq]['Names'][name] = []

            ret[seq]['Names'][name].append(qualifier)

    for seq, row in ret.items():
        consolidated_names = []
        for name, qualifiers in row['Names'].items():
            if len(qualifiers) == len(imgt_files):
                consolidated_name = name
            elif len(row['Names']) == 1 and len(qualifiers) == 1 and qualifiers[0] == 'post Sept 2021':
                consolidated_name = name
            else:
                consolidated_name = name + ' (' + ', '.join(qualifiers) + ')'

            consolidated_names.append(consolidated_name)

        ret[seq]['Gene Name'] = ', '.join(consolidated_names)
        if seq not in imgt_functional_only:
            ret[seq]['Gene Name'] += ' (P)'

    return ret

imgt = read_imgt()

cottrell_curated_file = '../cottrell_et_al/rhesus_consolidated_edit_20220512_corrected_names.csv'
cottrell_curated = {}
curated_v_recs = 0

with open(cottrell_curated_file, 'r') as fi:
    reader = csv.DictReader(fi)
    for row in reader:
        row['sequence'] = row['sequence'].upper()
        row['sequence_gapped'] = row['sequence_gapped'].upper()
        if row['sequence'] not in cottrell_curated:
            cottrell_curated[row['sequence']] = row
            if row['type'] == 'IGHV':
                curated_v_recs += 1
        elif len(row['sequence']):
            print(f"duplicated sequence in {cottrell_curated_file}: {row['sequence']}")
    cottrell_curated_headers = reader.fieldnames

print(f"Number of curated v_records: {curated_v_recs}")

def add_name_and_type(rec, label, type):
    if rec['type'] is None:
        rec['type'] = type
        rec['gene_label'] = type + '-' + label
    return

all_pubids = []
gene_set_labels = {}
curated_sequences = {}

headers = ['gene_label', 'type', 'functional', 'include_in_set', 'inference_type', 'species_subgroup', 'subgroup_type', 'kimdb', 'rhgldb', 'imgt', 'pubid', 'genbank', 'alt_names', 'notes', 'sequence', 'sequence_gapped']
cottrell_curated_headers = ['CL_rhgldb', 'dataset_count', 'total_VDJs', 'curation_notes']
headers = headers[:3] + cottrell_curated_headers + headers[3:]

consolidated_rows = []
curated_rows_added = 0

with open('macaca_mulatta_db.csv', 'r') as db, open('rhesus_consolidated.csv', 'w', newline='') as fo:
    reader = csv.DictReader(db)
    writer = csv.DictWriter(fo, fieldnames=headers)
    writer.writeheader()
    notfound = 0
    for row in reader:
        for seq in row['sequences'].split(','):
            rearranged = False
            unrearranged = False
            alt_names = []
            notes = []
            rec = {'gene_label': '', 'type': None, 'kimdb': '', 'rhgldb': '', 'imgt': '', 'genbank': '', 'pubid': '',
                   'inference_type': '', 'sequence': seq, 'sequence_gapped': seq, 'notes': '', 'functional': 'Y',
                   'species_subgroup': '', 'subgroup_type': '', 'alt_names': ''}
            genbanks = []
            if seq in kimdb:
                names = []
                for entry in kimdb[seq]:
                    names.append(entry['Gene Name'])
                    alt_names.append('kimdb:' + entry['Gene Name'])
                    add_name_and_type(rec, row['label'], entry['Type'])
                    rearranged = True
                    if entry['Gene Name'] in kimdb_accessions['inferred']:
                        genbanks.append(kimdb_accessions['inferred'][entry['Gene Name']])
                    if entry['Gene Name'] in kimdb_accessions['sanger']:
                        genbanks.append(kimdb_accessions['sanger'][entry['Gene Name']])
                        unrearranged = True
                    if entry['Gene Name'] in seq_notes['kimdb']:
                        notes.append(seq_notes['kimdb'][entry['Gene Name']])
                        used_seq_notes['kimdb'].append(entry['Gene Name'])
                rec['kimdb'] = ', '.join(names)

            if seq in imgt:
                rec['imgt'] = imgt[seq]['Gene Name']
                gene_type = imgt[seq]['Gene Name'][:4]
                add_name_and_type(rec, row['label'], gene_type)
            if seq in rhgldb:
                names = []
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
                    if entry['gene_name'] in seq_notes['rhgldb']:
                        notes.append(seq_notes['rhgldb'][entry['gene_name']])
                        used_seq_notes['rhgldb'].append(entry['gene_name'])

                    if entry['gene_confidence'] == '1' or entry['gene_confidence'] == '2':
                        unrearranged = True
                    else:
                        rearranged = True

                rec['rhgldb'] = ', '.join(names)
                rec['pubid'] = ','.join(pubids)
            if rec['type'] == 'IGHV':
                res, aa, numb_notes = number_ighv.gap_sequence(rec['sequence'], imgt_v_gapped, imgt_v_ungapped)
                rec['sequence_gapped'] = res
                if numb_notes:
                    notes.append(numb_notes)

            rec['notes'] = ', '.join(notes)
            rec['genbank'] = ','.join(genbanks)

            if rearranged and unrearranged:
                rec['inference_type'] = 'Both'
            elif rearranged:
                rec['inference_type'] = 'Rearranged'
            elif unrearranged:
                rec['inference_type'] = 'Unrearranged'

            rec['alt_names'] = ','.join(alt_names)

            if rec['gene_label']:
                if rec['type'] == 'IGHV':
                    if rec['sequence'] in cottrell_curated:
                        rec['include_in_set'] = 'Y'
                        curated_rows_added += 1
                        for el in cottrell_curated_headers:
                            rec[el] = cottrell_curated[rec['sequence']][el]
                    else:
                        rec['include_in_set'] = 'N'

                if 'include_in_set' in rec and rec['include_in_set'] == 'Y':
                    if rec['gene_label'] not in gene_set_labels:
                        gene_set_labels[rec['gene_label']] = 0
                    gene_set_labels[rec['gene_label']] += 1


                if rec['sequence'] not in curated_sequences:
                    curated_sequences[rec['sequence']] = 0
                curated_sequences[rec['sequence']] += 1

                consolidated_rows.append(rec)
            else:
                print('Orphan record: %s' % rec)

    print(f"Number of curated V records added: {curated_rows_added}")

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

    for db_name in ['kimdb', 'rhgldb']:
        notused = []
        for name in seq_notes[db_name]:
            if name not in used_seq_notes[db_name]:
                notused.append(name)

        if len(notused):
            print('Unused notes for %s: %s' % (db_name, ', '.join(notused)))

    for row in cottrell_curated.values():
        if row['type'] == 'IGHV' and row['gene_label'] not in gene_set_labels:
            print(f"curated label not added: {row['gene_label']}: {row['sequence']}")

#all_pubids = list(set(all_pubids))
#for pub in all_pubids:
#    print(pub)
