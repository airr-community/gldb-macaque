# unpack RhGLDB consolidated light chain records to FASTA

import csv
from receptor_utils import simple_bio_seq as simple

res = {
    'imgt': {},
    'rhgldb': {},
    'genbank': {},
}

with open('consolidated_lambda_LJ_20230419.csv', 'r') as fi:
    reader = csv.DictReader(fi)
    for row in reader:
        for db in ['genbank', 'rhgldb', 'imgt']:
            if row[db]:
                res[db][row[db]] = row['sequence']

for db in ['genbank', 'rhgldb', 'imgt']:
    simple.write_fasta(res[db], f'{db}.fasta')
