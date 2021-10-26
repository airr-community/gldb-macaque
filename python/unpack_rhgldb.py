# unpack RhGLDB records to FASTA

import csv
import simple_bio_seq as simple

heavy = {}
light = {}

with open('RhGLDB.csv', 'r') as fi:
    reader = csv.DictReader(fi)
    for row in reader:
        if 'H' in row['Type']:
            heavy[row['Gene Name']] = row['Gene Sequence']
        else:
            light[row['Gene Name']] = row['Gene Sequence']

simple.write_fasta(heavy, 'rhgldb_heavy.fasta')
simple.write_fasta(light, 'rhgldb_light.fasta')