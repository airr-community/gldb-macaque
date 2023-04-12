# compare two rhesus_consolidated databases

import argparse
from receptor_utils import simple_bio_seq as simple

parser = argparse.ArgumentParser(description='Extract reference files for nominated species')
parser.add_argument('db1', help='gapped imgt reference file')
parser.add_argument('db2', help='species name for IMGT file, e.g. Homo sapiens, Macaca mulatta')
args = parser.parse_args()

db1 = simple.read_csv(args.db1)
db2 = simple.read_csv(args.db2)
#simple.write_csv('unspaced'+args.db2, db2)


# make a sequence db and check that sequences are unique

def make_seq_db(db, dbname):
    seqs = {}
    for row in db:
        if row['sequence'] not in seqs:
            seqs[row['sequence']] = []

        seqs[row['sequence']].append(row['gene_label'])

    for k, v in seqs.items():
        if len(v) > 1:
            print(f"{dbname}: the same sequence is listed against multiple labels: {', '.join(v)}")

    return seqs

db1_seqs = make_seq_db(db1, "db1")
db2_seqs = make_seq_db(db2, "db2")

# report differences

for seq1, labels1 in db1_seqs.items():
    if seq1 not in db2_seqs:
        print(f"seq in db1 with labels {', '.join(labels1)} is not in db2")
    else:
        labels2 = db2_seqs[seq1]
        if len(set(labels1).difference(set(labels2))) > 0:
            print(f"seq in db1 with labels {', '.join(labels1)} has labels {', '.join(labels1)} in db2")

for seq2, labels2 in db2_seqs.items():
    if seq2 not in db1_seqs:
        print(f"seq in db2 with labels {', '.join(labels2)} is not in db1")


