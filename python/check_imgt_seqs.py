# check whether all functional imgt IGH sequences from 2021 are in the database, and also whether they are in the macaque_asc baseline

from receptor_utils import simple_bio_seq as simple

imgt_seqs = simple.read_fasta("d:/research/gldb-macaque/macaca_mulatta/IGH/imgt/Macaca_mulatta_IGH_29_apr_2022_Functional.fasta")
imgt_seqs = {k: v for k, v in imgt_seqs.items() if 'IGHV' in k}
imgt_inv = {v: k for k, v in imgt_seqs.items()}

db_recs = simple.read_csv("d:/research/gldb-macaque/macaca_mulatta/IGH/db/macaca_mulatta_db.csv")
found_imgt = []

for rec in db_recs:
    for seq in rec['sequences'].split(','):
        if seq in imgt_inv:
            found_imgt.append(imgt_inv[seq])

# list items that are in imgt_seqs.keys() but not in found_imgt

all_imgt = set(list(imgt_seqs.keys()))
missing_imgt = all_imgt - set(found_imgt)
print(missing_imgt)

base_seqs = simple.read_fasta("d:/research/macaque_asc/distributed/2024-11-15/IGH/V.fasta")
base_seqs_inv = {v: k for k, v in base_seqs.items()}

missing_imgt = []
found_imgt = []
for imgt_name, seq in imgt_seqs.items():
    if seq not in base_seqs_inv:
        missing_imgt.append(imgt_name)

print(f"there are {len(missing_imgt)} imgt sequences missing from the baseline")
base_recs = simple.read_csv("d:/research/macaque_asc/distributed/2024-11-15/IGH/crossref_with_digger_support.csv")

for miss in missing_imgt:
    miss_seq = imgt_seqs[miss]
    found = False
    for rec in base_recs:
        if miss_seq in rec['seq']:
            print(f"missing {miss} is a subseq of {rec['name']}")
            found = True
        if rec['seq'] in miss_seq:
            print(f"missing {miss} is a superseq of {rec['name']}")
            found = True
        if found:
            break
    if not found:
        print(f"no sub/superseqs of {miss} are in the baseline")


