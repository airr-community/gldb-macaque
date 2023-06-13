# use Chris's labels in the iglabel database

from receptor_utils import simple_bio_seq as simple

rhgldb = simple.read_csv('../cottrell_et_al/consolidated_lambda_LJ_20230419.csv')

rhgldb_labels = {rec['sequence'].upper(): rec['gene_label'].split('-')[1] for rec in rhgldb}

iglabeldb = simple.read_csv('macaca_mulatta_igl.csv')

res = []
for rec in iglabeldb:
    rec['label'] = rhgldb_labels[rec['longest_seq']]
    res.append(rec)

simple.write_csv('macaca_mulatta_igl_relabelled.csv', res)
