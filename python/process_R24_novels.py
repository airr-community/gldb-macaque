# process R24 novel reads - split by locus, check for presence in the database, write to the correct
# directory if not found

from receptor_utils import simple_bio_seq as simple


infile = 'IGH/bosinger_watson_R24/new_alleles_without_iglabel_3_Dec_2024.csv'
novels = simple.read_csv(infile)


for locus in ['IGH', 'IGK', 'IGL']:
    notes = []
    novels_to_add = {}
    novels_processed = {}
    ld = 'db' if locus == 'IGH' else locus.lower()
    db_recs = simple.read_csv(f'{locus}/db/macaca_mulatta_{ld}.csv')

    seq_lookup = {rec['longest_seq']: rec for rec in db_recs}
    all_db_longest = list(seq_lookup.keys())

    seqs_in_samples = simple.read_csv(f"IGH/bosinger_watson_R24/27_11_{locus}_allele_reference_table.csv")
    seqs_in_samples = {rec['seq']: rec for rec in seqs_in_samples}

    for novel in novels:
        if not novel['seq'] or novel['chain'] != locus:
            continue

        if novel['seq'] in all_db_longest:
            print(f"{novel['allele']} is existing {seq_lookup[novel['seq']]['label']}")
            novel['iglabel'] = seq_lookup[novel['seq']]['label']
            continue

        if novel['seq'] in novels_processed:
            print(f"{novel['allele']} matches {novels_processed[novel['seq']]}")
            continue

        note = {
            'novel': novel['allele'],
            'novel in AIRR-seq': '',
            'novel in genomic': '',
            'sub_or_super': '',
            'sub/super in AIRR-seq': '',
            'sub/super in genomic': '',
        }

        found = False
        for db_seq in all_db_longest:
            if novel['seq'] == db_seq:
                note['sub_or_super'] = 'existing seq'
                found = True
            elif novel['seq'] in db_seq:
                for subseq in seq_lookup[db_seq]['sequences'].split(','):
                    if subseq == novel['seq']:
                        note['sub_or_super'] = 'existing_sub'
                        if db_seq in seqs_in_samples:
                            note['sub/super in genomic'] = seqs_in_samples[db_seq]['sample_count_genomic']
                            note['sub/super in AIRR-seq'] = seqs_in_samples[db_seq]['sample_count_AIRRseq']
                        found = True
                        break
                if not found:
                    note['sub_or_super'] = 'new_sub'
                    if db_seq in seqs_in_samples:
                        note['sub/super in genomic'] = seqs_in_samples[db_seq]['sample_count_genomic']
                        note['sub/super in AIRR-seq'] = seqs_in_samples[db_seq]['sample_count_AIRRseq']
                    found = True
                    break
            elif db_seq in novel['seq']:
                note['sub_or_super'] = 'super'
                if db_seq in seqs_in_samples:
                    note['sub/super in genomic'] = seqs_in_samples[db_seq]['sample_count_genomic']
                    note['sub/super in AIRR-seq'] = seqs_in_samples[db_seq]['sample_count_AIRRseq']
                found = True
                break

        if found:
            if novel['seq'] in seqs_in_samples:
                note['novel in AIRR-seq'] = seqs_in_samples[novel['seq']]['sample_count_AIRRseq']
                note['novel in genomic'] = seqs_in_samples[novel['seq']]['sample_count_genomic']
            notes.append(note)

        if novel['seq'] in novels_processed:
            breakpoint()

        novels_to_add[novel['allele']] = novel['seq']
        novels_processed[novel['seq']] = novel['allele']

    simple.write_fasta(f'{locus}/bosinger_watson_R24/novels.fasta', novels_to_add)
    simple.write_csv(f'{locus}/bosinger_watson_R24/novels_notes.csv', notes)

    used_labels = []

    for novel in novels:
        if novel['iglabel'] and novel['iglabel'] in used_labels:
            print(f"{novel['iglabel']} is duplicated")
        used_labels.append(novel['iglabel'])

    simple.write_csv(f"{infile.replace('.csv', '_updated_with_iglabel.csv')}", novels)

