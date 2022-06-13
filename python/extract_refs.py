# Extract reference files for nominated species

import argparse
import simple_bio_seq as simple

parser = argparse.ArgumentParser(description='Extract reference files for nominated species')
parser.add_argument('imgt_file', help='gapped imgt reference file')
parser.add_argument('species_name', help='species name for IMGT file, e.g. Homo sapiens, Macaca mulatta')
parser.add_argument('-F', '--functional_only', action='store_true')
args = parser.parse_args()

segs = ['IGHV', 'IGHD', 'IGHJ']
refs = simple.read_imgt_fasta(args.imgt_file, [args.species_name], segs, functional_only=args.functional_only)

simple.write_fasta(refs[args.species_name]['IGHV'], '%s_IGHV_gapped.fasta' %args.species_name.replace(' ', '_'))

ungapped = {}

for seg in segs:
    ungapped[seg] = {}
    for id, seq in refs[args.species_name][seg].items():
        ungapped[seg][id] = seq.replace('.', '')

    simple.write_fasta(ungapped[seg], '%s_%s.fasta' % (args.species_name.replace(' ', '_'), seg))

