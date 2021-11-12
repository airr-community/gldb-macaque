wget http://www.imgt.org/download/GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-nt-WithGaps-F+ORF+inframeP
python ../../python/extract_refs.py IMGTGENEDB-ReferenceSequences.fasta-nt-WithGaps-F+ORF+inframeP "Macaca mulatta"
cat Macaca_mulatta_IGHV.fasta Macaca_mulatta_IGHD.fasta Macaca_mulatta_IGHJ.fasta >Macaca_mulatta_IGH.fasta
cd ../imgt_feb_2021
python ../../python/extract_refs.py imgt_ref_with_gaps_18_feb_2021.fasta "Macaca mulatta"
