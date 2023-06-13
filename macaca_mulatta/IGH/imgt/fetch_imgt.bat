wget --no-check-certificate http://www.imgt.org/download/GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-nt-WithGaps-F+ORF+inframeP
python ../../python/extract_refs.py IMGTGENEDB-ReferenceSequences.fasta-nt-WithGaps-F+ORF+inframeP "Macaca mulatta"
cat Macaca_mulatta_IGHV.fasta Macaca_mulatta_IGHD.fasta Macaca_mulatta_IGHJ.fasta >Macaca_mulatta_IGH_29_apr_2022.fasta
python ../../python/extract_refs.py -F IMGTGENEDB-ReferenceSequences.fasta-nt-WithGaps-F+ORF+inframeP "Macaca mulatta"
cat Macaca_mulatta_IGHV.fasta Macaca_mulatta_IGHD.fasta Macaca_mulatta_IGHJ.fasta >Macaca_mulatta_IGH_29_apr_2022_Functional.fasta
rm Macaca_mulatta_IGHV.fasta Macaca_mulatta_IGHD.fasta Macaca_mulatta_IGHJ.fasta
cd ../imgt_feb_2021
python ../../python/extract_refs.py imgt_ref_with_gaps_18_feb_2021.fasta "Macaca mulatta"
cat Macaca_mulatta_IGHV.fasta Macaca_mulatta_IGHD.fasta Macaca_mulatta_IGHJ.fasta >Macaca_mulatta_IGH_18_feb_2021.fasta
cd ../imgt_may_2020
python ../../python/extract_refs.py imgt_ref_2020_05_10.fasta "Macaca mulatta"
cat Macaca_mulatta_IGHV.fasta Macaca_mulatta_IGHD.fasta Macaca_mulatta_IGHJ.fasta >Macaca_mulatta_IGH_10_may_2020.fasta
