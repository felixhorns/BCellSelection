

python /local10G/rfhorns/Bcell/flu_highres/scripts/lookup_seq_str_uids.py /quakestor4/rfhorns/Bcell/flu_highres/figures/v2/data/Bcell_flu_high_res.sequences.csv 20201010000000500 > sequence_string_uids.txt


grep -e "^.99" < sequence_string_uids.txt > sequence_string_uids_V6_Full.txt

python /local10G/rfhorns/Bcell/flu_highres/scripts/align_ungapped.py /quakestor4/rfhorns/Bcell/flu_highres/figures/v2/data/Bcell_flu_high_res.sequence_strings.csv sequence_string_uids_V6_Full.txt > alignment_ungapped.fasta

python /local10G/rfhorns/Bcell/flu_highres/scripts/trees/lookup_germline.py /quakestor4/rfhorns/Bcell/flu_highres/figures/v2/data/Bcell_flu_high_res.sequences.csv 20201010000000500 ../../resources/Vsegments_20150201.fasta.cleaned ../../resources/Jsegments_20150201.fasta.cleaned > germline.fasta


/local10G/rfhorns/resources/muscle3.8.31/muscle -seqtype dna -maxiters 1 -diags -profile -in1 alignment_ungapped.fasta -in2 germline.fasta -out alignment_ungapped_germline.fasta

/local10G/rfhorns/resources/muscle3.8.31/muscle -seqtype dna -refine -in alignment_ungapped_germline.fasta -out alignment_ungapped_germline_refined.fasta -diags


