# DDM1_manuscript
scripts use to analyze ChIP-seq and BS-seq datasets for Lee, Ipsaro, Adams, Cahn et al. 2023

The main script 'final_script_DDM1.sh' is a wrapper that will analyze ChIP-seq samples using the script 'ChIPseq_primary.sh' (that will call 'ChIPseq_secondary.sh' and 'ChIPseq_tertiary.sh' scripts) and mC samples using the scripts 'MethylC_seq_primary.sh' and 'MethylC_seq_secondary.sh', before plotting.
