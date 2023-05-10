# DDM1_manuscript
scripts use to analyze ChIP-seq and BS-seq datasets for Lee, Ipsaro, Adams, Cahn et al. 2023

The main script 'final_script_DDM1.sh' is a wrapper that will analyze ChIP-seq samples using the script 'maizecode.sh' (that will call '..ChIP_sample' and '..ChIP_analysis.sh' scripts) and mC samples using the scripts 'MethylC_seq_primary.sh' and 'MethylC_seq_bismark_singlesample_PE.sh', before plotting.
