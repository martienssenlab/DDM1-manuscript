#!/bin/bash
#$ -V
#$ -cwd
#$ -pe threads 24
#$ -l m_mem_free=4G
#$ -l tmp_free=5G
#$ -o mCsinglesample.log
#$ -j y
#$ -N methylcsample

##### Script to run Bismark and analyze Whole-Genome bisulfite data on each sample
##### Called by MethylC_seq_primary.sh
##### It runs fastqc, filter reads, map with bismark/bowtie2, remove duplicates and extract methylation calls
##### It is for alignement of paired-end data and would need to use the SE script for single end

set -e -o pipefail

export threads=$NSLOTS
export limthreads=$((threads/3))
export minthreads=$((threads - 1))

export ref_dir=$1
export name=$2
export samplefilename=$3
export met=$4
export step=$5

printf "\n\n"
date
printf "\n"

if [[ ${met} == "Pico" ]]; then
	param1="-j $threads -u 10 -U 10 -q 10 -m 20"
	param2="--non_directional --maxins 1000"
	param3=""
elif [[ ${met} == "EM" ]]; then
	param1="-j $threads -m 20"
	param2="--maxins 1000 --maxins 1000" 
	param3="--ignore_r2 2"
else
	printf "Missing information about kit used\nDefaulting to not pico or EM!\n"
	param1="-j $threads -m 20"
	param2="--maxins 1000 --maxins 1000" 
	param3="--ignore_r2 2"
fi

if [[ ${step} == "trim" ]]; then
	printf "\nRunning fastQC for ${name}\n"
	fastqc -o reports/ fastq/${name}_R1.fastq.gz
	fastqc -o reports/ fastq/${name}_R2.fastq.gz
	printf "\nRunning Cutadapt for ${name}\n"
	cutadapt ${param1} -a AGATCGGAAGAGCACACGTCTGAAC -A AGATCGGAAGAGCGTCGTGTAGGGA -o fastq/trimmed_${name}_R1.fastq.gz -p fastq/trimmed_${name}_R2.fastq.gz fastq/${name}_R1.fastq.gz fastq/${name}_R2.fastq.gz |& tee reports/trimming_${name}.txt
	printf "\nRunning fastQC on trimmed files for ${name}\n"
	fastqc -o reports/ fastq/trimmed_${name}_R1.fastq.gz
	fastqc -o reports/ fastq/trimmed_${name}_R2.fastq.gz
	rm -f fastq/${name}_R1.fastq.gz
	rm -f fastq/${name}_R2.fastq.gz
	step="map"
fi

if [[ ${step} == "map" ]]; then
	R1="fastq/trimmed_${name}_R1.fastq.gz"
	R2="fastq/trimmed_${name}_R2.fastq.gz"
	shortR1="trimmed_${name}_R1"
	printf "\nAligning ${name} with bismark_bowtie2\n"
	bismark --genome ${ref_dir} ${param2} --local --multicore ${limthreads} --temp_dir=${TMPDIR} -o mapped/${name} --gzip --nucleotide_coverage -1 ${R1} -2 ${R2} |& tee reports/alignment_bismark_${name}.txt
	printf "\nDeduplicating with bismark\n"
	deduplicate_bismark -p --output_dir mapped/${name}/ -o ${name} --bam mapped/${name}/${shortR1}_bismark_bt2_pe.bam |& tee reports/deduplication_bismark_${name}.txt
	rm -f mapped/${name}/${shortR1}_bismark_bt2_pe.bam
	printf "\nCalling mC for ${name}"
	bismark_methylation_extractor -p --comprehensive -o methylcall/ ${param3} --gzip --multicore ${limthreads} --buffer_size 10G --cytosine_report --CX --genome_folder ${ref_dir} mapped/${name}/${name}.deduplicated.bam
	rm -f methylcall/C*context_${name}*
	rm -f methylcall/${name}*bismark.cov*
	printf "\nMaking final html report for ${name}\n"
	bismark2report -o final_report_${name}.html --dir reports/ --alignment_report mapped/${name}/trimmed_${name}_R1_bismark_bt2_PE_report.txt --dedup_report mapped/${name}/trimmed_${name}_R1_bismark_bt2_pe.deduplication_report.txt --splitting_report methylcall/${name}.deduplicated_splitting_report.txt --mbias_report methylcall/${name}.deduplicated.M-bias.txt --nucleotide_report mapped/${name}/trimmed_${name}_R1_bismark_bt2_pe.nucleotide_stats.txt
	
	printf "\nCalculting coverage stats for ${name}\n"
	zcat methylcall/${name}.deduplicated.CX_report.txt.gz | awk -v OFS="\t" -v s=${name} '{a+=1; b=$4+$5; g+=b; if (b>0) {c+=1; d+=b;} else f+=1; if (b>2) e+=1} END {print s,a,f/a*100,c/a*100,e/a*100,g/a,d/c,n}' >> reports/summary_coverage.txt
	printf "\nMaking bedGraph files of each context for ${name}\n"
	zcat methylcall/${name}.deduplicated.CX_report.txt.gz | awk -v OFS="\t" -v s=${name} '($4+$5)>0 {a=$4+$5; if ($6=="CHH") print $1,$2-1,$2,$4/a*100 > "methylcall/"s"_CHH.bedGraph"; else if ($6=="CHG") print $1,$2-1,$2,$4/a*100 > "methylcall/"s"_CHG.bedGraph"; else print $1,$2-1,$2,$4/a*100 > "methylcall/"s"_CG.bedGraph"}'
	for context in CG CHG CHH
	do
		printf "\nMaking bigwig files of ${context} context for ${name}\n"
		LC_COLLATE=C sort -k1,1 -k2,2n methylcall/${name}_${context}.bedGraph > methylcall/sorted_${name}_${context}.bedGraph
		bedGraphToBigWig methylcall/sorted_${name}_${context}.bedGraph ${ref_dir}/chrom.sizes methylcall/${name}_${context}.bw
	done
	rm -f methylcall/*${name}*bedGraph*
	printf "\nGetting some stats with samtools\n"
	samtools stats -@ ${minthreads} mapped/${name}/${name}.deduplicated.bam > reports/alignment_samtools_stats_${name}.txt
	step="stats"
fi

if [[ ${step} == "stats" ]]; then
	if grep -q "L09137.2" ${ref_dir}/chrom.sizes; then
		printf "\nCalculting conversion stats for ${name} using EMseq controls\n"
		tot=$(cat reports/alignment_bismark_${name}.txt | grep "Sequence pairs analysed in total:" | awk -v FS="\t" 'END {print $2}')
		map=$(cat reports/alignment_bismark_${name}.txt | grep "Number of paired-end alignments with a unique best hit:" | awk -v FS="\t" 'END {print $2}')
		zcat methylcall/${name}.deduplicated.CX_report.txt.gz | awk -v OFS="\t" -v s=${name} -v t=${tot} -v r=${map} '{b=$4+$5; if ($1=="J02459.1_48502") {m+=$4; n+=b;}; if ($1=="L09137.2" && $6=="CG") {o+=$4; p+=b;}} END {print s,t,r" ("r/t*100")",m/n*100,o/p*100}' >> reports/summary_controls_${samplefilename}.txt
	elif [ -s methylcall/control_${name}.deduplicated.CX_report.txt.gz ]; then
		printf "\nCalculting conversion stats for ${name} using EMseq controls\n"
		tot=$(cat reports/alignment_bismark_control_${name}.txt | grep "Sequence pairs analysed in total:" | awk -v FS="\t" 'END {print $2}')
		map=$(cat reports/alignment_bismark_control_${name}.txt | grep "Number of paired-end alignments with a unique best hit:" | awk -v FS="\t" 'END {print $2}')
		zcat methylcall/control_${name}.deduplicated.CX_report.txt.gz | awk -v OFS="\t" -v s=${name} -v t=${tot} -v r=${map} '{b=$4+$5; if ($1=="J02459.1_48502") {m+=$4; n+=b;}; if ($1=="L09137.2" && $6=="CG") {o+=$4; p+=b;}} END {print s,t,r" ("r/t*100")",m/n*100,o/p*100}' >> reports/summary_controls_${samplefilename}.txt
	fi
	if grep -E -q "Pt|ChrC|chrc" ${ref_dir}/chrom.sizes; then
		printf "\nCalculting coverage stats for ${name}\n"
		zcat methylcall/${name}.deduplicated.CX_report.txt.gz | awk -v OFS="\t" -v s=${name} '{a+=1; b=$4+$5; g+=b; if (b>0) {c+=1; d+=b;} else f+=1; if (b>2) e+=1; if ($1 == "Pt" || $1 == "ChrC" || $1 == "chrC") {m+=$4; n+=b;}} END {print s,a,f/a*100,c/a*100,e/a*100,g/a,d/c,m/n*100}' >> reports/summary_coverage_${samplefilename}.txt
	elif grep -E -q "J02459.1_48502" ${ref_dir}/chrom.sizes; then
		printf "\nCalculting coverage stats for ${name}\n"
		nonc=$(awk -v s=${name} '$1==s {print $5}' reports/summary_controls_${samplefilename}.txt)
		zcat methylcall/${name}.deduplicated.CX_report.txt.gz | awk -v OFS="\t" -v s=${name} -v n=${nonc} '{a+=1; b=$4+$5; g+=b; if (b>0) {c+=1; d+=b;} else f+=1; if (b>2) e+=1;} END {print s,a,f/a*100,c/a*100,e/a*100,g/a,d/c,n}' >> reports/summary_coverage_${samplefilename}.txt
	fi
fi

printf "Script finished!\n"
touch chkpts/${name}
