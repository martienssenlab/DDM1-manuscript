#!/bin/bash
#$ -V
#$ -cwd
#$ -pe threads 12
#$ -l m_mem_free=2G
#$ -l tmp_free=4G
#$ -o logDDM1.log
#$ -j y
#$ -N DDM1

set -e -o pipefail

printf "\n\n"
date
printf "\n"

export threads=$NSLOTS

#############################################################################
###################### To create the samplefile #############################
#############################################################################

# printf "
# ChIP_JI\tCol0\tWT\tDDM1\tRep1\t_S4_\t/grid/martienssen/data_nlsas_norepl/archive/data/FCAAAWYHNHV\tPE\tTAIR10\n
# ChIP_JI\tCol0\tWT\tDDM1\tRep2\t_S8_\t/grid/martienssen/data_nlsas_norepl/archive/data/FCAAAWYHNHV\tPE\tTAIR10\n
# ChIP_JI\tCol0\tddm1\tDDM1\tRep1\t_S20_\t/grid/martienssen/data_nlsas_norepl/archive/data/FCAAAWYHNHV\tPE\tTAIR10\n
# ChIP_JI\tCol0\tddm1\tDDM1\tRep2\t_S24_\t/grid/martienssen/data_nlsas_norepl/archive/data/FCAAAWYHNHV\tPE\tTAIR10\n
# ChIP_K\tCol0\tWT\tH3K27me1\tRep1\twt.k27me1.ip.r1\t/grid/martienssen/data_norepl/bberube/arabidopsis/ddm1/geo_submission/K27me1/raw\tSE\tTAIR10\n
# ChIP_K\tCol0\tWT\tH3K27me1\tRep2\twt.k27me1.ip.r2\t/grid/martienssen/data_norepl/bberube/arabidopsis/ddm1/geo_submission/K27me1/raw\tSE\tTAIR10\n
# ChIP_K\tCol0\tddm1\tH3K27me1\tRep1\tddm1.k27me1.ip.r1\t/grid/martienssen/data_norepl/bberube/arabidopsis/ddm1/geo_submission/K27me1/raw\tSE\tTAIR10\n
# ChIP_K\tCol0\tddm1\tH3K27me1\tRep2\tddm1.k27me1.ip.r2\t/grid/martienssen/data_norepl/bberube/arabidopsis/ddm1/geo_submission/K27me1/raw\tSE\tTAIR10\n
# ChIP_J4\tCol0\tWT\tH4K16ac Rep1\t_S1_\t/grid/martienssen/data_nlsas_norepl/archive/data/FCAAAWYHNHV\tPE\tTAIR10\n
# ChIP_J4\tCol0\tWT\tH4K16ac Rep2\t_S5_\t/grid/martienssen/data_nlsas_norepl/archive/data/FCAAAWYHNHV\tPE\tTAIR10\n
# ChIP_J4\tCol0\tddm1\tH4K16ac\tRep1\t_S17_\t/grid/martienssen/data_nlsas_norepl/archive/data/FCAAAWYHNHV\tPE\tTAIR10\n
# ChIP_J4\tCol0\tddm1\tH4K16ac\tRep2\t_S21_\t/grid/martienssen/data_nlsas_norepl/archive/data/FCAAAWYHNHV\tPE\tTAIR10\n
# ChIP_H\tCol0\tWT\tH3.3\tRep1\twt.htr5.ip.r1\t/grid/martienssen/data_norepl/bberube/arabidopsis/ddm1/geo_submission/HTR5/raw\tSE\tTAIR10\n
# ChIP_H\tCol0\tWT\tH3.3\tRep2\twt.htr5.ip.r2\t/grid/martienssen/data_norepl/bberube/arabidopsis/ddm1/geo_submission/HTR5/raw\tSE\tTAIR10\n
# ChIP_H\tCol0\tddm1\tH3.3\tRep1\tddm1.htr5.ip.r1\t/grid/martienssen/data_norepl/bberube/arabidopsis/ddm1/geo_submission/HTR5/raw\tSE\tTAIR10\n
# ChIP_H\tCol0\tddm1\tH3.3\tRep2\tddm1.htr5.ip.r2\t/grid/martienssen/data_norepl/bberube/arabidopsis/ddm1/geo_submission/HTR5/raw\tSE\tTAIR10\n
# ChIP_M\tCol0\tWT\tMGH3\tRep1\twt.mgh3.ip\t/grid/martienssen/home/jcahn/norepl/projects/fin_DDM1/MGH3_fastq_se SE\tTAIR10\n
# ChIP_M\tCol0\tddm1\tMGH3\tRep1\tddm1.mgh3.ip\t/grid/martienssen/home/jcahn/norepl/projects/fin_DDM1/MGH3_fastq_se\tSE\tTAIR10\n
# ChIP_JI\tCol0\tWT\tInput\tRep1\t_S25_\t/grid/martienssen/data_nlsas_norepl/archive/data/FCAAAWYHNHV\tPE\tTAIR10\n
# ChIP_JI\tCol0\tWT\tInput\tRep2\t_S26_\t/grid/martienssen/data_nlsas_norepl/archive/data/FCAAAWYHNHV\tPE\tTAIR10\n
# ChIP_JI\tCol0\tddm1\tInput\tRep1\t_S29_\t/grid/martienssen/data_nlsas_norepl/archive/data/FCAAAWYHNHV\tPE\tTAIR10\n
# ChIP_JI\tCol0\tddm1\tInput\tRep2\t_S30_\t/grid/martienssen/data_nlsas_norepl/archive/data/FCAAAWYHNHV\tPE\tTAIR10\n
# ChIP_K\tCol0\tWT\tInput\tRep1\twt.k27me1.h3.r1\t/grid/martienssen/data_norepl/bberube/arabidopsis/ddm1/geo_submission/K27me1/raw\tSE\tTAIR10\n
# ChIP_K\tCol0\tWT\tInput\tRep2\twt.k27me1.h3.r2\t/grid/martienssen/data_norepl/bberube/arabidopsis/ddm1/geo_submission/K27me1/raw\tSE\tTAIR10\n
# ChIP_K\tCol0\tddm1\tInput\tRep1\tddm1.k27me1.h3.r1\t/grid/martienssen/data_norepl/bberube/arabidopsis/ddm1/geo_submission/K27me1/raw\tSE\tTAIR10\n
# ChIP_K\tCol0\tddm1\tInput\tRep2\tddm1.k27me1.h3.r2\t/grid/martienssen/data_norepl/bberube/arabidopsis/ddm1/geo_submission/K27me1/raw\tSE\tTAIR10\n
# ChIP_J4\tCol0\tWT\tInput\tRep1\t_S2_\t/grid/martienssen/data_nlsas_norepl/archive/data/FCAAAWYHNHV\tPE\tTAIR10\n
# ChIP_J4\tCol0\tWT\tInput\tRep2\t_S6_\t/grid/martienssen/data_nlsas_norepl/archive/data/FCAAAWYHNHV\tPE\tTAIR10\n
# ChIP_J4\tCol0\tddm1\tInput\tRep1\t_S18_\t/grid/martienssen/data_nlsas_norepl/archive/data/FCAAAWYHNHV\tPE\tTAIR10\n
# ChIP_J4\tCol0\tddm1\tInput\tRep2\t_S22_\t/grid/martienssen/data_nlsas_norepl/archive/data/FCAAAWYHNHV\tPE\tTAIR10\n
# ChIP_H\tCol0\tWT\tInput\tRep1\twt.htr5.h3.r1\t/grid/martienssen/data_norepl/bberube/arabidopsis/ddm1/geo_submission/HTR5/raw\tSE\tTAIR10\n
# ChIP_H\tCol0\tWT\tInput\tRep2\twt.htr5.h3.r2\t/grid/martienssen/data_norepl/bberube/arabidopsis/ddm1/geo_submission/HTR5/raw\tSE\tTAIR10\n
# ChIP_H\tCol0\tddm1\tInput\tRep1\tddm1.htr5.h3.r1\t/grid/martienssen/data_norepl/bberube/arabidopsis/ddm1/geo_submission/HTR5/raw\tSE\tTAIR10\n
# ChIP_H\tCol0\tddm1\tInput\tRep2\tddm1.htr5.h3.r2\t/grid/martienssen/data_norepl/bberube/arabidopsis/ddm1/geo_submission/HTR5/raw\tSE\tTAIR10\n
# ChIP_M\tCol0\tWT\tInput\tRep1\twt.mgh3.in\t/grid/martienssen/home/jcahn/norepl/projects/fin_DDM1/MGH3_fastq_se\tSE\tTAIR10\n
# ChIP_M\tCol0\tddm1\tInput\tRep1\tddm1.mgh3.in\t/grid/martienssen/home/jcahn/norepl/projects/fin_DDM1/MGH3_fastq_se\tSE\tTAIR10\n
# ChIP_S1\tCol0\tWT_Stroud\tH3.3\tRep1\tSRR394078\tSRA\tSE\tTAIR10\n
# ChIP_S2\tCol0\tWT_Stroud\tH3.1\tRep1\tSRR394079\tSRA\tSE\tTAIR10\n
# ChIP_S1\tCol0\tWT_Stroud\tInput\tRep1\tSRR394080\tSRA\tSE\tTAIR10\n
# ChIP_S2\tCol0\tWT_Stroud\tInput\tRep1\tSRR394081\tSRA\tSE\tTAIR10\n
# ChIP_L\tCol0\tWT_borg\tMGH3\tRep1\tSRR7945256\tSRA\tPE\tTAIR10\n
# ChIP_L\tCol0\tWT_borg\tInput\tRep1\tSRR7945253\tSRA\tPE\tTAIR10\n
# " > final_samplefile.txt

#############################################################################
############ To prepare the samples with the Maizecode pipeline #############
#############################################################################

###### To run the maizecode pipeline (mapping + generating tracks + stats + general plots)

# pids=()
# qsub -sync y -N maizecode -o logDDM1.log scripts/MaizeCode.sh -f final_samplefile.txt -p /grid/martienssen/home/jcahn/nlsas/Genomes/Arabidopsis -a all -x -m DDM1 &
# pids+=("$!")

# wait ${pids}

# ref_dir="/grid/martienssen/home/jcahn/nlsas/Genomes/Arabidopsis/TAIR10"

#############################################################################
################# To make the ddm1 vs WT bigwig files #######################
#############################################################################

# for mark in DDM1 H3.3 H3K27me1 H4K16ac MGH3
# do
	# if [ -e ChIP/mapped/Col0_ddm1_${mark}_merged.bam ] && [ -e ChIP/mapped/Col0_WT_${mark}_merged.bam ]; then
		# printf "Making ddm1 vs WT ${mark} merged\n"
		# bamCompare -p ${threads} -bs 10 --scaleFactorsMethod "None" --normalizeUsing "CPM" -b1 ChIP/mapped/Col0_ddm1_${mark}_merged.bam -b2 ChIP/mapped/Col0_WT_${mark}_merged.bam -o ChIP/tracks/ddm1_vs_WT_${mark}_cpm.bw
	# elif [ -e ChIP/mapped/Col0_ddm1_${mark}_Rep1.bam ] && [ -e ChIP/mapped/Col0_WT_${mark}_Rep1.bam ]; then
		# printf "Making ddm1 vs WT ${mark} Rep1\n"
		# bamCompare -p ${threads} -bs 10 --scaleFactorsMethod "None" --normalizeUsing "CPM" -b1 ChIP/mapped/Col0_ddm1_${mark}_Rep1.bam -b2 ChIP/mapped/Col0_WT_${mark}_Rep1.bam -o ChIP/tracks/ddm1_vs_WT_${mark}_Rep1_cpm.bw
	# fi
# done

#############################################################################
######### Figure5 and FigureS7 - Heatmaps at DMRs and Chr5 TEs ##############
#############################################################################

# awk '$3 != "WT_Stroud"' final_analysis_samplefile.txt > subset_analysis_samplefile.txt

# if [ ! -d ./figures_manuscript ]; then
	# mkdir ./figures_manuscript
	# mkdir ./figures_manuscript/draft
	# mkdir ./figures_manuscript/data
# fi

# ### To use the log(IP vs Input) for heatmaps

# samplelist1=()
# label_list1=()
# ordered_mark1=()
# while read type ref sample mark paired ref_dir
# do
	# name="${sample}_${mark}"
	# merged="ChIP/tracks/${ref}_${name}_merged.bw"
	# rep="ChIP/tracks/${ref}_${name}_Rep1.bw"
	# if [ -s ${merged} ]; then
		# samplelist1+=("${merged}")
	# else
		# samplelist1+=("${rep}")
	# fi
	# label_list1+=("${name}")
	# if [[ ! "${ordered_mark1[@]}" =~ "${mark}" ]]; then
		# ordered_mark1+=("${mark}")
	# fi
# done < subset_analysis_samplefile.txt
# output1="subset_DDM1"

# ref_dir="/grid/martienssen/home/jcahn/nlsas/Genomes/Arabidopsis/TAIR10"

# ## for DMRs
# cat revertant_DMRs_TAIR10.bed stable_DMRs_TAIR10.bed > all_DMRs_TAIR10.bed
# bedtools shuffle -noOverlapping -i revertant_DMRs_TAIR10.bed -g ${ref_dir}/chrom.sizes -excl all_DMRs_TAIR10.bed > shuffled_regions.bed
# nb=$(wc -l revertant_DMRs_TAIR10.bed | awk '{print $1}')
# shuf -n ${nb} ChIP/tracks/TAIR10_protein_coding_genes.bed > random_genes.bed
# dmrnb1=$(wc -l stable_DMRs_TAIR10.bed | awk '{print $1}')
# dmrnb2=$(wc -l revertant_DMRs_TAIR10.bed | awk '{print $1}')
# dmrnb3=$(wc -l shuffled_regions.bed | awk '{print $1}')
# dmrnb4=$(wc -l random_genes.bed | awk '{print $1}')
# dmrsfiles="stable_DMRs_TAIR10.bed revertant_DMRs_TAIR10.bed shuffled_regions.bed random_genes.bed"
# dmrsregionlab="StableDMRs(${dmrnb1}) RevertantDMRs(${dmrnb2}) Random_Regions(${dmrnb3}) Random_Genes(${dmrnb4})"

# # ## For chr5 TEs
# awk '$1=="Chr5"' combined/TSS/TAIR10_all_tes.bed > combined/TSS/TAIR10_chr5_tes.bed
# chrnb1=$(wc -l combined/TSS/TAIR10_chr5_tes.bed | awk '{print $1}')
# chrfile="combined/TSS/TAIR10_chr5_tes.bed"
# chrregionlab="Chr5_TEs(${chrnb1})"

# for type in DMRs chr5_TEs
# do
	# case "${type}" in 
		# DMRs)	tmpnameout="DMRs_and_controls"
				# files="${dmrsfiles}"
				# regionlab="${dmrsregionlab}"
				# colorlist="#3342FF #33E0FF #FFC300 #C70039";;
		# chr5_TEs)	tmpnameout="Chr5_TEs"
					# files="${chrfile}"
					# regionlab="${chrregionlab}"
					# colorlist="#3342FF";;
		# none)	break;;
	# esac
	# param="-bs 10 -a 5000 -b 5000 -m 2000";;
	# samplelist="${samplelist1[*]}"
	# labellist="${label_list1[*]}"
	# orderedmark="${ordered_mark1[*]}"
	# nameout="${tmpnameout}_${output1}_${version}"
	# nbmeta="${#samplelist1[@]}";;
	# printf "\nComputing ${nameout} matrix\n"
	# printf "\nsamples: ${samplelist}\n"
	# printf "\nlabels: ${labellist}\nnb: ${nbmeta}\n"
	# printf "\nordered marks: ${orderedmark}\n"
	# printf "\nfiles: ${files}\n"
	# computeMatrix scale-regions --missingDataAsZero --skipZeros -R ${files} -S ${samplelist} ${param} -p ${threads} -o figures_manuscript/data/matrix_${nameout}.gz
	# printf "\nGetting scale values\n"
	# plotProfile -m figures_manuscript/data/matrix_${nameout}.gz -out figures_manuscript/draft/${nameout}_temp_profile.pdf --samplesLabel ${labellist} --averageType mean --outFileNameData figures_manuscript/data/values_y_${nameout}.txt
	# rm -f figures_manuscript/draft/${nameout}_temp_profile.pdf
	# computeMatrixOperations dataRange -m figures_manuscript/data/matrix_${nameout}.gz > figures_manuscript/data/values_z_${nameout}.txt
	# zmins1=()
	# zmaxs1=()
	# ymins1=()
	# ymaxs1=()
	# for mark in ${orderedmark}
	# do
		# ## To get heatmap scale
		# zmini1=$(grep "${mark}" figures_manuscript/data/values_z_${nameout}.txt | awk 'BEGIN {m=999999} {a=$5; if (a<m) m=a;} END {print m}')
		# zmaxi1=$(grep "${mark}" figures_manuscript/data/values_z_${nameout}.txt | awk 'BEGIN {m=-999999} {a=$6; if (a>m) m=a;} END {print m}')
		# num=$(grep "${mark}" figures_manuscript/data/values_z_${nameout}.txt | wc -l)
		# test=$(awk -v a=${zmini1} -v b=${zmaxi1} 'BEGIN {if (a==0 && b==0) c="yes"; else c="no"; print c}')
		# if [[ ${test} == "yes" ]]; then
			# zmini1="0"
			# zmaxi1="0.005"
		# fi
		# for i in $(seq 1 ${num})
		# do
			# zmins1+=("${zmini1}")
			# zmaxs1+=("${zmaxi1}")
		# done
		
		# ## To get metaplot mean scale
		# ymini1=$(grep "${mark}" figures_manuscript/data/values_y_${nameout}.txt | awk '{m=$3; for(i=3;i<=NF;i++) if ($i<m) m=$i; print m}' | awk 'BEGIN {m=99999} {if ($1<m) m=$1} END {if (m<0) a=m*1.2; else a=m*0.7; print a}')
		# ymaxi1=$(grep "${mark}" figures_manuscript/data/values_y_${nameout}.txt | awk '{m=$3; for(i=3;i<=NF;i++) if ($i>m) m=$i; print m}' | awk 'BEGIN {m=-99999} {if ($1>m) m=$1} END {if (m<0) a=m*0.8; else a=m*1.2; print a}')
		# num=$(grep "${mark}" figures_manuscript/data/values_y_${nameout}.txt | cut -f1 | uniq | wc -l)
		# test=$(awk -v a=${ymini1} -v b=${ymaxi1} 'BEGIN {if (a==0 && b==0) c="yes"; else c="no"; print c}')
		# if [[ ${test} == "yes" ]]; then
			# ymini1="0"
			# ymaxi1="0.01"
		# fi
		# for i in $(seq 1 ${num})
		# do
			# ymins1+=("${ymini1}")
			# ymaxs1+=("${ymaxi1}")
		# done
	# done
	# printf "\nPlotting heatmap and profile for ${nameout} scaled\n"
	# plotHeatmap -m figures_manuscript/data/matrix_${nameout}.gz -out figures_manuscript/draft_fin/${nameout}_heatmap_keep.pdf --sortRegions keep --samplesLabel ${labellist} --regionsLabel ${regionlab} --interpolationMethod 'bilinear' --zMin ${zmins1[@]} --zMax ${zmaxs1[@]} --colorMap "RdYlBu_r" --whatToShow "heatmap and colorbar"			
	# plotHeatmap -m figures_manuscript/data/matrix_${nameout}.gz -out figures_manuscript/draft_fin/${nameout}_heatmap_sorted.pdf --sortRegions descend --sortUsing mean --sortUsingSamples '1' --samplesLabel ${labellist} --regionsLabel ${regionlab} --interpolationMethod 'bilinear' --yMin ${ymins1[@]} --yMax ${ymaxs1[@]} --zMin ${zmins1[@]} --zMax ${zmaxs1[@]} --colorMap "RdYlBu_r" --averageType mean
	# plotHeatmap -m figures_manuscript/data/matrix_${nameout}.gz -out figures_manuscript/draft_fin/${nameout}_heatmap.pdf --sortRegions descend --sortUsing mean --sortUsingSamples '1' --samplesLabel ${labellist} --regionsLabel ${regionlab} --interpolationMethod 'bilinear' --yMin ${ymins1[@]} --yMax ${ymaxs1[@]} --zMin ${zmins1[@]} --zMax ${zmaxs1[@]} --colorMap "RdYlBu_r" --whatToShow "heatmap and colorbar"
	# plotProfile -m figures_manuscript/data/matrix_${nameout}.gz -out figures_manuscript/draft_fin/${nameout}_profile.pdf --plotType 'lines' --averageType 'mean' --regionsLabel ${regionlab} --yMin ${ymins1[@]} --yMax ${ymaxs1[@]} --averageType mean --colors ${colorlist} --numPlotsPerRow ${nbmeta}
# done


#############################################################################
############### Figure1 - Coverage and squiggle plots #######################
#############################################################################

# ID="allchrs"
# ref_dir="/grid/martienssen/home/jcahn/nlsas/Genomes/Arabidopsis/TAIR10"
# rm -f figures_manuscript/data/${ID}_samplefile.txt

# for mark in DDM1 H3K27me1 H3.3 H4K16ac MGH3
# do
	# pathtobw="./ChIP/tracks"
	# if [[ "${mark}" == "MGH3" ]]; then
		# name1="Col0_WT_${mark}_Rep1"
		# name2="ddm1_vs_WT_${mark}_Rep1"
	# else
		# name1="Col0_WT_${mark}_merged"
		# name2="ddm1_vs_WT_${mark}"
	# fi
	# bw1="${pathtobw}/${name1}.bw"
	# bw2="${pathtobw}/${name2}_cpm.bw"
	# printf "WT\t${mark}\t${mark}\t${bw1}\n" >> figures_manuscript/data/${ID}_samplefile.txt
	# printf "ddm1_vs_WT\t${mark}_cpm\tddm1_vs_WT_${mark}_cpm\t${bw2}\n" >> figures_manuscript/data/${ID}_samplefile.txt
# done

# filelist=()
# labellist=()
# while read sample group label bw
# do
	# name="${sample}_${group}"
	# filelist+=("${bw}")
	# labellist+=("${name}")
	# printf "${name} added to filelist\n"
# done < figures_manuscript/data/${ID}_samplefile.txt

# printf "Summarize bigwigs in binsize of 50 kb\n"
# multiBigwigSummary bins -b ${filelist[@]} -l ${labellist[@]} -p ${threads} --chromosomesToSkip "ChrC ChrM" -bs=50000 -out figures_manuscript/data/${ID}.npz --outRawCounts figures_manuscript/data/locus_${ID}.tab
# rm -f figures_manuscript/data/${ID}.npz

# printf "Making squigly plots with R\n"
# sed 's/#//g' figures_manuscript/data/locus_allchrs.tab | sed "s/'//g" > figures_manuscript/data/locus_allchrs.txt

# printf "Summarize bigwigs in binsize of 150 kb\n"
# multiBigwigSummary bins -b ${filelist[@]} -l ${labellist[@]} -p ${threads} --chromosomesToSkip "ChrC ChrM" -bs=150000 -out figures_manuscript/data/${ID}.npz --outRawCounts figures_manuscript/data/locus_${ID}_150kb.tab
# rm -f figures_manuscript/data/${ID}.npz

# printf "Making squigly plots with R\n"
# sed 's/#//g' figures_manuscript/data/locus_allchrs_150kb.tab | sed "s/'//g" > figures_manuscript/data/locus_allchrs_150kb.txt

# qsub -N DDM1 -V -cwd -sync y -pe threads 1 -l m_mem_free=1G -l tmp_free=5G -j y -o logDDM1.log <<-'EOF'
	# #!/usr/bin/env Rscript
	
	# library(readr)
	# library(ggplot2)
	# library(dplyr)
	# library(tidyr)
	# library(ggpubr)
	# library(data.table)

	# tot1<-read.delim("figures_manuscript/data/locus_allchrs.txt", header=TRUE) %>%
		# gather("Mark","Coverage",-chr,-start,-end)
	# tot1$Mark<-as.factor(tot1$Mark)
	
	# tot2<-read.delim("figures_manuscript/data/locus_allchrs_150kb.txt", header=TRUE) %>%
		# gather("Mark","Coverage",-chr,-start,-end)
	# tot2$Mark<-as.factor(tot2$Mark)
		
	# Plot.coverage.chr5<-function(table) {
		# subset1<-c("WT_DDM1","WT_H3K27me1","WT_H3.3","WT_MGH3")
		# tab1<-filter(table, Mark %in% subset1, chr %in% c("Chr5"))
		# tab1$Mark<-gsub("WT_","",as.character(tab1$Mark))
		# tab1<-rowwise(tab1) %>%
			# mutate(comb=ifelse(Coverage<0, paste0(Mark,"_neg"),paste0(Mark,"_pos")))
		
		# tab1$Mark<-factor(tab1$Mark, levels=c("H3.3","MGH3","H3K27me1","DDM1"))

		# plot1<-ggplot(tab1, aes(start,Coverage)) + geom_bar(stat="identity",aes(fill=comb)) +
				# facet_grid(rows=vars(Mark), scales="free_y", switch = "y") +
				# scale_fill_manual(values=c(DDM1_pos="blue",DDM1_neg="grey",H3.3_pos="red",H3.3_neg="grey",H3K27me1_pos="darkgreen",H3K27me1_neg="grey",MGH3_pos="orange",MGH3_neg="grey"),
						# breaks=c("H3.3_pos","MGH3_pos","H3K27me1_pos","DDM1_pos"),
						# labels=c("HTR5-RFP","MGH3","H3K27me1","DDM1")) +
				# labs(y="log2 (IP/Input)",title="Chromosome 5",fill="") +
				# geom_hline(yintercept=0, size=0.2, linetype="dashed") +
				# theme(axis.title.y=element_text(size=10), 
					# axis.title.x=element_blank(),
					# axis.text.x=element_blank(),
					# axis.ticks.x=element_blank(), panel.border = element_rect(fill=NA),
					# panel.spacing = unit(0,'lines'),
					# panel.background = element_rect(fill = "white", colour = "black"),
					# strip.background = element_blank(),
					# strip.text = element_blank(),
					# legend.key=element_blank())
		# print(plot1)
	# }
			
	# Plot.squiggle.allchrs.withDDM1<-function(table) {	
		# subset2<-c("ddm1_vs_WT_H3K27me1_cpm","ddm1_vs_WT_H3.3_cpm","ddm1_vs_WT_MGH3_cpm","ddm1_vs_WT_DDM1_cpm")
		# tab2<-filter(table, Mark %in% subset2, chr %in% c("Chr1","Chr2","Chr3","Chr4","Chr5"))
		# tab2$Mark<-gsub("ddm1_vs_WT_","",as.character(tab2$Mark))
		# tab2$Mark<-gsub("_cpm","",as.character(tab2$Mark))
		# tab2$Mark<-factor(tab2$Mark, levels=c("H3.3","MGH3","H3K27me1","DDM1"))
	
		# chromsize<-data.table(chr=c("Chr1","Chr2","Chr3","Chr4","Chr5"),Centro=c(15086046,3607930,13799418,3956022,11725025))
		# tab3<-merge(tab2,chromsize,by=c("chr"))
  
		# plot<-ggplot(tab3) + geom_line(aes(x=start, y=Coverage, color=Mark),size=0.5) +
				# facet_grid(~chr, scales="free_x", space="free_x", switch="x") +
				# scale_color_manual(values=c(H3.3="red",H3K27me1="darkgreen",MGH3="orange",DDM1="blue")) +
				# labs(y="log2 (ddm1 vs WT)",color="") +
				# expand_limits(y=-1.3) +
				# geom_hline(yintercept=0, size=0.2, linetype="dashed") +
				# geom_pointrange(aes(xmin=start,xmax=end,x=Centro),y=-1.2,size=1) +
				# theme(axis.title.y=element_text(size=10), 
					# axis.title.x=element_blank(),
					# axis.text.x=element_blank(),
					# axis.ticks.x=element_blank(), panel.border = element_rect(fill=NA),
					# panel.spacing = unit(0,'lines'),
					# panel.background = element_rect(fill = "white", colour = "black"),
					# strip.background = element_rect(fill = 'white', colour = 'black'),
					# legend.key=element_blank())
		# print(plot)
	# }
		
	# plot1<-ggarrange(Plot.coverage.chr5(tot1), Plot.squiggle.allchrs.withDDM1(tot2), labels=c("B","C"), common.legend = T, legend = "right", nrow = 2, heights=c(3,2))	
	# pdf("figures_manuscript/draft_fin/Figure_1_genome_50_150kb_withDDM1.pdf",width=10,height=6)
	# print(plot1)
	# dev.off()
	
# EOF

#############################################################################
################### To prepare DNA methylation data #########################
#############################################################################

# cd /grid/martienssen/home/jcahn/nlsas/projects/mC_DDM1

# printf "mC_A\tWT\tRep1\t305054.1\tFCHW5JFAFXY\tNo\tPE\nmC_A\thira\tRep1\t305054.2\tFCHW5JFAFXY\tNo\tPE\nmC_A\tddm1\tRep1\t305054.3\tFCHW5JFAFXY\tNo\tPE\nmC_A\tddm1hira\tRep1\t305054.4\tFCHW5JFAFXY\tNo\tPE\nmC_A\tWT\tRep2\tWT_S1\tFCAAC3TYMM5\tNo\tPE\nmC_A\thira\tRep2\thira_S2\tFCAAC3TYMM5\tNo\tPE\nmC_A\tddm1\tRep2\tddm1_S3\tFCAAC3TYMM5\tNo\tPE\nmC_A\tddm1hira\tRep2\tddm1hira_S4\tFCAAC3TYMM5\tNo\tPE\nmC_B\tWT\tRep1\tWT1_S5\tFCAAC3TYMM5\tNo\tPE\nmC_B\tatrx\tRep1\tatrx1_S6\tFCAAC3TYMM5\tNo\tPE\nmC_B\tddm1\tRep1\tddm11_S7\tFCAAC3TYMM5\tNo\tPE\nmC_B\tddm1atrx\tRep1\tddm1atrx1_S8\tFCAAC3TYMM5\tNo\tPE\nmC_B\tWT\tRep2\tWT2_S9\tFCAAC3TYMM5\tNo\tPE\nmC_B\tatrx\tRep2\tatrx2_S10\tFCAAC3TYMM5\tNo\tPE\nmC_B\tddm1\tRep2\tddm12_S11\tFCAAC3TYMM5\tNo\tPE\nmC_B\tddm1atrx\tRep2\tddm1atrx2_S12\tFCAAC3TYMM5\tNo\tPE\n" > DDM1_samplefile.txt

# qsub -sync y -N met_DDM1mC -o mC_DDM1.log /grid/martienssen/home/jcahn/data/Scripts/MethylC_seq_primary.sh TAIR10 DDM1_samplefile.txt No &
# pid="$!"
# wait ${pid}

#############################################################################
####################### To make all bigwig files ############################
#############################################################################

# ref_dir="/grid/martienssen/home/jcahn/nlsas/Genomes/Arabidopsis/Col0"

# for sample in WT_A ddm1_A hira_A ddm1hira_A WT_B ddm1_B atrx_B ddm1atrx_B
# do
	# printf "${sample}\n"
	# zcat methylcall/${sample}_Rep*.deduplicated.CX_report.txt.gz | sort -k1,1 -k2,2n | awk -v OFS="\t" '{print $1,$2-1,$2,$3,$4,$5,$6,$7}' > methylcall/temp_${sample}.bed
	# bedtools merge -d -1 -o distinct,sum,sum,distinct,distinct -c 4,5,6,7,8 -i methylcall/temp_${sample}.bed > methylcall/temp2_${sample}.bed
	# cat methylcall/temp2_${sample}.bed | awk -v OFS="\t" -v s=${sample} '($5+$6)>0 {a=$5+$6; if ($7=="CHH") print $1,$2,$3,$5/a*100 > "methylcall/"s"_CHH.bedGraph"; else if ($7=="CHG") print $1,$2,$3,$5/a*100 > "methylcall/"s"_CHG.bedGraph"; else if ($7=="CG") print $1,$2,$3,$5/a*100 > "methylcall/"s"_CG.bedGraph"}'
	# # cat methylcall/temp2_${sample}.bed | awk -v OFS="\t" -v s=${sample} '($5+$6)>0 {a=$5+$6; if ($8=="CAG" || $8=="CTG") print $1,$2,$3,$5/a*100 > "methylcall/"s"_CWG.bedGraph"; else if ($8~"CAA" || $8~"CTA") print $1,$2,$3,$5/a*100 > "methylcall/"s"_CWA.bedGraph"}'
	
	# rm -f methylcall/temp*
	# for context in CG CHG CHH CWG CWA
	# do
		# printf "\nMaking bigwig files of ${context} context for ${sample}\n"
		# LC_COLLATE=C sort -k1,1 -k2,2n methylcall/${sample}_${context}.bedGraph > methylcall/sorted_${sample}_${context}.bedGraph
		# bedGraphToBigWig methylcall/sorted_${sample}_${context}.bedGraph ${ref_dir}/chrom.sizes methylcall/${sample}_${context}.bw
	# done
	# rm -f methylcall/*.bedGraph
		
	### To add a coverage cutoff
	# printf "${sample}\n"
	# zcat methylcall/${sample}_Rep*.deduplicated.CX_report.txt.gz | sort -k1,1 -k2,2n | awk -v OFS="\t" '($4+$5)>3 {print $1,$2-1,$2,$3,$4,$5,$6,$7}' > methylcall/temp_${sample}.bed
	# bedtools merge -d -1 -o distinct,sum,sum,distinct,distinct -c 4,5,6,7,8 -i methylcall/temp_${sample}.bed > methylcall/temp2_${sample}.bed
	# cat methylcall/temp2_${sample}.bed | awk -v OFS="\t" -v s=${sample} '($5+$6)>6 {a=$5+$6; if ($7=="CHH") print $1,$2,$3,$5/a*100 > "methylcall/"s"_CHH.bedGraph"; else if ($7=="CHG") print $1,$2,$3,$5/a*100 > "methylcall/"s"_CHG.bedGraph"; else if ($7=="CG") print $1,$2,$3,$5/a*100 > "methylcall/"s"_CG.bedGraph"}'
	# # cat methylcall/temp2_${sample}.bed | awk -v OFS="\t" -v s=${sample} '($5+$6)>0 {a=$5+$6; if ($8=="CAG" || $8=="CTG") print $1,$2,$3,$5/a*100 > "methylcall/"s"_CWG.bedGraph"; else if ($8~"CAA" || $8~"CTA") print $1,$2,$3,$5/a*100 > "methylcall/"s"_CWA.bedGraph"}'
	
	# rm -f methylcall/temp*
	# for context in CG CHG CHH
	# do
		# printf "\nMaking bigwig files of ${context} context for ${sample}\n"
		# LC_COLLATE=C sort -k1,1 -k2,2n methylcall/${sample}_${context}.bedGraph > methylcall/sorted_${sample}_${context}.bedGraph
		# bedGraphToBigWig methylcall/sorted_${sample}_${context}.bedGraph ${ref_dir}/chrom.sizes methylcall/${sample}_${context}_min.bw
	# done
	# rm -f methylcall/*.bedGraph
	
	# ### To make stranded bw
		
	# printf "${sample}\n"
	# zcat methylcall/${sample}_Rep*.deduplicated.CX_report.txt.gz | sort -k1,1 -k2,2n | awk -v OFS="\t" '{print $1,$2-1,$2,$3,$4,$5,$6,$7}' > methylcall/temp_${sample}.bed
	# bedtools merge -d -1 -o distinct,sum,sum,distinct,distinct -c 4,5,6,7,8 -i methylcall/temp_${sample}.bed > methylcall/temp2_${sample}.bed
	# for strand in plus minus
	# do
		# case "${strand}" in 
			# plus)	sign="+";;
			# minus)	sign="-";;
		# esac
		# awk -v n=${sign} '$4==n' methylcall/temp2_${sample}.bed | awk -v OFS="\t" -v s=${sample} -v d=${strand} '($5+$6)>0 {a=$5+$6; if ($7=="CHH") print $1,$2,$3,$5/a*100 > "methylcall/"s"_CHH_"d".bedGraph"; else if ($7=="CHG") print $1,$2,$3,$5/a*100 > "methylcall/"s"_CHG_"d".bedGraph"; else if ($7=="CG") print $1,$2,$3,$5/a*100 > "methylcall/"s"_CG_"d".bedGraph"}'
		# for context in CG CHG CHH
		# do
			# printf "\nMaking bigwig files of ${context} context ${strand} strand for ${sample}\n"
			# LC_COLLATE=C sort -k1,1 -k2,2n methylcall/${sample}_${context}_${strand}.bedGraph > methylcall/sorted_${sample}_${context}_${strand}.bedGraph
			# bedGraphToBigWig methylcall/sorted_${sample}_${context}_${strand}.bedGraph ${ref_dir}/chrom.sizes methylcall/${sample}_${context}_${strand}.bw
		# done
	# done
	
	# ### To make one file per replicate
	
	# for rep in Rep1 Rep2
	# do
		# printf "${sample} ${rep}\n"
		# zcat methylcall/${sample}_${rep}.deduplicated.CX_report.txt.gz | awk -v OFS="\t" '($4+$5)>0 {a=$4+$5; print $1,$2-1,$2,$4/a*100}' > methylcall/temp_${sample}_${rep}.bedGraph
		# zcat methylcall/${sample}_${rep}.deduplicated.CX_report.txt.gz | awk -v OFS="\t" '($4+$5)>3 {a=$4+$5; print $1,$2-1,$2,$4/a*100}' > methylcall/temp_${sample}_${rep}_min3.bedGraph
		# printf "\nMaking bigwig file for ${sample}\n"
		# LC_COLLATE=C sort -k1,1 -k2,2n methylcall/temp_${sample}_${rep}.bedGraph > methylcall/sorted_${sample}_${rep}.bedGraph
		# bedGraphToBigWig methylcall/sorted_${sample}_${rep}.bedGraph ${ref_dir}/chrom.sizes methylcall/${sample}_${rep}_CX.bw
		# LC_COLLATE=C sort -k1,1 -k2,2n methylcall/temp_${sample}_${rep}_min3.bedGraph > methylcall/sorted_${sample}_${rep}_min3.bedGraph
		# bedGraphToBigWig methylcall/sorted_${sample}_${rep}_min3.bedGraph ${ref_dir}/chrom.sizes methylcall/${sample}_${rep}_CX_min3.bw
	# done
	# rm -f methylcall/temp*
	# rm -f methylcall/*.bedGraph
# done

# cd /grid/martienssen/home/jcahn/norepl/projects/manuscript_DDM1/analysis

#############################################################################
######################## Figure2 - mC Metaplots #############################
#############################################################################

# pathtobw="/grid/martienssen/home/jcahn/nlsas/projects/mC_DDM1/methylcall"

# # for type in A_Rep1 B_Rep2
# for type in A B
# do
	# case "${type}" in
		# A*)	uniqsample=("WT_${type}" "ddm1_${type}" "hira_${type}" "ddm1hira_${type}")
			# uniqlabel=("WT_${type}" "ddm1_${type}" "hira_${type}" "ddm1hira_${type}")
			# colors=("blue" "red" "orange" "grey");;
		# B*)	uniqsample=("WT_${type}" "ddm1_${type}" "atrx_${type}" "ddm1atrx_${type}")
			# uniqlabel=("WT_${type}" "ddm1_${type}" "atrx_${type}" "ddm1atrx_${type}")
			# colors=("blue" "red" "orange" "grey");;
	# esac
	# for context in CG CHG CHH
	# do
		# samplelist1=()
		# samplelist2=()
		# for sample in ${uniqsample[*]}
		# do
			# samplelist1+=("${pathtobw}/${sample}_${context}.bw")
			# samplelist2+=("${pathtobw}/${sample}_${context}_min.bw")
		# done
		# filename="TEs_${context}_${type}"
		# file="combined/TSS/TAIR10_all_tes.bed"
		# title="TEs"
		# nb=$(wc -l ${file} | awk '{print $1}')
		# regionlab="${title} (${nb})"
		# param="-bs 50 -a 1000 -b 1000 -m 5000"
		
		# printf "\nComputing matrix for ${filename}\n"
		# computeMatrix scale-regions -R ${file} -S ${samplelist1[@]} ${param} -p ${threads} -o figures_manuscript/data/matrix_${filename}.gz --smartLabels
		# printf "\nPlotting profile for ${filename}\n"
		# plotProfile -m figures_manuscript/data/matrix_${filename}.gz -out figures_manuscript/draft_fin/profile_${filename}.pdf --plotType 'lines' --averageType 'mean' --perGroup --colors ${colors[@]} --regionsLabel "${regionlab}"
			
		# printf "\nComputing matrix for ${filename} mincov\n"
		# computeMatrix scale-regions -R ${file} -S ${samplelist2[@]} ${param} -p ${threads} -o figures_manuscript/data/matrix_${filename}_mincov.gz --smartLabels
		# printf "\nPlotting profile for ${filename} mincov\n"
		# plotProfile -m figures_manuscript/data/matrix_${filename}_mincov.gz -out figures_manuscript/draft_fin/profile_${filename}_mincov.pdf --plotType 'lines' --averageType 'mean' --perGroup --colors ${colors[@]} --regionsLabel "${regionlab}"
	# done
# done

#############################################################################
############ Figure S10 - to plot correlation of replicates #################
#############################################################################

# pathtobw="/grid/martienssen/home/jcahn/nlsas/projects/mC_DDM1/methylcall"

# uniqsample=("WT_A" "ddm1_A" "hira_A" "ddm1hira_A" "WT_B" "ddm1_B" "atrx_B" "ddm1atrx_B")

# # for context in CG CHG CHH
# for context in CX CX_min3
# do
	# for sample in ${uniqsample[*]}
	# do
		# printf "${sample} ${context}\n"
		# label=${sample%%_*}
		# multiBigwigSummary bins -b ${pathtobw}/${sample}_Rep1_${context}.bw ${pathtobw}/${sample}_Rep2_${context}.bw -o figures_manuscript/data/${sample}_${context}.npz --labels Rep1 Rep2 -p ${threads} --chromosomesToSkip "ChrC ChrM"
		# plotCorrelation --corData figures_manuscript/data/${sample}_${context}.npz --corMethod pearson --whatToPlot scatterplot --plotFile figures_manuscript/draft_fin/pearson_${sample}_${context}.pdf --plotTitle "${label} ${context}"
		# printf "\n\n"
	# done
# done

# pathchip="/grid/martienssen/home/jcahn/norepl/projects/manuscript_DDM1/analysis/ChIP/tracks"

# while read type ref sample mark paired ref_dir
# do
	# name="${sample}_${mark}"
	# if [ -e ${pathchip}/${ref}_${name}_Rep1.bw ] && [ -e ${pathchip}/${ref}_${name}_Rep2.bw ]; then
		# printf "${name}\n"
		# multiBigwigSummary bins -b ${pathchip}/${ref}_${name}_Rep1.bw ${pathchip}/${ref}_${name}_Rep2.bw -o figures_manuscript/data/${name}.npz --labels Rep1 Rep2 -p ${threads} --chromosomesToSkip "ChrC ChrM"
		# plotCorrelation --corData figures_manuscript/data/${name}.npz --corMethod pearson --whatToPlot scatterplot --plotFile figures_manuscript/draft_fin/pearson_${name}_noout.pdf --plotTitle "${sample} ${mark}" --removeOutliers
	# fi
# done < subset_analysis_samplefile.txt

# printf "WT MGH3 Borg\n"
# multiBigwigSummary bins -b ${pathchip}/Col0_WT_MGH3_Rep1.bw ${pathchip}/Col0_WT_borg_MGH3_Rep1.bw -o figures_manuscript/data/WT_MGH3_1.npz --labels "WT MGH3 (this study)" "MGH3 rep1 (Borg et al. 2020)" -p ${threads} --chromosomesToSkip "ChrC ChrM"
# plotCorrelation --corData figures_manuscript/data/WT_MGH3_1.npz --corMethod pearson --whatToPlot scatterplot --plotFile figures_manuscript/draft_fin/pearson_open_borg_WT_MGH3_rep1.pdf --plotTitle "WT MGH3 correlation"
# printf "ddm1 MGH3 Borg\n"
# multiBigwigSummary bins -b ${pathchip}/Col0_ddm1_MGH3_Rep1.bw ${pathchip}/Col0_WT_borg_MGH3_Rep1.bw -o figures_manuscript/data/ddm1_MGH3_1.npz --labels "ddm1 MGH3 (this study)" "MGH3 rep1 (Borg et al. 2021)" -p ${threads} --chromosomesToSkip "ChrC ChrM"
# plotCorrelation --corData figures_manuscript/data/ddm1_MGH3_1.npz --corMethod pearson --whatToPlot scatterplot --plotFile figures_manuscript/draft_fin/pearson_open_borg_ddm1_MGH3_rep1.pdf --plotTitle "ddm1 MGH3 correlation"

# printf "WT H3.3 Stroud\n"
# multiBigwigSummary bins -b ${pathchip}/Col0_WT_H3.3_merged.bw ${pathchip}/Col0_WT_Stroud_H3.3_Rep1.bw -o figures_manuscript/data/WT_H3.3.npz --labels "H3.3 (this study)" "H3.3 (Stroud et al. 2012)" -p ${threads} --chromosomesToSkip "ChrC ChrM"
# plotCorrelation --corData figures_manuscript/data/WT_H3.3.npz --corMethod pearson --whatToPlot scatterplot --plotFile figures_manuscript/draft_fin/pearson_open_stroud_H3.3.pdf --plotTitle "WT H3.3 correlation"
# printf "WT H3.1 Stroud\n"
# multiBigwigSummary bins -b ${pathchip}/Col0_WT_H3K27me1_merged.bw ${pathchip}/Col0_WT_Stroud_H3.1_Rep1.bw -o figures_manuscript/data/WT_H3.1.npz --labels "H3K27me1 (this study)" "H3.1 (Stroud et al. 2012)" -p ${threads} --chromosomesToSkip "ChrC ChrM"
# plotCorrelation --corData figures_manuscript/data/WT_H3.1.npz --corMethod pearson --whatToPlot scatterplot --plotFile figures_manuscript/draft_fin/pearson_open_stroud_H3.1.pdf --plotTitle "WT H3K27me1/H3.1"

# rm -f figures_manuscript/data/*.npz

#############################################################################
################ To call DMRs between ddm1hira and ddm1 #####################
#############################################################################

# qsub -N DMRs -V -cwd -sync y -pe threads 10 -l m_mem_free=2G -l tmp_free=4G -j y -o logDDM1.log <<-'EOF'
	# #!/usr/bin/env Rscript

	# library(DMRcaller)
	# library(dplyr)

	# tair10_chrs<-GRanges(Rle(c("Chr1","Chr2","Chr3","Chr4","Chr5")), ranges=IRanges(start=c(1,1,1,1,1), end=c(30427671,19698289,23459830,18585056,26975502)))

	# ddm1hirarep1<-readBismarkPool("DMRs/CX_reports/ddm1hira_A_Rep1.deduplicated.CX_report.txt.gz")
	# ddm1hirarep2<-readBismarkPool("DMRs/CX_reports/ddm1hira_A_Rep2.deduplicated.CX_report.txt.gz")
	# ddm1rep1<-readBismark("DMRs/CX_reports/ddm1_A_Rep1.deduplicated.CX_report.txt.gz")
	# ddm1rep2<-readBismark("DMRs/CX_reports/ddm1_A_Rep2.deduplicated.CX_report.txt.gz")
	# ddm1hirapool<-poolTwoMethylationDatasets(ddm1hirarep1,ddm1hirarep2)
	# ddm1pool<-poolTwoMethylationDatasets(ddm1rep1,ddm1rep2)
	
	# ###
	
	# call.DMR.bin<-function(seqcontext,sample1,sample2,lab) {
		# if (seqcontext=="CG") {
			# min=0.3
		# } else if (seqcontext=="CHG") {
			# min=0.2
		# } else {
			# min=0.1
		# }
		# DMRs<-computeDMRs(sample1, sample2, regions=tair10_chrs, context=seqcontext, method="bins", binSize=100, test="score", pValueThreshold=0.01, minCytosinesCount=5, minProportionDifference=min, minGap=200, minSize=50, minReadsPerCytosine=3, cores=10)
		# DMRsmerged<-mergeDMRsIteratively(DMRs, minGap=200, respectSigns=TRUE, sample1, sample2, context=seqcontext, minProportionDifference=min, minReadsPerCytosine=3, pValueThreshold=0.01, test="score")	
		# tab<-data.frame(Chr=seqnames(DMRsmerged),Start=start(DMRsmerged)-1,End=end(DMRsmerged),ddm1hira=elementMetadata(DMRsmerged)[,3],ddm1=elementMetadata(DMRsmerged)[,6], Pvalue=elementMetadata(DMRsmerged)[,10]) %>%
			# mutate(Delta=ddm1hira-ddm1) %>%
			# arrange(desc(Delta),Pvalue)
		# tab1<-unique(tab)
		# tab2<-format(tab, scientific = FALSE)
		# write.table(tab2,paste0("DMRs/DMRs_",seqcontext,"_",lab,"_bins.txt"),sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
	# }
		
	# for (context in c("CG","CHG","CHH")) {
		# call.DMR.bin(context,ddm1hirapool,ddm1pool,"pool")
	# }

# EOF

#############################################################################
################ Figure2 - Violin plots H3.3 at DMRs ########################
#############################################################################


# bwfiles=()
# labels=()
# for sample in WT ddm1
# do
	# bwfiles+=("ChIP/tracks/Col0_${sample}_H3.3_merged.bw")
	# labels+=("${sample}_H3.3_vs_Input")
# done

# rm -f DMRs/all_DMRs_bins.bed
# for context in CG CHG CHH
# do	
	# ### on DMRs
	# awk '$7>0' DMRs/DMRs_${context}_pool_bins.txt | sort -k6,6n > DMRs/filtered_DMRs_${context}_pool_bins.txt	
	# awk -v OFS="\t" -v c=${context} '{print $1,$2,$3,"DMRs_"c"_"NR,".",".",c}' DMRs/filtered_DMRs_${context}_pool_bins.txt >> DMRs/all_DMRs_bins.bed
# done
# sort -k1,1 -k2,2n DMRs/all_DMRs_bins.bed > DMRs/sorted_DMRs_bins.bed
# bedtools merge -i DMRs/sorted_DMRs_bins.bed -o distinct -c 7 > DMRs/merged_DMRs_bins.bed

# cat DMRs/merged_DMRs_bins.bed > DMRs/regions_for_violins.bed

# # awk -v OFS="\t" '$4=="CG,CHG,CHH" {print $1,$2,$3}' DMRs/merged_DMRs_bins.bed > DMRs/temp_TE_like.bed
# # bedtools shuffle -noOverlapping -i DMRs/temp_TE_like.bed -g chrom.sizes -excl DMRs/merged_DMRs_bins.bed | awk -v OFS="\t" '{print $1,$2,$3,"Shuffled"}' > DMRs/shuffled_TE_like.bed

# cat DMRs/shuffled_TE_like.bed >> DMRs/regions_for_violins.bed

# multiBigwigSummary BED-file -p ${threads} -b ${bwfiles[@]} --labels ${labels[@]} --BED DMRs/regions_for_violins.bed --outRawCounts DMRs/H3_on_DMRs_and_rand.tab -o DMRs/H3_on_DMRs_and_rand.npz

# printf "Making violin plots with R\n"
# sed 's/#//g' DMRs/H3_on_DMRs_and_rand.tab | sed "s/'//g" > DMRs/H3_on_DMRs_and_rand.txt

# Rscript --vanilla - <<-'EOF'
	# #!/usr/bin/env Rscript
	# library(readr)
	# library(ggplot2)
	# library(dplyr)
	# library(tidyr)
	# library(stringr)
	# library(ggpubr)

	# dmrs<-read.table("DMRs/regions_for_violins.bed", header=FALSE,
                 # col.names = c("chr","start","end","Context"))

	# tab<-read.table("DMRs/H3_on_DMRs_and_rand.txt", header=TRUE) %>%
		# merge(dmrs, by=c("chr","start","end")) %>%
		# mutate(DMR=paste0(chr,"_",start)) %>%
		# select(-chr,-start,-end) %>%
		# gather(Sample,Value,-DMR,-Context)
	# tab$Sample<-as.factor(tab$Sample)
	
	# temp<-read.table("DMRs/H3_on_DMRs_and_rand.txt", header=TRUE) %>%
			# merge(dmrs, by=c("chr","start","end")) %>%
			# rowwise() %>%
			# mutate(Groups=ifelse(Context=="Shuffled","Random",
					# ifelse(Context=="CG","CG",
                       # ifelse(Context=="CHG" | Context=="CHH", "non-CG",
                              # ifelse(Context=="CG,CHG", "DDM1 targets",
                                     # ifelse(Context=="CG,CHG,CHH", "TE-like", "Other")))))) %>%
			# group_by(Groups) %>%
			# filter(Groups!="Other") %>%
			# summarise(Tot=n()) %>%
			# mutate(Label=paste0(Groups," (n=",Tot,")")) %>%
			# select(-Tot)

	# tabfin<-mutate(tab, Groups=ifelse(Context=="Shuffled","Random",
						# ifelse(Context=="CG","CG",
                             # ifelse(Context=="CHG" | Context=="CHH", "non-CG",
                                    # ifelse(Context=="CG,CHG", "DDM1 targets",
                                           # ifelse(Context=="CG,CHG,CHH", "TE-like", "Other")))))) %>%
			# merge(temp, by=c("Groups"))
				
	# tabfin$Groups<-factor(tabfin$Groups, levels=c("CG","non-CG","DDM1 targets","TE-like","Random"))
	# tabfin<-arrange(tabfin, Groups)
	# tabfin$Label<-factor(tabfin$Label, levels=unique(tabfin$Label))

	# tabfin$Sample<-factor(tabfin$Sample, levels=c("WT_H3.3_vs_Input","ddm1_H3.3_vs_Input"),
                      # labels=c("WT","ddm1"))

	# plot1<-ggplot(tabfin, aes(Sample,Value)) +
		# geom_violin(alpha=0.3, aes(color=Sample, fill=Sample),
              # draw_quantiles = c(0.5), linewidth=0.5, show.legend=FALSE) +
		# theme_classic() +
		# scale_color_manual(values = c("WT"="blue","ddm1"="red")) +
		# scale_fill_manual(values = c("WT"="blue","ddm1"="red")) +
		# labs(x="",y="H3.3 log2(IP/Input)") +
		# stat_compare_means(method = "t.test", paired=FALSE, size=2,
                     # label="p.signif",ref.group = "WT") +
		# facet_grid(~Label) +
		# theme(strip.text.x = element_text(size = 6))
	
	# pdf("figures_manuscript/revisions/H3.3_at_DMRs_v1.pdf",width=5,height=3)
	# print(plot1)
	# dev.off()

# EOF

#############################################################################
############################ Figure2 - browser shots ########################
#############################################################################

# rm -f all_DMRs_TAIR10.bed
# awk -v OFS="\t" '{print $1,$2,$3,"revertant",".","."}' revertant_DMRs_TAIR10.bed > all_DMRs_TAIR10.bed
# awk -v OFS="\t" '{print $1,$2,$3,"stable",".","."}' stable_DMRs_TAIR10.bed >> all_DMRs_TAIR10.bed

# # for samp in browser_subset_A browser_subset_B
# for samp in browser_subset_A
# do
	# nameout="${samp}"
	# samplefile="${samp}.txt"
	# rm -f figures_manuscript/data/${nameout}_samplefile.txt
	
	# filelist=()
	# labellist=()
	# while read data line tissue sample paired ref_dir
	# do
		# case "${tissue}" in
			# WT)		backcolor="blue";;
			# ddm1)	backcolor="red";;
			# ddm1hira) 	backcolor="grey";;
			# ddm1atrx) 	backcolor="pink";;
			# hira) 	backcolor="orange";;
			# atrx) backcolor="purple";;
		# esac
		# case "${data}" in
			# ChIP*) 	datatype="ChIP"
					# pathtobw="/grid/martienssen/home/jcahn/norepl/projects/manuscript_DDM1/analysis/ChIP/tracks"
					# name="${line}_${tissue}"
					# label="${tissue}_${sample}";;
			# mC*)	datatype="mC"
				# pathtobw="/grid/martienssen/home/jcahn/nlsas/projects/mC_DDM1/methylcall"
				# group=${data#mC_}
				# name="${tissue}_${group}"
				# context="${sample}"
				# label="${tissue}";;
		# esac
		# ref=${ref_dir##*/}
		# if [[ ${datatype} == ChIP ]]; then
			# case "${sample}" in 
				# H3K9me2)	trackcolor="#EE616E"
							# fillcolor1="#ebd8da"
							# fillcolor2="#ed939b";;
				# DDM1)	trackcolor="#0728ed"
						# fillcolor1="#afb8ed"
						# fillcolor2="#2541e8";;
				# H3.3)	trackcolor="#ff0000"
						# fillcolor1="#ff9999"
						# fillcolor2="#ff0000";;
				# H3K27me1)	trackcolor="#e825c7"
						# fillcolor1="#de9ed3"
						# fillcolor2="#e65ace";;
				# H4K20me1)	trackcolor="#4a2d11"
							# fillcolor1="#cf9b69"
							# fillcolor2="#784d23";;
				# H4K16ac)	trackcolor="#eb7f17"
							# fillcolor1="#edb177"
							# fillcolor2="#ed9037";;
				# H*)	trackcolor="#15b3d6"
							# fillcolor1="#a3cbd4"
							# fillcolor2="#52bad1";;
			# esac
			# if [ -e ${pathtobw}/${name}_${sample}_merged.bw ]; then
				# bw="${pathtobw}/${name}_${sample}_merged.bw"
			# elif [ -e ${pathtobw}/${name}_${sample}_Rep1.bw ]; then
				# bw="${pathtobw}/${name}_${sample}_Rep1.bw"
			# else
				# bw="${pathtobw}/${name}_${sample}.bw"
			# fi
			# printf "${name}\t${sample}\t${label}\t${bw}\t${backcolor}\t${trackcolor}\t${fillcolor1}\t${fillcolor2}\tno\n" >> figures_manuscript/data/${nameout}_samplefile.txt
			# filelist+=("${bw}")
			# labellist+=("${name}_${sample}")
		# elif [[ ${datatype} == mC ]]; then
			# case "${context}" in 
				# CG) trackcolor="#5a5255"
					# fillcolor1="#5a5255"
					# fillcolor2="#5a5255"
					# max=100;;
				# CHG)	trackcolor="#1b85b8"
						# fillcolor1="#1b85b8"
						# fillcolor2="#1b85b8"
						# max=100;;
				# CHH)	trackcolor="#559e83"
						# fillcolor1="#559e83"
						# fillcolor2="#559e83"
						# max=100;;
			# esac
			# if [ -e ${pathtobw}/${name}_${context}.bw ]; then
				# bw1="${pathtobw}/${name}_${context}_plus.bw"
				# bw2="${pathtobw}/${name}_${context}_minus.bw"
			# else
				# bw1="${pathtobw}/${name}_Rep1_${context}_plus.bw"
				# bw2="${pathtobw}/${name}_Rep1_${context}_minus.bw"
			# fi
			# printf "${name}\t${context}_plus\t${name}_${context}_plus\t${bw1}\t${backcolor}\t${trackcolor}\t${fillcolor1}\t${fillcolor2}\tno\t0;${max}\n" >> figures_manuscript/data/${nameout}_samplefile.txt
			# printf "${name}\t${context}_minus\t${name}_${context}_minus\t${bw2}\t${backcolor}\t${trackcolor}\t${fillcolor1}\t${fillcolor2}\tno\t0;-${max}\n" >> figures_manuscript/data/${nameout}_samplefile.txt
			# filelist+=("${bw1} ${bw2}")
			# labellist+=("${name}_${context}_plus ${name}_${context}_minus")
		# fi
	# done < ${samplefile}

	# rm -f figures_manuscript/data/${nameout}_loci.txt
	# printf "Chr5:11503749:11513009\tpanelA_${nameout}_max\t1\n" >> figures_manuscript/data/${nameout}_loci.txt
	# # printf "Chr5:11651769:11659707\tpanelB_${nameout}_max\t1\n" >> figures_manuscript/data/${nameout}_loci.txt
	# # printf "Chr5:11493749:11523009\tpanelA_large_${nameout}_max\t1\n" >> figures_manuscript/data/${nameout}_loci.txt
	# # printf "Chr5:11641769:11669707\tpanelB_large_${nameout}_max\t1\n" >> figures_manuscript/data/${nameout}_loci.txt
	# # awk -v OFS="\t" -v n=${nameout} '$4=="CG,CHG,CHH" {a=$2-1000; b=$3+1000; print $1":"a":"b,"DMR_"NR"_"n,1}' DMRs/merged_DMRs_bins.bed | sort -k5,5nr | awk 'NR<=10' >> figures_manuscript/data/${nameout}_loci.txt

	# ref_dir="/grid/martienssen/home/jcahn/nlsas/Genomes/Arabidopsis/TAIR10"
	
	# while read locus ID bins hstart hwid
	# do
		# printf "Preparing files for browser on ${ID}\n"
		# printf "${locus}\n" | awk -F"[:]" -v OFS="\t" '{print $1,$2,$3}' > figures_manuscript/data/locus_${ID}.bed
		# chr=$(awk '{print $1}' figures_manuscript/data/locus_${ID}.bed)
		# start=$(awk '{print $2}' figures_manuscript/data/locus_${ID}.bed)
		# end=$(awk '{print $3}' figures_manuscript/data/locus_${ID}.bed)
		# region="${locus}"
		# printf "Summarize bigwigs in binsize of ${bins} bp on ${region}\n"
		# multiBigwigSummary bins -b ${filelist[@]} -l ${labellist[@]} -r ${region} -p ${threads} -bs=${bins} -out figures_manuscript/data/${ID}.npz --outRawCounts figures_manuscript/data/locus_${ID}.tab
		# rm -f figures_manuscript/data/${ID}.npz

		# while read sample group label bw backcolor trackcolor fillcolor1 fillcolor2 invert
		# do
			# name="${sample}_${group}"
			# printf "Making bw for ${name}\n"
			# col=($(awk -v ORS=" " -v t=${name} 'NR==1 {for(i=1;i<=NF;i++) if ($i~t) print i}' figures_manuscript/data/locus_${ID}.tab))
			# if [[ ${invert} == "invert" ]]; then
				# awk -v OFS="\t" -v a=${col} 'NR>1 {if ($a == "nan") b=0; else b=-$a; print $1,$2,$3,b}' figures_manuscript/data/locus_${ID}.tab | bedtools sort -g ${ref_dir}/chrom.sizes > figures_manuscript/data/${ID}_${name}.bedGraph
			# else
				# awk -v OFS="\t" -v a=${col} 'NR>1 {if ($a == "nan") b=0; else b=$a; print $1,$2,$3,b}' figures_manuscript/data/locus_${ID}.tab | bedtools sort -g ${ref_dir}/chrom.sizes > figures_manuscript/data/${ID}_${name}.bedGraph
			# fi
			# bedGraphToBigWig figures_manuscript/data/${ID}_${name}.bedGraph ${ref_dir}/chrom.sizes figures_manuscript/data/${name}_locus.bw
		# done < figures_manuscript/data/${nameout}_samplefile.txt
	
		# printf "Name\tGroup\tLabel\tBackcolor\tTrackcolor\tFillcolor1\tFillcolor2\tYmin\tYmax\n" > figures_manuscript/data/filenames_${ID}.txt
		# while read sample group label bw backcolor trackcolor fillcolor1 fillcolor2 invert lims
		# do
			# if [[ ${lims} != "" ]]; then
				# ylimmin=${lims%%;*}
				# ylimmax=${lims##*;}
			# else
				# ylimmin=$(cat figures_manuscript/data/${ID}_*${group}.bedGraph | awk 'BEGIN {a=9999} {if ($4<a) a=$4;} END {if (a<0) b=a*1.2; else b=a*0.8; print b}' )
				# ylimmax=$(cat figures_manuscript/data/${ID}_*${group}.bedGraph | awk 'BEGIN {a=-9999} {if ($4>a) a=$4;} END {if (a>0) b=a*1.2; else b=a*0.8; print b}' )
			# fi
			# printf "${sample}\t${group}\t${label}\t${backcolor}\t${trackcolor}\t${fillcolor1}\t${fillcolor2}\t${ylimmin}\t${ylimmax}\n" >> figures_manuscript/data/filenames_${ID}.txt
		# done < figures_manuscript/data/${nameout}_samplefile.txt
		# rm -f figures_manuscript/data/${ID}_*.bedGraph
	
		# awk -v OFS="\t" '$3 ~ "gene"' ~/nlsas/Genomes/Arabidopsis/TAIR10_GFF3_genes_transposons.gff | awk -F"[=;]" -v OFS="\t" '{print $1,$2}' | awk -v OFS="\t" '{print $1,$4-1,$5,$10,".",$7}' > ChIP/tracks/TAIR10_all_features.bed
		# bedtools intersect -a ChIP/tracks/TAIR10_all_features.bed -b figures_manuscript/data/locus_${ID}.bed | awk -v OFS="\t" '{if ($6!="+" && $6!="-") $6="*"; print $0}' > figures_manuscript/data/features_in_locus_${ID}.bed
		# bedtools intersect -a DMRs/merged_DMRs_bins.bed -b figures_manuscript/data/locus_${ID}.bed | awk -v OFS="\t" '{if ($6!="+" && $6!="-") $6="*"; print $0}' > figures_manuscript/data/DMRs_in_locus_${ID}.bed
		
		# printf "\nPlotting browser on ${ID}\n\n"
		
		# if [[ htsart != "" ]]; then
			# arg="${ID} ${hstart} ${hwid}"
		# else
			# arg="${ID}"
		# fi
		# printf "${arg}\n"
		# Rscript --vanilla - <<-'EOF' ${arg}
			# #!/usr/bin/env Rscript
			# library(Gviz)
			# library(GenomicFeatures)
			# library(rtracklayer)
			
			# args<-commandArgs(trailingOnly = TRUE)
			
			# ID<-args[1]
			# filenames<-read.delim(paste0("figures_manuscript/data/filenames_",ID,".txt"), header=TRUE)
			# genes<-import(paste0("figures_manuscript/data/features_in_locus_",ID,".bed"))
			# dmrs<-import(paste0("figures_manuscript/data/DMRs_in_locus_",ID,".bed"))
			
			# options(ucscChromosomeNames=FALSE)

			# tracksize<-c(1,1)
			# tracklist<-list()
			# for ( i in c(1:nrow(filenames)) ) {
				# name<-filenames$Name[i]
				# group<-filenames$Group[i]
				# backcolor<-filenames$Backcolor[i]
				# trackcolor<-filenames$Trackcolor[i]
				# fillcolor1<-filenames$Fillcolor1[i]
				# fillcolor2<-filenames$Fillcolor2[i]
				# ymin<-filenames$Ymin[i]
				# ymax<-filenames$Ymax[i]
				# ymintick<-sign(ymin)*((floor(abs(ymin)*100)/100))
				# ymaxtick<-sign(ymax)*((floor(abs(ymax)*100)/100))
				# tracksize<-c(tracksize,1)
				# sample<-paste0(name,"_",group)
				# filename<-paste0("figures_manuscript/data/",sample,"_locus.bw")
				# print(paste0("Importing bw for ",sample))
				# bw<-import(filename)
				# print(paste0("Creating track for ",sample))
				# track<-DataTrack(bw, type="polygon", baseline=0, name=sample, background.title = backcolor, col=trackcolor, fill.mountain=c(fillcolor1,fillcolor2), col.baseline="grey50", ylim=c(ymin,ymax), yTicksAt=c(ymintick,ymaxtick), rotation.title=0, cex.title=0.5, lwd=0.01, hjust.title=1)
				# tracklist<-append(tracklist, track)
			# }

			# axistrack<-GenomeAxisTrack(scale=0.1, labelPos = "above")
			# genetrack<-AnnotationTrack(genes, name = "Annotations", fill = "lightblue", fontcolor.group = "darkblue", shape="fixedArrow", arrowHeadWidth=5, rotation.title=0, cex.title=0.5, lwd=0.1, col="darkblue",
				# just.group="above", group=genes$name, groupAnnotation = "group")
			# ddm1track<-AnnotationTrack(ddm1, name = "ddm1 DMRs", revertant="lightblue", stable="stableblue", shape="box", rotation.title=0, cex.title=0.5, lwd=0.1,
				# group=ddm1$name, groupAnnotation = "group")
			# dmrtrack<-AnnotationTrack(dmrs, name = "DMRs", fill="red", shape="box", rotation.title=0, cex.title=0.5, lwd=0.1)
			# tracks<-append(list(axistrack, genetrack), tracklist,list(dmrtrack))
			# tracksize<-c(tracksize,0.5)
			# pdf(paste0("figures_manuscript/draft_fin/",ID,".pdf"),paper="a4")
			# plotTracks(tracks, sizes=tracksize, title.width=2.5)
			# dev.off()

		# EOF
		# rm -f figures_manuscript/data/*${ID}*
	# done < figures_manuscript/data/${nameout}_loci.txt
# done


#######################################################################################################################################################################################

printf "Script finished!\n"

#######################################################################################################################################################################################
