#!/usr/bin/bash

## This pipeline was used to generate genome catalogs.
## Be sure you have installed drep, gtdbtk by conda(recommanded) and fasttree before use the pipeline.

### Set your software path, CPU use, input files and output path ###
#software and cpus
conda_drep="conda activate drep"
conda_gtdbtk="conda activate gtdbtk"
gtdbdb=/path/of/the/gtdb-db
## Please note that different versions of gtdb-db requires different gtdbtk versions to run, refer to https://ecogenomics.github.io/GTDBTk/installing/index.html for more information.
fasttree=/path/of/your/fasttree
drep_cpu=32
gtdbtk_cpu=8

#drep input file and parameters
self_genome_path_file=/path/of/your/genomes/path/file
self_genome_info_file=/path/of/your/genomes/quality/file
public_genome_path_file=/path/of/public/genomes/path/file
public_genome_info_file=/path/of/public/genomes/quality/file
suffix=fasta #suffix of genome, maybe fasta or fa
self_s_ani=0.99 #ANI threshold for drep of self data to form secondary clusters, 0.99 or 0.95
pub_s_ani=0.95 #ANI threshold for drep of merged public data to form secondary clusters, 0.99 or 0.95
com=50 #completeness
con=10 #contamination
pa=0.9 #ANI threshold to form primary (MASH) clusters
cov=0.3 #Minmum level of overlap between genomes when doing secondary comparisons

#output path
out=/path/of/output/
###########################################################

######################### Main ############################
#1. drep99 of self data
eval $conda_drep
if [ ! -d $out/1_drep99 ];then mkdir -p $out/1_drep99;fi
dRep dereplicate $out/1_drep99 --genomeInfo $self_genome_info_file -g $self_genome_path_file -p $drep_cpu -comp $com -con $con -pa $pa --S_ani $self_s_ani --cov_thresh $cov && \
ls $out/1_drep99/dereplicated_genomes > $out/1_drep99/drep99.list && \
grep -Ff $out/1_drep99/drep99.list $self_genome_path_file > $out/1_drep99/drep99.representative.path && \
grep -Ff $out/1_drep99/drep99.list $self_genome_info_file > $out/1_drep99/drep99.representative.genome.info && \
touch $out/1.drep99.done

#2. drep95 with public data
if [ ! -f $out/1.drep99.done ]
then
	echo "drep99 maybe not finished, please check and try again."
	exit
else
	eval $conda_drep
	if [ ! -d $out/2_drep95 ];then mkdir -p $out/2_drep95;fi
	cat $public_genome_path_file $out/1_drep99/drep99.representative.path > $out/drep95_input.path
	cat $public_genome_info_file $out/1_drep99/drep99.representative.genome.info > $out/drep95_input.genome.info
	dRep dereplicate $out/2_drep95 --genomeInfo $out/drep95_input.genome.info -g $out/drep95_input.path -p $drep_cpu -comp $com -con $con -pa $pa --S_ani $pub_s_ani --cov_thresh $cov && \
	ls $out/2_drep95/dereplicated_genomes > $out/2_drep95/drep95.list && \
	grep -Ff $out/2_drep95/drep95.list $out/drep95_input.path > $out/2_drep95/drep95.representative.path && \
	grep -Ff $out/2_drep95/drep95.list $out/drep95_input.genome.info > $out/2_drep95/drep95.representative.genome.info && \
	touch $out/2.drep95.done
fi
#3. gtdb
if [ ! -f $out/2.drep95.done ]
then
	echo "drep95 maybe not finished, please check and try again."
	exit
else
	eval $conda_gtdbtk
	export GTDBTK_DATA_PATH=$gtdbdb
	if [ ! -d $out/3_gtdb ];then mkdir -p $out/3_gtdb;fi
	gtdbtk classify_wf --genome_dir $out/2_drep95/dereplicated_genomes/ --out_dir $out/3_gtdb -x $suffix --cpus $gtdbtk_cpu --skip_ani_screen && touch $out/3.gtdb.done
fi

#4. fastree
if [ ! -f $out/3.gtdb.done ]
then
	echo "gtdb maybe not finished, please check and try again."
	exit
else
	if [ ! -d $out/4_fasttree ];then mkdir -p $out/4_fasttree;fi
	if [ -f $out/3_gtdb/align/gtdbtk.ar53.user_msa.fasta.gz ]
	then
		gunzip -c $out/3_gtdb/align/gtdbtk.ar53.user_msa.fasta.gz > $out/4_fasttree/ar53.fa
		chmod -R 755 $out/4_fasttree/ar53.fa
		$fasttree $out/4_fasttree/ar53.fa > $out/4_fasttree/ar53.tree && touch $out/4.fasttree.ar53.done
	fi
	if [ -f $out/3_gtdb/align/gtdbtk.bac120.user_msa.fasta.gz ]
	then
		gunzip -c $out/3_gtdb/align/gtdbtk.bac120.user_msa.fasta.gz > $out/4_fasttree/bac120.fa
		chmod -R 755 $out/4_fasttree/bac120.fa
		$fasttree $out/4_fasttree/bac120.fa > $out/4_fasttree/bac120.tree && touch $out/4.fasttree.bac120.done
	fi
fi
