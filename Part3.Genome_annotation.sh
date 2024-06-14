#!/usr/bin/bash

## this pipeline is used to predict CDS and annotate functions of them.

####################### Parameters ########################
# input, output and bin path
input= #input fasta file
out= #output dir
prefix= #prefix of genome
bin=bin/ #path of scripts used in this progess

#prokka
conda_prokka="conda activate prokka"
centre=X
prokka_cpu=8

#kofamscan
conda_kofamscan="conda activate kofamscan"
kofamscan_config= #config.yaml path
threshold_scale=0.7

#rgi
rgi= #rgi path
card_json= #CARD database json path
alignment_tool=BLAST #DIAMOND or BLAST
rgi_cpu=8

#CRISPRCasTyper
conda_cctyper="conda activate cctyper"
cctyper= #cctyper path
cctyper_db= #cctyper database path
cctyper_cpu=8

#acafinder
conda_acafinder="conda activate acafinder"
acafinder= #AcaFind_runner.py path
vibrant_db= #VIBRANT database path
acr_protein_database= #The Acr proteins that will be used search for Acas, default are the published Acrs + AcrHub predicted Acrs + 2500 high confident Acr prediction of AcrCatalog
published_acahmm= #HMM for all 12 publsihed Aca proteins, recommended to use the default hmm provided from acafinder
phamdir= #Directory of all pfam hmm files with .dat files and other binaries
all_protein_length_in_acraca_operon=1000 #max proten length in Acr-Aca operon when length of Acr homolog < 200aa
intergenic_dist_in_acraca_operon=500 #Maximum Intergenic distance in Acr-Aca operon
acr_aca_inbetweengenes=15 #Maximum number of genes allowed between Aca and Acr proteins + 1 (e.g if the input is 4, then maximum 3 genes are allowed between the potental Aca genes to its closest Acr homolog)
acafinder_cpu=8

#MobileElementFinder
conda_mefinder="conda activate mefinder"
mefiner_cpu=8

#OGT_prediction
external_tools= #Create a file called external_tools.txt listing bedtools, tRNAscan-SE, barrnap, and prodigal tab separated from the absolute path for each executable.
conda_ogt="conda activate OGT"
ogt= #prediction_pipeline.py
taxonomy= #genome's taxonomy info(GTDBtk format)
model= # OGT models path

#antismash
conda_antismash="conda activate antismash"
taxon=bacteria #bacteria or fungi
genefinding_tool=prodigal #prodigal(bacteria) or glimmerhmm(fungi)
antismash_db= #antismash database path
antismash_cpu=8

#Attention,LSTM,BERT
conda_amp="conda activate AMP"
c_amp_models=c_AMP_Models/ #c_AMPs-prediction used models' path
apd3=c_AMP_Models/APD3.fa #known AMP seqs downloaded from APD3 database, or you can download the lastest version in https://aps.unmc.edu/downloads

###########################################################

######################### Main ############################
real_out=`realpath $out`

#1. prokka
eval $conda_prokka
if [ ! -d $real_out/1_prokka/$prefix ];then mkdir -p $real_out/1_prokka/$prefix;fi
prokka $input --addgenes --centre $centre  --compliant --prefix $prefix --force --cpu $prokka_cpu --outdir $real_out/1_prokka/$prefix && touch $real_out/1_prokka/prokka.$prefix.done

#2. kofamscan
if [ ! -f $real_out/1_prokka/prokka.$prefix.done ]
then
	echo "$prefix prokka maybe not finished,please check and try again."
	exit
elif [ ! -f $real_out/1_prokka/$prefix/$prefix.faa ]
then
	echo "$real_out/1_prokka/$prefix/$prefix.faa could not found, please check your protein file and set file path according to previous."
	exit
else
	eval $conda_kofamscan
	if [ ! -d $real_out/2_kofamscan/$prefix/tmp ];then mkdir -p $real_out/2_kofamscan/$prefix/tmp;fi
	exec_annotation $real_out/1_prokka/$prefix/$prefix.faa -c $kofamscan_config -o $real_out/2_kofamscan/$prefix/$prefix.kofam.result -T $threshold_scale --tmp-dir $real_out/2_kofamscan/$prefix/tmp && \
	rm -rf $real_out/2_kofamscan/$prefix/tmp && touch $real_out/2_kofamscan/kofamscan.$prefix.done
fi

#3. rgi
if [ ! -f $real_out/1_prokka/prokka.$prefix.done ]
then
	echo "$prefix prokka maybe not finished,please check and try again."
	exit
elif [ ! -f $real_out/1_prokka/$prefix/$prefix.faa ]
then
	echo "$real_out/1_prokka/$prefix/$prefix.faa could not found, please check your protein file and set file path according to previous."
	exit
else
	if [ ! -d $real_out/3_rgi/$prefix ];then mkdir -p $real_out/3_rgi/$prefix;fi
	cd $real_out/3_rgi/$prefix
	$rgi load -i $card_json --local && \
	$rgi main -i $real_out/1_prokka/$prefix/$prefix.faa -o $prefix.anno -a $alignment_tool --include_loose --clean --local -t protein -n $rgi_cpu && \
	cd - && touch $real_out/3_rgi/rgi.$prefix.done
fi

#4. CRISPRCasTyper
if [ ! -f $real_out/1_prokka/prokka.$prefix.done ]
then
	echo "$prefix prokka maybe not finished,please check and try again."
	exit
elif [ ! -f $real_out/1_prokka/$prefix/$prefix.fna ]
then
	echo "$real_out/1_prokka/$prefix/$prefix.fna could not found, please check your fasta file and set file path according to previous."
	exit
else
	eval $conda_cctyper
	if [ ! -d $real_out/4_cctyper/ ];then mkdir -p $real_out/4_cctyper/;fi
	if [ -d $real_out/4_cctyper/$prefix ];then rm -rf $real_out/4_cctyper/$prefix;fi
	$cctyper --db $cctyper_db $real_out/1_prokka/$prefix/$prefix.fna $real_out/4_cctyper/$prefix -t $cctyper_cpu && touch $real_out/4_cctyper/cctyper.$prefix.done
fi

#5. acafinder
if [ ! -f $real_out/1_prokka/prokka.$prefix.done ]
then
	echo "$prefix prokka maybe not finished,please check and try again."
	exit
elif  [ -f $real_out/1_prokka/$prefix/$prefix.fna ] && [ -f $real_out/1_prokka/$prefix/$prefix.gff ] && [ -f $real_out/1_prokka/$prefix/$prefix.faa ]
then
	eval $conda_acafinder
	if [ ! -d $real_out/5_acafinder/$prefix ];then mkdir -p $real_out/5_acafinder/$prefix;fi
	export VIBRANT_DATA_PATH=$vibrant_db
	export CCTYPER_DB=$cctyper_db
	python $acafinder --FNA_file $real_out/1_prokka/$prefix/$prefix.fna --GFF_file $real_out/1_prokka/$prefix/$prefix.gff --FAA_file $real_out/1_prokka/$prefix/$prefix.faa -o $real_out/5_acafinder/$prefix \
	-r $acr_protein_database -z $phamdir -y $published_acahmm -l $all_protein_length_in_acraca_operon -i $intergenic_dist_in_acraca_operon -b $intergenic_dist_in_acraca_operon -d $acafinder_cpu && \
	touch $real_out/5_acafinder/acafinder.$prefix.done
else
	echo "$prefix fna, gff or protein was missing at least one, please check the files path."
	exit
fi

#6. MobileElementFinder
if [ ! -f $real_out/1_prokka/prokka.$prefix.done ]
then
	echo "$prefix prokka maybe not finished,please check and try again."
	exit
elif [ ! -f $real_out/1_prokka/$prefix/$prefix.fna ]
then
	echo "$real_out/1_prokka/$prefix/$prefix.fna could not found, please check your fasta file and set file path according to previous."
	exit
else
	eval $conda_mefinder
	if [ ! -d $real_out/6_MobileElementFinder/$prefix ];then mkdir -p $real_out/6_MobileElementFinder/$prefix;fi
	if [ ! -d $real_out/6_MobileElementFinder/${prefix}_tmp ];then mkdir -p $real_out/6_MobileElementFinder/${prefix}_tmp;fi
	mefinder find -c $real_out/1_prokka/$prefix/$prefix.fna $real_out/6_MobileElementFinder/$prefix/$prefix --temp-dir $real_out/6_MobileElementFinder/${prefix}_tmp -g -t $mefiner_cpu && \
	rm -rf $real_out/6_MobileElementFinder/${prefix}_tmp && touch $real_out/6_MobileElementFinder/MobileElementFinder.$prefix.done
fi

#7. OGT
if [ ! -f $real_out/1_prokka/prokka.$prefix.done ]
then
	echo "$prefix prokka maybe not finished,please check and try again."
	exit
elif [ ! -f $real_out/1_prokka/$prefix/$prefix.fna ]
then
	echo "$real_out/1_prokka/$prefix/$prefix.fna could not found, please check your fasta file and set file path according to previous."
	exit
else
	eval $conda_ogt
	if [ ! -d $real_out/7_OGT/$prefix ];then mkdir -p $real_out/7_OGT/$prefix;fi
	cd $real_out/7_OGT/$prefix && \
	if [ ! -f external_tools.txt ];then ln -s $external_tools;fi && \
	if [ ! -d genomes/$prefix ];then mkdir -p genomes/$prefix;fi && \
	gzip -c $real_out/1_prokka/$prefix/$prefix.fna > genomes/$prefix/$prefix.fa.gz && \
	echo -e "$prefix.fa.gz\t$prefix" > genomes_info.txt && \
	header="species	superkingdom	phylum	class	order	family" && \
	tax=`echo $taxonomy | sed 's/__;/__None;/g' | sed 's/d__\(.*\);p__\(.*\);c__\(.*\);o__\(.*\);f__\(.*\);g__.*/\1\t\2\t\3\t\4\t\5/'` && \
	echo -e "$header\n$prefix\t$tax" > tax_info.txt && \
	python $ogt $model genomes_info.txt tax_info.txt && \
	cd - && touch $real_out/7_OGT/OGT.$prefix.done
fi

#8. antismash
if [ ! -f $real_out/1_prokka/prokka.$prefix.done ]
then
	echo "$prefix prokka maybe not finished,please check and try again."
	exit
elif [ ! -f $real_out/1_prokka/$prefix/$prefix.fna ] || [ ! -f $real_out/1_prokka/$prefix/$prefix.gff ]
then
	echo "$prefix fna or gff was missing at least one, please check the files path."
	exit
else
	eval $conda_antismash
	if [ ! -d $real_out/8_antismash/$prefix ];then mkdir -p $real_out/8_antismash/$prefix;fi
	antismash -c $antismash_cpu --taxon $taxon --genefinding-tool $genefinding_tool --cb-general --cb-subclusters --cb-knownclusters --asf --pfam2go \
	--output-dir $real_out/8_antismash/$prefix --genefinding-gff3 $real_out/1_prokka/$prefix/$prefix.gff $real_out/1_prokka/$prefix/$prefix.fna \
	--html-title $prefix --databases $antismash_db && \
	ls $real_out/8_antismash/$prefix/*region*gbk | while read id;do gbk=`echo $id | awk -F  "/" '{print $(NF)}'`;real=`realpath $id`;echo -e "$prefix\t$gbk\t$real";done > $real_out/8_antismash/$prefix.gbk.path && \
	touch $real_out/8_antismash/antismash.$prefix.done
fi

#9. Attention, LSTM and BERT
if [ ! -f $real_out/8_antismash/antismash.$prefix.done ]
then
	echo "$prefix antismash maybe not finished,please check and try again."
	exit
else
	eval $conda_amp
	export PYTORCH_PRETRAINED_BERT_CACHE="$c_amp_models/used"
	if [ ! -d $real_out/9_AMP/$prefix/predict ];then mkdir -p $real_out/9_AMP/$prefix/predict;fi
	perl $bin/p3_1_get_gbk_info.pl $real_out/8_antismash/$prefix.gbk.path $real_out/9_AMP/$prefix/1_GBK_info.xls && \
	perl $bin/p3_2_get_fa.pl $real_out/9_AMP/$prefix/1_GBK_info.xls $real_out/9_AMP/$prefix/2_core_peptides.fa && \
	perl $bin/p3_3_format.pl $real_out/9_AMP/$prefix/2_core_peptides.fa none > $real_out/9_AMP/$prefix/predict/Core_peptides.format && \
	python $bin/p3_3_prediction_attention.py $real_out/9_AMP/$prefix/predict/Core_peptides.format $real_out/9_AMP/$prefix/predict/Core_peptides.att && \
	python $bin/p3_3_prediction_lstm.py $real_out/9_AMP/$prefix/predict/Core_peptides.format $real_out/9_AMP/$prefix/predict/Core_peptides.lstm && \
	python $bin/p3_3_prediction_bert.py $real_out/9_AMP/$prefix/2_core_peptides.fa $real_out/9_AMP/$prefix/predict/Core_peptides.bert && \
	perl $bin/p3_3_merge.pl $real_out/9_AMP/$prefix/predict/Core_peptides.att $real_out/9_AMP/$prefix/predict/Core_peptides.lstm $real_out/9_AMP/$prefix/predict/Core_peptides.bert \
	$real_out/9_AMP/$prefix/2_core_peptides.fa $real_out/9_AMP/$prefix/1_GBK_info.xls $real_out/9_AMP/$prefix/predict/Core_peptides_predict_result.xls && \
	perl $bin/p3_4_filter.pl $real_out/9_AMP/$prefix/predict/Core_peptides_predict_result.xls $apd3 $real_out/9_AMP/$prefix/3_AMP_info.xls && \
	sed '1d' $real_out/9_AMP/$prefix/3_AMP_info.xls | cut -f 1 > $real_out/9_AMP/$prefix/AMP.list && \
	perl $bin/p3_5_extract_seq.pl -f $real_out/9_AMP/$prefix/2_core_peptides.fa -l $real_out/9_AMP/$prefix/AMP.list -o $real_out/9_AMP/$prefix/3_AMP.fa && rm $real_out/9_AMP/$prefix/AMP.list && \
	touch $real_out/9_AMP/AMP.$prefix.done
fi
