#!/usr/bin/bash

## this pipeline is used to construct geneset and conduct the KEGG functional annotation and PETase identification.

####################### Parameters ########################
#$1 
#Sample assemble list for MetaGeneMark 
#$2
#A fasta file with IsPETase reference sequences
#$3
#Output directory, this is now set to the same for all steps. In the settings below, you can set them differently for each individual steps.
###########################################################

output_dir=$3
#################### Settings for metaGeneMark
meta_genemark_list=$1
meta_genemark_path=/home/your_meta_gene_mark_path/meta_genemark
meta_genemark_Outdir=$output_dir
pl_path=./bin
pl_gmm=p4_deal_gmm.pl
pl_transform=p4_transform.pl
pl_codon_table=p4_codon-table.11

#################### Settings for mmseqs
mmseqs_path=/home/your_mmseqs_path/bin
mmseqs_Outdir=$output_dir
mmseqs_output_filename=clust.all.gene_prot.fa
mmseqs_output_filename_cat=all.clust_split.gene_prot.fa
mmseqs_output_filename_final=GOPC.geneset.pep.fa

#################### Settings for kofamscan
kofamscan_env_path=/home/your_kofamscan_env_path
kofamscan_Outdir=$output_dir
kofamscan_input_filename=$mmseqs_output_filename_final\_rep_seq.fasta
kofamscan_input_filename_nospace=nospace_$kofamscan_input_filename
kofamscan_output_filename=kofamscan_out.geneset.pep.fa

#################### Settings for diamond
diamond_path=/home/your_diamond_path
diamond_ref_fastafile=$2
diamond_Outdir=$output_dir
diamond_ouput_filename=blast.your_diamond_output.m8
filtered_candidates_Outdir=$output_dir
filtered_output_filename=blast.filter.fa

#################### Settings for SignalP
signalP_path=/home/your_signalP_path/bin

#################### Test if all inputs needed extst.
if [ -z "$meta_genemark_list" ]
then
        echo "Please specify the sample assemble list for MetaGeneMark as the 1st input. :)"
        exist_list=false
else
        exist_list=true
fi
if [ -z "$diamond_ref_fastafile" ]
then
	echo "Please specify a fasta file with reference sequences as the 2nd input. :)"
        exist_ref=false
else
        exist_ref=true
fi
if [ -z "$output_dir" ]
then
        echo "Please specify an output directory as the 3rd input. :)"
        exist_out=false
else
	if [ ! -e "$output_dir" ]
	then
		mkdir $output_dir
		echo "$output_dir does not exist, created!"
	fi
        exist_out=true
fi
#################### Test if all tools needed exist.
if [ ! -e "$meta_genemark_path/gmhmmp" ]
then
         echo "Please specify the path of MetaGeneMark v3.38 to replace the value of meta_genemar_path. :)"
        exist_meta=false
else
        exist_meta=true
fi
if [ ! -e "$pl_path/$pl_gmm" ] || [ ! -e "$pl_path/$pl_transfrom" ] || [ ! -e "$pl_path/$pl_codon_table" ] 
then
        echo "Please download the .pl scripts and codon table (deal_gmm.pl, transform.pl and codon-table.11) from our github repository and put them in ./bin. Of course they can be put in other directory, then the value of pl_path need to be changed accordingly. :)"
        exist_pl=false
else
        exist_pl=true
fi
if [ ! -e "$mmseqs_path/mmseqs" ]
then
	echo "Please specify the path of MMseqs2 easy-linclust v12.113e3 to replace the value of mmseqs_path. :)"
        exist_mmseqs=false
else
        exist_mmseqs=true
fi
if [ ! -e "$kofamscan_env_path/bin/exec_annotation" ]
then
	echo "Please specify the path of kofamscan v1.3.0, KEGG database v87.0 to replace the value of kofamscan_env_path. if kofamscan has not been installed yes, please do so. :)"
        exist_kofamscan=false
else
	export PATH=$kofamscan_env_path/bin/:$PATH
        exist_kofamscan=true
fi
if [ ! -e "$diamond_path/diamond" ]
then
	echo "Please specify the path of DIAMOND v0.8.23.85 to replace the value of diamond_path. :)"
        exist_diamond=false
else
        exist_diamond=true
fi
if [ ! -e "$signalP_path/signalp" ]
then
	echo "Please specify the path of SignalP v5.0b to replace the value of signalP_path. :)"
        exist_singlep=false
else
        exist_singlep=true
fi


if $exist_list && $exist_ref && $exist_out && $exist_meta && $exist_pl && $exist_mmseqs && $exist_kofamscan && $exist_diamond && $exist_singlep
then
echo "Congratulations! Everything needed exists, start running!"
#MetaGeneMark v3.38, input sample assemble list
echo "["$(date +"%y-%m-%d %T")"] Running MetaGeneMark"
while IFS= read -r line
do
    {
	meta_genemark_sample_name=$(cut -d " " -f 1 <<< $line)
	meta_genemark_fasta_path=$(cut -d " " -f 2 <<< $line)
	meta_genemark_outputfile=$meta_genemark_sample_name.gmm
	pl_output_filename=$meta_genemark_sample_name
	$meta_genemark_path/gmhmmp -m $meta_genemark_path/MetaGeneMark_v1.mod -o $meta_genemark_Outdir/$meta_genemark_outputfile $meta_genemark_fasta_path
	$pl_path/$pl_gmm $meta_genemark_Outdir/$meta_genemark_outputfile $meta_genemark_fasta_path 100 GL $meta_genemark_Outdir/$pl_output_filename
	$pl_path/$pl_transform $pl_path/$pl_codon_table $meta_genemark_Outdir/$meta_genemark_sample_name.gene_nucl.fa $meta_genemark_Outdir/$meta_genemark_sample_name.gene_prot.fa
    } &
done < $meta_genemark_list
wait

#split 30 into parts
cat $meta_genemark_Outdir/*.gene_prot.fa >$meta_genemark_Outdir/all.gene_prot.fa
split -n 30 -d $meta_genemark_Outdir/all.gene_prot.fa $meta_genemark_Outdir/all.gene_prot.fa.

#MMseqs2 easy-linclust v12.113e3
echo "["$(date +"%y-%m-%d %T")"] Running MMseqs2"
for i in {00..30}
do
	$mmseqs_path/mmseqs easy-linclust $meta_genemark_Outdir/all.gene_prot.fa.$i $mmseqs_Outdir/$mmseqs_output_filename.$i $mmseqs_Outdir/temp.$i  --threads 30 --cov-mode 1 -c 0.99 --min-seq-id 0.95 &
done
wait

#merge MMseqs2  easy-linclust
cat $mmseqs_Outdir/$mmseqs_output_filename.*rep* >$mmseqs_Outdir/$mmseqs_output_filename_cat
$mmseqs_path/mmseqs easy-linclust $mmseqs_Outdir/$mmseqs_output_filename_cat $mmseqs_Outdir/$mmseqs_output_filename_final $mmseqs_Outdir/temp --threads 30 --cov-mode 1 -c 0.99 --min-seq-id 0.95

#remove spaces in file
sed -e 's/\s\+/_/g' $mmseqs_Outdir/$kofamscan_input_filename > $mmseqs_Outdir/$kofamscan_input_filename_nospace
#kofamscan v1.3.0, KEGG database v87.0
echo "["$(date +"%y-%m-%d %T")"] Running kofamscan"
$kofamscan_env_path/bin/exec_annotation -c $kofamscan_env_path/config.yaml -o $kofamscan_Outdir/$kofamscan_output_filename $mmseqs_Outdir/$kofamscan_input_filename_nospace --tmp-dir $kofamscan_Outdir/temp -T 0.7

#DIAMOND v0.8.23.85, E-value < 1e-5
echo "["$(date +"%y-%m-%d %T")"] Running DIAMOND"
$diamond_path/diamond --makedb -in $diamond_ref_fastafile -db $diamond_Outdir/$diamond_ref_fastafile.dmnd
$diamond_path/diamond blastp --evalue 1e-5 -q $mmseqs_Outdir/$kofamscan_input_filename_nospace -d $diamond_Outdir/$diamond_ref_fastafile.dmnd -o $diamond_Outdir/$diamond_ouput_filename --sensitive --max-target-seqs 10 --salltitles --outfmt 6 --threads 16

#get candidates
perl -e 'open(IN,"$diamond_Outdir/$diamond_ouput_filename");while(<IN>){chomp;@a=split/\t/;$h{$a[0]}{$a[1]}++;}close IN;$/=">";<>;while(<>){chomp;@a=split/\n/,$_,2;next if($a[0] !~ /Complete/);$a[0]=(split/\s+/,$a[0])[0];$l=length($a[1]);next if($l<250 || $l>350);@k=sort keys %{$h{$a[0]}};$k=join"\t",@k;print">$a[0]\t$k\n$a[1]";}' $mmseqs_Outdir/$kofamscan_input_filename_nospace >$filtered_candidates_Outdir/$filtered_output_filename

#SignalP v5.0b
echo "["$(date +"%y-%m-%d %T")"] Running SignalP"
#for archaea
$signalP_path/signalp -batch 20000 -fasta $filtered_candidates_Outdir/$filtered_output_filename -format short -org arch -prefix gram_arch -tmp $filtered_candidates_Outdir/tmp_arch
#for bacteria
$signalP_path/signalp -batch  20000 -fasta $filtered_candidates_Outdir/$filtered_output_filename -format short -org gram+ -prefix gram_pos -tmp $filtered_candidates_Outdir/tmp_pos
$signalP_path/signalp -batch  20000 -fasta $filtered_candidates_Outdir/$filtered_output_filename -format short -org gram- -prefix gram_neg -tmp $filtered_candidates_Outdir/tmp_neg
else
echo "Oops, something is missing!"
fi
