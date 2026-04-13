#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use FindBin '$Bin';
use Cwd;
use Cwd 'abs_path';
use File::Basename;
use Time::Local;

my ($MAGs,$Antismash,$Phylums,$outdir,$L2norm_threshold,$Endpoint,$Nboot,$SE,$Novetly_threshold,$project,$queue,$sproject,$squeue);
my ($run_s1,$run_s2,$run_s3,$run_s4,$run_s5,$run_s6);
my $para = $ARGV[0];
get_params("$para");

## 加载全局变量
my $pipeline_bin = "$Bin/bin";
require("$pipeline_bin/config.pl");
our($Rscript,$perl,$python,$condaac,$condaacc,$condadeac,$bigslice,$bigslice2,$bigfam);

## 读取系统时间
my ($min,$hour,$day,$month,$year)=(localtime)[1,2,3,4,5];
($min,$hour,$day,$month,$year) = (sprintf("%02d", $min),sprintf("%02d", $hour),sprintf("%02d", $day),sprintf("%02d", $month + 1),$year + 1900);

## 准备log文件
my $run_start_time="$year$month$day$hour$min";
my $run_cmd_log="run_cmd\_$run_start_time.log";
my $LOG=open_OUT1($run_cmd_log);

## 准备log和sh的输出文件夹
# $phylum_list ||= "ALL_PHYLUM";
$outdir ||= "./Output";
$outdir = abs_path($outdir);
Mkdir($outdir);
my $cwd = getcwd();
my $sd = Mkdir("$cwd/sh_$run_start_time");

################## Main #######################
#1. prepare for bigslice
prepare_bigslice() if $run_s1 eq "y";

#2. run bigslice
run_bigslice() if $run_s2 eq "y";;

#3. clustering by l2norm distance
run_ed_clustering() if $run_s3 eq "y";;
#run_cos_clustering();
#run_cos_clustering_0928();

#4. circlize
run_circlize() if $run_s4 eq "y";;

#5. rarefaction
run_rarefaction() if $run_s5 eq "y";;

#6. novetly GCFs
run_novetly() if $run_s6 eq "y";;


################## Sub Functions #######################
sub prepare_bigslice{
    my $job = "prepare_bigslice_input";
    my $sh_prepare_bigslice = "$sd/prepare_bigslice_input.sh";
    my $outdir_sub = Mkdir("$outdir/1_Bigslice_Input");
    # my $cmd = "$perl $pipeline_bin/get_gbk.pl -p $Phylums -a $Antismash -m $MAGs -o $outdir_sub";
    my $cmd = "$perl $pipeline_bin/prepare_bigslice_input.pl $MAGs $Antismash $outdir_sub";
    my @out_files = ("$outdir_sub/datasets.tsv","$outdir_sub/target_MAGs.list");
    print_cmd($cmd,$sh_prepare_bigslice,\@out_files);
    run_cmd($job,$sh_prepare_bigslice,1,16,4,$project,$queue);
}

sub run_bigslice{
    my $job = "bigslice";
    my $sh_bigslice = "$sd/bigslice.sh";
    my $outdir_sub = "$outdir/2_Bigslice_Output";
    my $outdirr_sub = "$outdir/2_Bigslice_Output/tsv_info";
    my $cmd1 = "$bigslice2 -i $outdir/1_Bigslice_Input $outdir_sub --threshold 9999 --complete -t 16";
    my @out_files = ("$outdir_sub/result/data.db");
    print_cmd($cmd1,$sh_bigslice,\@out_files);
	my $cmd2 = "$bigslice2 --export-tsv $outdirr_sub $outdir_sub";
	@out_files = ("$outdirr_sub/bgc_metadata.tsv", "$outdirr_sub/gcf_membership.tsv","$outdirr_sub/run_metadata.tsv");
	print_cmd($cmd2,$sh_bigslice,\@out_files);
    run_cmd($job,$sh_bigslice,2,16,300,$sproject,$squeue);
}

sub run_ed_clustering{
    my $job = "ed_clustering";
    my $sh_clustering = "$sd/ed_clustering.sh";
    my $outdir_sub = "$outdir/2_Bigslice_Output";
	my $cmd1 = $condaac;
    my $cmd2 = "$python $pipeline_bin/perform_l2norm_clustering.py $outdir_sub $L2norm_threshold $outdir_sub/GCF_models.tsv $outdir_sub/BGC_membership.tsv";
	my $cmd3 = "$Rscript $pipeline_bin/get_bgc_info.R $outdir_sub/result/data.db $outdir_sub/BGC_info.txt";
	my $cmd4 = "$perl $pipeline_bin/get_gcf_bgc_info.pl -m $outdir/1_Bigslice_Input/target_MAGs.list -b $outdir_sub/BGC_info.txt -g $outdir_sub/BGC_membership.tsv -o $outdir_sub/GCF_info.txt";
	my $cmd5 = $condadeac;
    my @out_files = ("$outdir_sub/GCF_info.txt");
    print_cmd($cmd1,$sh_clustering,\@out_files);
    print_cmd($cmd2,$sh_clustering,\@out_files);
    print_cmd($cmd3,$sh_clustering,\@out_files);
    print_cmd($cmd4,$sh_clustering,\@out_files);
    print_cmd($cmd5,$sh_clustering,\@out_files);
    run_cmd($job,$sh_clustering,5,16,64,$sproject,$squeue);
}

sub run_cos_clustering{
    my $job = "run_cos_clustering";
    my $sh_clustering = "$sd/cos_clustering.sh";
    my $outdir_sub = "$outdir/2_Bigslice_Output";
	Mkdir("$outdir_sub/temp");
	my $cmd1 = $condaac;
    my $cmd2 = "$python $pipeline_bin/export_matrixs.py $outdir_sub/result $outdir_sub/result $L2norm_threshold $outdir_sub/temp/BGC_features.tsv $outdir_sub/temp/GCF_Models.tsv";
	my $cmd3 = "$perl $pipeline_bin/prepare_for_cal_cos.pl $outdir_sub/temp/BGC_features.tsv $outdir_sub/temp/GCF_Models.tsv $outdir_sub/temp/BGC_features_input.tsv $outdir_sub/temp/GCF_models_input.tsv";
	my $cmd4 = "$python $pipeline_bin/cal_cos_distance.py $outdir_sub/temp/BGC_features_input.tsv $outdir_sub/temp/GCF_models_input.tsv $outdir_sub/Cos_Similarity.txt";
	my $cmd5 = "$perl $pipeline_bin/get_min_cos_dis.pl $outdir_sub/temp/GCF_header.list $outdir_sub/Cos_Similarity.txt $outdir_sub/Cos_Distance_Info.txt $outdir_sub/Cos_Distance_to_GCF_Models.txt";
	my $cmd6 = "$perl $pipeline_bin/get_gcf_bgc_info.pl -m $outdir/1_Bigslice_Input/target_MAGs.list -b $outdir_sub/BGC_info.txt -g $outdir_sub/Cos_Distance_to_GCF_Models.txt -o $outdir_sub/cos_GCF_info.txt";
	my $cmd7 = $condadeac;
    my @out_files = ("$outdir_sub/cos_GCF_info.txt");
    print_cmd($cmd1,$sh_clustering,\@out_files);
    print_cmd($cmd2,$sh_clustering,\@out_files);
    print_cmd($cmd3,$sh_clustering,\@out_files);
    print_cmd($cmd4,$sh_clustering,\@out_files);
    print_cmd($cmd5,$sh_clustering,\@out_files);
    print_cmd($cmd6,$sh_clustering,\@out_files);
    print_cmd($cmd7,$sh_clustering,\@out_files);
    run_cmd($job,$sh_clustering,7,16,64);
}

sub run_cos_clustering_0928{
    my $job = "run_cos_clustering";
    my $sh_clustering = "$sd/cos_clustering.sh";
    my $outdir_sub = "$outdir/2_Bigslice_Output";
	Mkdir("$outdir_sub/temp");
	my $cmd1 = $condaac;
    my $cmd2 = "$python $pipeline_bin/extract_bgc_features_matrix.py $outdir/2_Bigslice_Output $outdir_sub/temp/BGC_features.tsv";
	my $cmd3 = "$perl $pipeline_bin/prepare_bgc_for_cal_cos.pl $outdir_sub/temp/BGC_features.tsv $outdir_sub/temp/BGC_features_input.tsv";
	my $cmd4 = "$python $pipeline_bin/cos_cluster.py $outdir_sub/temp/BGC_features_input.tsv $outdir_sub/temp/Trees_Info.txt";
	my $cmd5 = "$perl $pipeline_bin/get_cos_cluster.pl $outdir_sub/temp/Trees_Info.txt $outdir_sub/BGC_membership.txt";
	my $cmd6 = "$Rscript $pipeline_bin/get_bgc_info.R $outdir_sub/result/data.db $outdir_sub/BGC_info.txt";
	my $cmd7 = "$perl $pipeline_bin/get_gcf_bgc_info.pl -m $outdir/1_Bigslice_Input/target_MAGs.list -b $outdir_sub/BGC_info.txt -g $outdir_sub/BGC_membership.txt -o $outdir_sub/GCF_info.txt";
	my $cmd8 = $condadeac;
    my @out_files = ("$outdir_sub/GCF_info.txt");
    print_cmd($cmd1,$sh_clustering,\@out_files);
    print_cmd($cmd2,$sh_clustering,\@out_files);
    print_cmd($cmd3,$sh_clustering,\@out_files);
    print_cmd($cmd4,$sh_clustering,\@out_files);
    print_cmd($cmd5,$sh_clustering,\@out_files);
    print_cmd($cmd6,$sh_clustering,\@out_files);
    print_cmd($cmd7,$sh_clustering,\@out_files);
    print_cmd($cmd8,$sh_clustering,\@out_files);
    run_cmd($job,$sh_clustering,8,16,64);
}

sub run_circlize{
    my $job = "circlize";
    my $sh_circlize = "$sd/circlize.sh";
    my $outdir_sub = Mkdir("$outdir/3_Circlize");
	my $cmd1 = $condaac;
	my $cmd2 = "$perl $pipeline_bin/get_share_gcfs.pl $outdir/2_Bigslice_Output/GCF_info.txt $outdir_sub";
	my $cmd3 = "$perl $pipeline_bin/get_circlize_input.pl $outdir_sub";
	my $cmd4 = "$perl $pipeline_bin/get_circlize_info.pl $outdir_sub";
	my $cmd5 = "$Rscript $pipeline_bin/plot_circlize.R $outdir_sub/Circlize_Input/Chr_len.txt $outdir_sub/Circlize_Input/Links_input.txt $outdir_sub/Circlize.pdf";
	my $cmd6 = $condadeac;
	my @out_files = ("$outdir_sub/Circlize.pdf","$outdir_sub/GCF_counts_per_phylim_total.txt","$outdir_sub/GCF_counts_per_phylim_uniq.txt");
	print_cmd($cmd1,$sh_circlize,\@out_files);
	print_cmd($cmd2,$sh_circlize,\@out_files);
	print_cmd($cmd3,$sh_circlize,\@out_files);
	print_cmd($cmd4,$sh_circlize,\@out_files);
	print_cmd($cmd5,$sh_circlize,\@out_files);
	print_cmd($cmd6,$sh_circlize,\@out_files);
	run_cmd($job,$sh_circlize,6,16,16,$project,$queue);
}

sub run_rarefaction{
    my $job = "rarefaction_analysis";
    my $sh_rarefaction = "$sd/rarefaction.sh";
    my $outdir_sub = Mkdir("$outdir/4_Rarefaction");
	my $cmd1 = $condaac;
	my $cmd2 = "$perl $pipeline_bin/get_inext_input.pl $outdir/1_Bigslice_Input $outdir/2_Bigslice_Output $outdir_sub";
	my $cmd3 = "$Rscript $pipeline_bin/iNEXT.R $outdir_sub $Endpoint $Nboot $SE";
	my $cmd4 = "$perl $pipeline_bin/get_inext_output_info.pl $outdir_sub";
	my $cmd5 = $condadeac;
	my @out_files = ("$outdir_sub/Rarefaction_result.pdf","$outdir_sub/iNEXT_Info.txt");
	print_cmd($cmd1,$sh_rarefaction,\@out_files);
	print_cmd($cmd2,$sh_rarefaction,\@out_files);
	print_cmd($cmd3,$sh_rarefaction,\@out_files);
	print_cmd($cmd4,$sh_rarefaction,\@out_files);
	print_cmd($cmd5,$sh_rarefaction,\@out_files);
	run_cmd($job,$sh_rarefaction,5,16,32,$project,$queue);
}

sub run_novetly{
    my $job = "novetly";
    my $sh_novetly = "$sd/novetly.sh";
    my $outdir_sub = Mkdir("$outdir/5_Novetly");
	Mkdir("$outdir_sub/temp");
	my $cmd1 = "$python $pipeline_bin/extract_bgc_features_matrix.py $bigfam $bigfam/Big_fam_BGC_feature.txt";
	my @out_files = ("$bigfam/Big_fam_BGC_feature.txt");
	print_cmd($cmd1,$sh_novetly,\@out_files);
	my $cmd2 = "$perl $pipeline_bin/prepare_bgc_for_cal_cos.pl $bigfam/Big_fam_BGC_feature.txt $bigfam/Big_fam_BGC_feature_input.txt";
	@out_files = ("$bigfam/Big_fam_BGC_feature_input.txt");
	print_cmd($cmd2,$sh_novetly,\@out_files);
	my $cmd3 = "$python $pipeline_bin/extract_bgc_features_matrix.py $outdir/2_Bigslice_Output $outdir_sub/temp/All_BGC_feature.txt";
	my $cmd4 = "$perl $pipeline_bin/Novelty/0_sort_hmm_order.pl $bigfam/Big_fam_BGC_feature.txt $outdir_sub/temp/All_BGC_feature.txt $outdir_sub/temp/All_BGC_feature_sorted.txt";
	my $cmd5 = "$perl $pipeline_bin/prepare_bgc_for_cal_cos.pl $outdir_sub/temp/All_BGC_feature_sorted.txt $outdir_sub/0_All_BGC_feature_sorted_input.txt";
	my $cmd6 = "$python $pipeline_bin/Novelty/1_cal_dis.py $outdir_sub/0_All_BGC_feature_sorted_input.txt $bigfam/Big_fam_BGC_feature_input.txt $outdir_sub/1_All_BGC_Dist.txt";
	my $cmd7 = "$perl $pipeline_bin/Novelty/2_get_min_cos_dis.pl $outdir/2_Bigslice_Output/GCF_info.txt $outdir_sub/1_All_BGC_Dist.txt $outdir_sub/2_All_BGC_Min_Dist.txt $outdir_sub/3_All_GCF_Novelty.txt $Novetly_threshold";
	@out_files = ("$outdir_sub/1_All_BGC_Dist.txt","$outdir_sub/2_All_BGC_Min_Dist.txt","$outdir_sub/3_All_GCF_Novelty.txt");
	print_cmd($cmd3,$sh_novetly,\@out_files);
	print_cmd($cmd4,$sh_novetly,\@out_files);
	print_cmd($cmd5,$sh_novetly,\@out_files);
	print_cmd($cmd6,$sh_novetly,\@out_files);
	print_cmd($cmd7,$sh_novetly,\@out_files);
	run_cmd($job,$sh_novetly,7,16,128,$project,$queue);
}


sub Mkdir {
    my $path = shift;
    `mkdir -p $path` unless (-d $path);
    return $path;
}

sub print_cmd1 {
	my ($cmd,$out_sh,$out_file) = @_;
	
	my @out_file=@$out_file;
	
	open(OUT,">>$out_sh");
	if($cmd=~/^cd / || $cmd=~/^ln /){
		print OUT "$cmd\n";
	}
	else{
		my $tmp = $cmd;
        my $flag=0;
        foreach my $out_file(@out_file){
            unless(-e $out_file){
                $flag=1;
            }
        }
        if($flag==0){
            $tmp="# $cmd";
        }
		print OUT "$tmp\n";
	}
	close OUT;
	`sleep 0.1s`;
}

sub print_cmd {
        my ($cmd,$out_sh,$out_file) = @_;
        my @out_file=@$out_file;
        open(OUT,">>$out_sh");
        print OUT "$cmd\n";
        close OUT;
	`sleep 0.1s`;
}


sub run_cmd {
	my ($job,$sh,$ls,$p,$vf,$project,$queue)=@_;
	if(-e $sh){
		my ($start,$start_h) = get_start_time();
		my $ret = 0;
		#my $cmd = "sh $sh";
		my $cmd = "perl /dellfsqd2/ST_OCEAN/USER/zhouchanghao/script/qsub/qsub_pipline.pl -b $ls -l vf=${vf}G,p=$p -P $project -q $queue -N $job -m 100 $sh";
		print $LOG "#$job:\n";
		print $LOG "$sh\n";
		print STDERR "Start $sh at $start_h ...\n";
		
		$ret=system($cmd);

		my ($end,$end_h) = get_end_time();
		my $dur = get_duration($start,$end);
		my $status = "Done";
		if ($ret != 0) {
			$status = "Stop";
		}
		
		print $LOG "#Run status: $status\n\n";
		print STDERR "$status at $end_h.\n";
		print STDERR "Duration: $dur\n\n";
		if ($ret != 0) {
			# exit(-1);
			exit;
		}
	}
}

sub get_duration {
	my ($start,$end)=@_;
	my $sec = $end - $start;
	my $days = int($sec/(24*60*60));
	my $hours = ($sec/(60*60))%24;
	my $mins = ($sec/60)%60;
	my $secs = $sec%60;
	
	return "$days days $hours hours $mins minutes $secs seconds";
}

sub get_start_time {
	my $start = time();
	my $start_h = localtime();
	return ($start,$start_h);
}

sub get_end_time {
	my $end = time();
	my $end_h = localtime();
	return ($end,$end_h);
}

sub open_OUT1 {
	my ($file) = @_;
	my $OUT = ();
	if($file=~/\.gz/){
		open($OUT,"| gzip >$file") or die "could not write $file $!";
	}else{
		open($OUT,">$file") or die "could not write $file $!";
	}	
	return $OUT;
}

sub get_params{
	my $paramfile_tmp=shift;
	my %paramfile_hash=();
	my @faults=();
	open (INF, "$paramfile_tmp" ) || die "cannot open parameter file $paramfile_tmp: $!\n";
	while(<INF>){
		next if ($_ =~ /^\#/);
		next unless ($_ =~ /=/);
		chomp $_;
		my($key, $value) = split('=', $_,2);
		$key=~s/\s//g;
		$value=~s/\s+$//;
		$value=~s/^\s+//;
		$paramfile_hash{$key} = $value;
	}
	close INF;
	$MAGs = $paramfile_hash{MAGs};
	$Antismash = $paramfile_hash{Antismash};
	$Phylums = $paramfile_hash{Phylums};
	$outdir = $paramfile_hash{Output};
	$L2norm_threshold = $paramfile_hash{L2norm_threshold};
	$Endpoint = $paramfile_hash{Endpoint};
	$Nboot = $paramfile_hash{Nboot};
	$SE = $paramfile_hash{SE};
	$Novetly_threshold = $paramfile_hash{Novetly_threshold};
        $run_s1 = $paramfile_hash{RUN_PREPARE};
        $run_s2 = $paramfile_hash{RUN_BIGSLICE};
        $run_s3 = $paramfile_hash{RUN_ED_CLUSTERING};
        $run_s4 = $paramfile_hash{RUN_CIRCLIZE};
        $run_s5 = $paramfile_hash{RUN_RAREFACTION};
	$run_s6 = $paramfile_hash{RUN_NOVELTY};
	$project = $paramfile_hash{PROJECT};
	$queue = $paramfile_hash{QUEUE};
	$sproject = $paramfile_hash{SUPER_PROJECT};
	$squeue = $paramfile_hash{SUPER_QUEUE};
}
