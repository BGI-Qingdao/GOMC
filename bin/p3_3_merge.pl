use strict;
#perl $0 <att.result> <lstm.result> <bert.result> <orf file> <output filtered orf list>

my ($att,$lstm,$bert,$orf,$gbk_file,$out) = @ARGV;

my %info;
open ATT, $att or die $!;
while(<ATT>){
	chomp;
	$info{att}{$.} = $_;
}
close ATT;
open LSTM, $lstm or die $!;
while(<LSTM>){
	chomp;
	$info{lstm}{$.} = $_;
}
close LSTM;
open BERT, $bert or die $!;
while(<BERT>){
	chomp;
	$info{bert}{$.} = $_;
}
close BERT;

open GBK, $gbk_file or die $!;
my %ii;
while(<GBK>){
	chomp;
	next if $. == 1;
	my @tmp = split("\t");
	next if $tmp[-1] eq "-";
	my $gbk = ($tmp[1]=~/(.*)\.gbk/)[0];
	my $value = join("\t", $tmp[0],$gbk,@tmp[2..6]);
	$ii{$tmp[0]}{$gbk} = $value;
}
close GBK;

open ORF, $orf or die $!;
open OUT, ">$out";
print OUT "AMP	Sample	GBK	GBK_size(Kb)	GBK_type	GBK_subtype	Leader_seq	Core_seq	Attention	LSTM	Bert\n";
while(<ORF>){
	chomp;
	if (/^>/){
		my $num = int($./2)+1;
		my $id = ($_=~/>(.*)/)[0];
		my @tmp = split("\\|", $id);
		my @values = split("\t", $ii{$tmp[0]}{$tmp[1]});
		my @leaders = split(",", $values[-2]);
		my @cores = split(",", $values[-1]);
		my $site = ($tmp[-1]=~/core_peptide_(.*)/)[0]-1;
		$values[-2] = $leaders[$site];
		$values[-1] = $cores[$site];
		my $line = join("\t", $id,@values,$info{att}{$num},$info{lstm}{$num},$info{bert}{$num});
		print OUT "$line\n";

		# if ($info{att}{$num}>0.5 && $info{lstm}{$num}>0.5 && $info{bert}{$num}>0.5){
		# 	my $id = ($_=~/>(.*)/)[0];
		# 	print OUT "$id\t$info{att}{$num}\t$info{lstm}{$num}\t$info{bert}{$num}\n";
		# }
	}
}
close ORF;
close OUT;