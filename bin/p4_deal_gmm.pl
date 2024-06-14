#!/usr/bin/perl -w
use strict;

die "usage: $0 [gmm] [fasta] [length] [mark] [prefix]\n"
unless (@ARGV == 5);

my ($gmm, $fa, $min_len, $mark, $pre) = @ARGV;


my(@info, $seqid, $i, @temp, $flag, $tmp, %genes, $count, $seq, $len);

open I, $gmm or die "failed to open file $gmm $!\n";
$/ = "FASTA definition line: ";
<I>;
while(<I>) {
	chomp;
	@info = split /\n/;
	$info[0] =~ /^(\S+)\s*.*/;
	$seqid = $1;

	for($i = 4; $i <= $#info; $i++){
		$flag=0;
		@temp = split /\s+/, $info[$i];
		next if @temp < 6;
		unshift @temp, "" if @temp == 6;
		if($temp[3] =~ /\<(\d+)$/){
			$flag = 1;
			$temp[3] = $1;
			if($temp[4] =~ /\>(\d+)$/){
				$flag = 2;
				$temp[4] = $1;
			}
		}elsif($temp[4] =~ /\>(\d+)$/){
			$flag = 3;
			$temp[4] = $1;
		}

		$tmp = $temp[3]."\t".$temp[4]."\t".$temp[2]."\t0\t";
		if(0 == $flag){
			$tmp .= "Complete";
		}elsif(1 == $flag){
			if($temp[2] eq '+'){
				$tmp .= "Lack 5\'-end";
			}else{
				$tmp .= "Lack 3\'-end";
			}
		}elsif(2 == $flag){
			$tmp .= "Lack both ends";
		}elsif(3 == $flag){
			if($temp[2] eq '+') {
				$tmp .= "Lack 3\'-end";
			}else{
				$tmp .= "Lack 5\'-end";
			}
		}
		push(@{$genes{$seqid}}, $tmp);
	}
}
$/ = "\n";
close I;

open I, $fa or die "failed to open file $fa $!\n";
open O1,">$pre.gff" or die "failed to open file $pre.gff $!\n";
print O1 "##gff-version 3\n";
open O2,">$pre.gene_nucl.fa" or die "failed to open file $pre.gene_nucl.fa $!\n";

$/ = ">";
<I>;
$count=1;

while(<I>){
	chomp;
	@info = split /\n/;
	$seqid = shift(@info);
	@temp = split /\s+/, $seqid;
	$seqid = $temp[0];
	$seq = join "", @info;
	$len = length($seq);

	next unless defined $genes{$seqid};
	print O1 "##sequence-region $seqid 1 $len\n";
	for($i = 0; $i< @{$genes{$seqid}}; $i++){
		$tmp = &deal_gff($genes{$seqid}[$i], $count, $seqid);
		print O1 $tmp;
		my $gene_len;
		($tmp, $gene_len) = &deal_seq($seq, $genes{$seqid}[$i], $count, $seqid);
		print O2 $tmp if $gene_len >= $min_len;
		$count++;
	}
}
$/ = "\n";
close I;
close O1;
close O2;


sub deal_gff {
	my($infos, $num, $id) = @_;
	my(@arry, $temps);
	@arry = split /\t/, $infos;
	$temps = $id."\tMetaGeneMark\tgene\t".$arry[0]."\t".$arry[1]."\t\.\t";
	$temps .= $arry[2]."\t\.\t";
	$temps .= sprintf("ID=GL%07d;Name=GL%07d;Note=", $num, $num);
	$temps .= $arry[-1]."\n";
	return $temps;
}

sub deal_seq {
	my($line, $infos, $num, $id) = @_;
	my($temps, @arry);
	@arry = split /\t/, $infos;
	$flag = 0;
	$id = sprintf(">%s%07d  [gene]  locus=%s:%d:%d:%s\[", $mark, $num, $id, $arry[0], $arry[1], $arry[2]);
	$id .= $arry[-1]."\]\n";
	$temps = substr($line, $arry[0] - 1, $arry[1] - $arry[0] + 1);
	$temps =~ tr/atcgn/ATCGN/;
	if($arry[2] eq '-'){
		$temps = reverse $temps;
		$temps =~ tr/ATCGatcg/TAGCtagc/;
	}
	$temps = $id.$temps."\n";
	return ($temps, $arry[1] - $arry[0] + 1);
}

