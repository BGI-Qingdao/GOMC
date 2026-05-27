use strict;
use Getopt::Long;

sub usage {
	die "
	Description: get GCF info

	Usage: perl $0 [Options]
	Options:
		-m	the relationship between mags and species, the format should be 'species\tmags\tmags_folder', which had been established in previous analysis and named 'target_MAGs.list'
		-b	the relationship between bgcs and mags, which had been established above and named 'BGC_info.txt'
		-g	the relationship between bgcs and gcfs, which had been established above and named 'BGC_membership.tsv'
		-o	the output path

	Example: perl $0 -m ./target_MAGs.list -b output_0804/BGC_info.txt -g output_0804/BGC_membership.tsv -o output/GCF_info.txt

	Contact: zhouchanghao\@genomics.cn
	Last Update: 2022/08/04
	";
}

my ($target_mags,$bgc_info,$bgc_member,$out);

&GetOptions(
	'm=s' => \$target_mags,
	'b=s' => \$bgc_info,
	'g=s' => \$bgc_member,
	'o=s' => \$out
);


unless (-e $target_mags || -e $bgc_info || -e $bgc_member || -e $out){
	die usage();
}

#################### Main #############################
my($sp_count,$bgc_count,$gcf_count);

my %mag_info;
my %species_count;
open MAG_INFO, $target_mags or die $!;
while (<MAG_INFO>){
	chomp;
	my @tmp = split("\t");
	$mag_info{$tmp[1]} = $tmp[0];
}
close MAG_INFO;

my %bgc_info;
open BGC_INFO, $bgc_info or die $!;
while (<BGC_INFO>){
	chomp;
	next if $.==1;
	my @tmp = split("\t");
	$bgc_info{$tmp[0]} = $mag_info{$tmp[1]} if exists $mag_info{$tmp[1]};
}
close BGC_INFO;

my %final_info;
open MEMBERSHIP, $bgc_member or die $!;
<MEMBERSHIP>;
while(<MEMBERSHIP>){
	chomp;
	next if $.==1;
	my @tmp = split("\t");
#	next if $tmp[1]==0;
	$final_info{$tmp[1]}{$tmp[0]} = $bgc_info{$tmp[0]};
}
close MEMBERSHIP;

my %stat;
open OUT, ">$out";
print OUT "GCF\tBGC\tPhylum\n";
foreach my $gcf (sort {$a<=>$b} keys %final_info){
	my $per_gcf_bgc_count=0;
	foreach my $bgc (sort {$a<=>$b} keys $final_info{$gcf}){
		my $line = join("\t", $gcf,$bgc,$final_info{$gcf}{$bgc});
		print OUT "$line\n";
		$bgc_count++;
		$per_gcf_bgc_count++;
	}
	$gcf_count++;
	$stat{$per_gcf_bgc_count}++;
}
close OUT;
open COUNT, ">$out.stat.txt";
print COUNT "Total BGCs: $bgc_count\n";
print COUNT "Total GCFs: $gcf_count\n";
print COUNT "And the details were shown below:\n\n";
print COUNT "Numbers_of_BGC_per_GCF\tNumbers_of_GCF\n";
foreach my $count (sort {$a<=>$b} keys %stat){
	print COUNT "$count\t$stat{$count}\n";
}
close COUNT;
