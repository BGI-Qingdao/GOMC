use strict;
use Getopt::Long;

# sub usage {
# 	die "
# 	Description: prepare files for bigslice

# 	Usage: perl $0 [Options]
# 	Options:
# 		-p	prepare folder path
# 		-o	output folder path

# 	Example: perl $0 -p prepare -o bigslice_input

# 	Contact: zhouchanghao\@genomics.cn
# 	Last Update: 2022/07/26
# 	";
# }

# my ( $output,$dataset_name,$prepare_dir);
my ( $phylum_list,$mags,$antismash_o,$output);

&GetOptions(
	'o=s' => \ $output,
	'm=s' => \ $mags,
	'a=s' => \ $antismash_o,
	'p=s' => \$phylum_list,
);

################# Main ###################
# #创建目录
unless (-e " $output/dataset"){
	`mkdir -p $output/dataset`;
}
# unless (-e " $output/taxonomy"){
# 	`mkdir -p  $output/taxonomy`;
# }

# my(%list, $antismash_path,)

my $partial;
if ($phylum_list eq "ALL"){
	$partial = "no";
}else{
	$partial = "yes";
}

#读取门列表
my %list;
if ($partial eq "yes"){
	open LIST, $phylum_list or die $!;
	while(<LIST>){
		chomp;
		$list{$_}=$_;
	}
	close LIST;
}

open MAGS, $mags or die $!;
open OUT, ">$output/target_MAGs.list";
# open TAXONOMY, "> $output/taxonomy/taxonomy_dataset.tsv";
# print TAXONOMY "# Genome folder\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\tOrganism\n";
while (<MAGS>){
	chomp;
	my @tmp = split("\t");
	my ($file, $dataset_line, $dataset_dir);
	if ($partial eq "yes"){
		if (exists $list{$tmp[1]}){
			$file = "$antismash_o/$tmp[0]";
			$dataset_line = join("\t", $tmp[1],$tmp[0],$file);
			print OUT "$dataset_line\n";

			#创建各基因组文件夹并拷贝gbk到该文件夹下
			my $dataset_dir = " $output/dataset/$tmp[0]";
			`mkdir $dataset_dir`;
			`cp $file/*region*gbk $dataset_dir/`;
		}
	}else{
		$file = "$antismash_o/$tmp[0]";
		$dataset_line = join("\t", $tmp[1],$tmp[0],$file);
		print OUT "$dataset_line\n";

		#创建各基因组文件夹并拷贝gbk到该文件夹下
		my $dataset_dir = " $output/dataset/$tmp[0]";
		`mkdir $dataset_dir`;
		`cp $file/*region*gbk $dataset_dir/`;
	}
}
close MAGS;
close OUT;

open TSV, "> $output/datasets.tsv" or die $!;
print TSV "# Dataset name\tPath to folder\tPath to taxonomy\tDescription\n";
print TSV "dataset\tdataset/\t\tnothing\n";
close TSV;