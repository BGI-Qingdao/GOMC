use strict;

#perl $0 <mag_info.txt> <antismash_result_path/> <out_path/>

my ($mag_info,$antismash_result_path,$out_path) = @ARGV;

`ln -s $antismash_result_path $out_path/dataset`;

open TSV, ">$out_path/datasets.tsv";
print TSV "# Dataset name\tPath to folder\tPath to taxonomy\tDescription\n";
print TSV "dataset\tdataset/\t\tnothing\n";
close TSV;

open IN, $mag_info;
open OUT, ">$out_path/target_MAGs.list";
while(<IN>){
	chomp;
	my @tmp = split("\t");
	print OUT "$tmp[1]\t$tmp[0]\t$out_path/dataset/$tmp[0]\n";
}
close IN;
close OUT;