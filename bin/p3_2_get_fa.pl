use strict;

my ($in, $out) = @ARGV;
#perl $0 <input gbk info file> <output core peptide fasta file>

open IN, $in;
open OUT, ">$out";
while(<IN>){
	chomp;
	next if $. == 1;
	my @tmp = split("\t");
	next if $tmp[-1] eq "-";
	my $mag = $tmp[0];
	my $type = $tmp[3];
	my $subtype = $tmp[4];
	my $gbk = ($tmp[1]=~/(.*)\.gbk/)[0];
	my @cores = split(",", $tmp[-1]);
	my $num = 0;
	foreach my $core(@cores){
		$num++;
		my $id = join("|", $mag,$gbk,$type,$subtype,"core_peptide_$num");
		print OUT ">$id\n$core\n";
	}
}
close IN;
close OUT;