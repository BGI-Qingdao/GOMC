use strict;

my ($in, $apd3, $out) = @ARGV;

open ADP, $apd3 or die $!;
my %known;
my %v_names;
while(<ADP>){
	chomp;
	if (/^>/){
		my $name = ($_ =~/>(.*)/)[0];
		my $num = int($./2+1);
		$v_names{$num} = $name;
	}else{
		my $num = $./2;
		$known{$_} = $v_names{$num};
	}
}
close ADP;

open ALL, $in or die $!;
open OUT, ">$out";
while(<ALL>){
	chomp;
	if ($. == 1){
		print OUT "$_\tAPD3_description\n";
	}else{
		my @tmp = split("\t");
		if ($tmp[-3]>0.5 && $tmp[-2]>0.5 && $tmp[-1]>0.5){
			my $apd = "-";
			$apd = "$known{$tmp[-4]}" if exists $known{$tmp[-4]};
			print OUT "$_\t$apd\n";
		}
	}
}
close ALL;
