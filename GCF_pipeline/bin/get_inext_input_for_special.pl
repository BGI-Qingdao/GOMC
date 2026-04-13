use strict;
open BGC, "$ARGV[1]/BGC_info.txt" or die $!;
my %bgcs;
while(<BGC>){
	chomp;
	next if $. == 1;
	my @tmp = split("\t");
	push @{$bgcs{$tmp[1]}}, $tmp[0];
}
close BGC;

my %phylums;
open MAG, "target_mags.txt" or die $!;
while (<MAG>){
	chomp;
	my @tmp = split("\t");
	$tmp[0] = (split(";",$tmp[0]))[1];
	if (exists $bgcs{$tmp[1]}){
		push @{$phylums{$tmp[0]}}, @{$bgcs{$tmp[1]}};
	}
}
close MAG;

open GCF, "$ARGV[1]/GCF_info.txt" or die $!;
my %gcfs;
while(<GCF>){
	chomp;
	next if $. == 1;
	my @tmp = split("\t");
	push @{$gcfs{$tmp[0]}}, $tmp[1];
}
close GCF;

`mkdir -p "$ARGV[2]/iNEXT_Input"` unless -e "$ARGV[2]/iNEXT_Input";
`mkdir -p "$ARGV[2]/temp"` unless -e "$ARGV[2]/temp";
my @r_list;
foreach my $phylum (sort keys %phylums){
	push @r_list, "$phylum=NULL";
	open OUT, ">$ARGV[2]/iNEXT_Input/$phylum.txt";
	my $line = join("\t", "\t",@{$phylums{$phylum}});
	print OUT "$line\n";
	foreach my $gcf (sort keys %gcfs){
		print OUT "$gcf";
		foreach my $bgc (@{$phylums{$phylum}}){
			if (grep /^$bgc$/, @{$gcfs{$gcf}}){
				print OUT "\t1";
			}else{
				print OUT "\t0";
			}
		}
		print OUT "\n";
	}
	close OUT;
}