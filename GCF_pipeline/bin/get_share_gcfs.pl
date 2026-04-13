use strict;

my @phylums;
my %gcf;
open IN, "$ARGV[0]" or die $!;
while (<IN>){
	chomp;
	next if $. == 1;
	my @tmp = split("\t");
	my $phylum = (split(";", $tmp[2]))[1];
	unless (grep /^$phylum$/, @{$gcf{$tmp[0]}}){
		push @{$gcf{$tmp[0]}}, $phylum;
	}
	unless (grep /^$phylum$/, @phylums){
		push @phylums, $phylum;
	}
}
close IN;

`mkdir -p $ARGV[1]/temp` unless -e "$ARGV[1]/temp";
my (%count_uniq, %count_common);
open LINK, ">$ARGV[1]/temp/link.txt";
# foreach my $cluster (sort keys %gcf){
foreach my $cluster (sort {$a<=>$b} keys %gcf){
	if ($#{$gcf{$cluster}} == 0){
		$count_uniq{${$gcf{$cluster}}[0]}++;
	}else{
		foreach my $phylum (@{$gcf{$cluster}}){
			$count_common{$phylum}++;
		}
		foreach my $query (0..$#{$gcf{$cluster}}){
			foreach my $sub ($query+1..$#{$gcf{$cluster}}){
				my $line = join("\t", $cluster,${$gcf{$cluster}}[$query],${$gcf{$cluster}}[$sub]);
				print LINK "$line\n";
			}
		}
	}
}
close LINK;

open UNIQ, ">$ARGV[1]/temp/GCF_counts_per_phylim_uniq.txt";
print UNIQ "Phylum\tNumber(Uniq)\n";
open TOTAL, ">$ARGV[1]/temp/GCF_counts_per_phylim_total.txt";
print TOTAL "Phylum\tNumber(Total)\n";
foreach my $phylum (@phylums){
	my ($line_uniq,my $line_total);
	#统计每个phylum中uniq的GCF数目
	if (exists $count_uniq{$phylum}){
		$line_uniq = join("\t", $phylum,$count_uniq{$phylum});
	}else{
		$line_uniq = join("\t", $phylum,0);
	}
	#统计每个phylum中common的GCF数目
	if (exists $count_common{$phylum}){
		if (exists $count_uniq{$phylum}){
			$line_total = join("\t", $phylum,$count_uniq{$phylum}+$count_common{$phylum});
		}else{
			$line_total = join("\t", $phylum,$count_common{$phylum});
		}
	}else{
		if (exists $count_uniq{$phylum}){
			$line_total = join("\t", $phylum,$count_uniq{$phylum});
		}else{
			$line_total = join("\t", $phylum,0);
		}
	}
	print UNIQ "$line_uniq\n";
	print TOTAL "$line_total\n";
}
close UNIQ;
close TOTAL;