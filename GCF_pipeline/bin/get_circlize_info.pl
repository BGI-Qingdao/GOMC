use strict;
my $file_path = $ARGV[0];
my %phylums;
open TOTAL, "$file_path/temp/GCF_counts_per_phylim_total.txt" or die $!;
while(<TOTAL>){
    chomp;
    next if $. == 1;
    my @tmp = split("\t");
    $phylums{$tmp[0]}{total} = $tmp[1];
}
close TOTAL;

open UNIQ, "$file_path/temp/GCF_counts_per_phylim_uniq.txt" or die $!;
while(<UNIQ>){
    chomp;
    next if $. == 1;
    my @tmp = split("\t");
    $phylums{$tmp[0]}{uniq} = $tmp[1];
}
close UNIQ;

open OUT, ">$file_path/Circlize_Info.txt";
print OUT "Phylum\tTotal GCFs\tUniq GCFs\tShared GCFs\n";
foreach my $phylum (sort keys %phylums){
    my $total = $phylums{$phylum}{total};
    my $uniq = $phylums{$phylum}{uniq};
    my $shared = $total-$uniq;
    my $line = join("\t", $phylum,$total,$uniq,$shared);
    print OUT "$line\n";
}