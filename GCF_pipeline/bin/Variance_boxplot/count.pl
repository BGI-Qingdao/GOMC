use strict;

my %info;
open IN, "2022_GCF_info_all.txt" or die $!;
while(<IN>){
    chomp;
    my @tmp  =split("\t");
    my $p = (split(";",$tmp[2]))[1];
    unless (grep /^$tmp[0]$/, @{$info{$p}}){
        push @{$info{$p}}, $tmp[0];
    }
}
close IN;

open OUT, ">2022_result.txt";
foreach my $p (sort keys %info){
    my $num = $#{$info{$p}}+1;
    print OUT "$p\t$num\n";
}
close OUT;