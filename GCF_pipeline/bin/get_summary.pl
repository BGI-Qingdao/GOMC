use strict;

my ($bgc_count, $gcf_count, %count, %phylums);
open GCF_INFO, "$ARGV[0]/2_Bigslice_Output/GCF_info.txt" or die $!;
while(<GCF_INFO>){
    chomp;
    next if $. == 1;
    my @tmp = split("\t");
    $tmp[2] = (split(";",$tmp[2]))[1];
    unless (grep /^$tmp[1]$/, @{$count{bgc}}){
        push @{$count{bgc}}, $tmp[1];
    }
    unless (grep /^$tmp[1]$/, @{$phylums{$tmp[2]}{bgc}}){
        push @{$phylums{$tmp[2]}{bgc}}, $tmp[1];
    }
    unless (grep /^$tmp[0]$/, @{$count{gcf}}){
        push @{$count{gcf}}, $tmp[0];
    }
    unless (grep /^$tmp[0]$/, @{$phylums{$tmp[2]}{gcf}}){
        push @{$phylums{$tmp[2]}{gcf}}, $tmp[0];
    }
}
close GCF_INFO;

open RAREFACTION, "$ARGV[0]/4_Rarefaction/iNEXT_Info.txt" or die $!;
while(<RAREFACTION>){
    chomp;
    next if $. == 1;
    my @tmp = split("\t");
    $tmp[0] =~ s/\./-/g;
    $phylums{$tmp[0]}{rarefaction} = $tmp[2];
}

open OUT, ">$ARGV[0]/summary.txt";
print OUT "Phylum\tBGC number\tGCF number\tPotential GCFs\n";
foreach my $phylum (sort keys %phylums){
    my $bgc = $#{$phylums{$phylum}{bgc}}+1;
    my $gcf = $#{$phylums{$phylum}{gcf}}+1;
    my $rarefaction = $phylums{$phylum}{rarefaction};
    my $line = join ("\t", $phylum,$bgc,$gcf,$rarefaction);
    print OUT "$line\n";
}
$bgc_count = $#{$count{bgc}}+1;
$gcf_count = $#{$count{gcf}}+1;
print OUT "Total\t$bgc_count\t$gcf_count\n";
