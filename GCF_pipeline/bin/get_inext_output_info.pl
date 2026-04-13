use strict;
my %info;
# foreach my $file (glob "$ARGV[0]/iNEXT/iNEXT_temp/*.txt"){
foreach my $file (glob "$ARGV[0]/temp/*.txt"){
    my ($phylum,$observed,$total,$potential);
    open FILE, "$file" or die $!;
    while (<FILE>){
        chomp;
        my @tmp = split("\t");
        if ($.==1){
            $phylum = ($tmp[0] =~ /(.*)\.t/)[0];
        }
        $observed = $tmp[3] if $tmp[1] eq "observed";
        $total = $tmp[3] if ($tmp[3] >= $total && $tmp[1] eq "extrapolated");
    }
    close FILE;
    $potential = int($total-$observed);
    @{$info{$phylum}} = ($observed,$potential);    
}
open OUT, ">$ARGV[0]/iNEXT_Info.txt";
print OUT "Phylum\tObserved_GCFs\tPotential_GCFs\n";
foreach my $phylum (sort keys %info){
    my $line = join("\t", $phylum,@{$info{$phylum}});
    print OUT "$line\n";
}
close OUT;