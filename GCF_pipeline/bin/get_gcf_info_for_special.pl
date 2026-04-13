use strict;
my @phylums;
open IN, $ARGV[1] or die $!; 
while(<IN>){
    chomp;
    # push @phylums, $_;
    my $phy = (split("\t"))[0];
    push @phylums, $phy;
}
close IN;

open GCF, "$ARGV[0]/2_Bigslice_Output/GCF_info.txt" or die $!;
open OUT, ">GCF_Info.txt";
while(<GCF>){
    chomp;
    my @tmp = split("\t");
    if ($. == 1){
        print OUT "$_\n";
    }else{
        my $phylum = (split(";", $tmp[2]))[1];
        if (grep /^$phylum$/, @phylums){
            print OUT "$_\n";
        }
    }
}
close GCF;
close OUT;