use strict;
open IN, "$ARGV[0]/temp/link.txt";
open OUT, ">$ARGV[0]/temp/link_simplify.txt";
my $gcf = 0;
while (<IN>){
    chomp;
    my @tmp = split("\t");
    my $gcf_in = join(":", $tmp[0],$tmp[1]);
    if ($gcf_in ne $gcf){
        print OUT "$_\n";
        $gcf = $gcf_in;
    }
}
close IN;
close OUT;