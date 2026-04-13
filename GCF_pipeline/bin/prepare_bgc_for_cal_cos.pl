use strict;
open BGC, "$ARGV[0]" or die $!;
open BGC_OUT, ">$ARGV[1]";

my @gcf_models;
while(<BGC>){
    chomp;
    # next if $. == 1;
    if ($. == 1){
        @gcf_models = split("\t", $_);
        shift @gcf_models;
    }else{
        my @tmp = split("\t");
        shift @tmp;
        my $line = join("\t", @tmp);
        print BGC_OUT "$line\n";
    }
}
close BGC;
close BGC_OUT;

