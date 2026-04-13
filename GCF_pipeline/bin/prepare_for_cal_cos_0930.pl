use strict;
open BGC, "$ARGV[0]" or die $!;
open BGC_OUT, ">$ARGV[2]/BGC_features_input.tsv";

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

open GCF, "$ARGV[1]" or die $!;
open GCF_OUT, ">$ARGV[3]/GCF_models_input.tsv";
while(<GCF>){
    chomp;
    next if $. == 1;
    my @tmp = split("\t");
    shift @tmp;
    my $line = join("\t", @tmp);
    print GCF_OUT "$line\n";
}
close GCF;
close GCF_OUT;

open GCF_HEADER, ">$ARGV[2]/GCF_header.list";
my $line = join("\t", @gcf_models);
print GCF_HEADER "$line\n";
close GCF_HEADER;