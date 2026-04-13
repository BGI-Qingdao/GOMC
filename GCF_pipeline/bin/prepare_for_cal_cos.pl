use strict;
open BGC, "$ARGV[0]" or die $!;
open BGC_OUT, ">$ARGV[2]/BGC_features_input.tsv";
while(<BGC>){
    chomp;
    next if $. == 1;
    my @tmp = split("\t");
    shift @tmp;
    my $line = join("\t", @tmp);
    print BGC_OUT "$line\n";
}
close BGC;
close BGC_OUT;

my @gcf_models;
open GCF, "$ARGV[1]" or die $!;
open GCF_OUT, ">$ARGV[3]";
while(<GCF>){
    chomp;
    next if $. == 1;
    next if $. == 2;
    my @tmp = split("\t");
    my $gcf_model = shift @tmp;
    my $line = join("\t", @tmp);
    print GCF_OUT "$line\n";
    push @gcf_models, $gcf_model;
}
close GCF;
close GCF_OUT;

open GCF_HEADER, ">$ARGV[2]/GCF_header.list";
my $line = join("\t", @gcf_models);
print GCF_HEADER "$line\n";
close GCF_HEADER;