use strict;
my %genome_info;
# open IN, "/dellfsqd2/ST_OCEAN/USER/zhouchanghao/metagenome/get_gcf/20220905/Output/Variance/2022_total_bin_metadata" or die $!;
open IN, $ARGV[0] or die $!;
while(<IN>){
    chomp;
    # my @tmp = split("\t",$_,15);
    # $genome_info{$tmp[0]} = $tmp[13];
    my @tmp = split("\t");
    $genome_info{$tmp[0]} = $tmp[1];
}
close IN;

my %bgc_info;
# open BGC, "/dellfsqd2/ST_OCEAN/USER/zhouchanghao/metagenome/get_gcf/20220905/Output/2_Bigslice_Output/BGC_info.txt" or die $!;
open BGC, "$ARGV[1]/2_Bigslice_Output/BGC_info.txt" or die $!;
while (<BGC>){
    chomp;
    next if $. == 1;
    my @tmp = split("\t");
    $bgc_info{$tmp[0]} = $genome_info{$tmp[1]};
}
close BGC;

# open GCF, "/dellfsqd2/ST_OCEAN/USER/zhouchanghao/metagenome/get_gcf/20220905/Output/2_Bigslice_Output/GCF_info.txt" or die $!;
open GCF, "$ARGV[1]/2_Bigslice_Output/GCF_info.txt" or die $!;
open OUT, ">./1_GCF_info_all.txt";
while(<GCF>){
    chomp;
    next if $. == 1;
    my @tmp = split("\t");
    $tmp[2] = $bgc_info{$tmp[1]};
    my $line = join("\t", @tmp);
    print OUT "$line\n";
}
close GCF;
close OUT;
