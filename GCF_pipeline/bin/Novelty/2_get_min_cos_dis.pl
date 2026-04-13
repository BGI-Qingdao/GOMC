use strict;
use List::Util qw/min/;

if (@ARGV != 5){
print "[ERROR] Missing parameters. Please input 1 GCF info file, 1 BGC dist file, and set the BGC min dist result file, GCF dist result file and threshold\n";
print "[Example] perl $0 <GCF_info.txt> <All_dist_result.txt> <BGC_min_dist_result.txt> <GCF_dist_result.txt> <0.2>\n";
exit;
}

my %gcf;
open GCF, "$ARGV[0]" or die $!;
while(<GCF>){
chomp;
next if $. == 1;
my @tmp = split("\t");
$gcf{$tmp[1]} = $tmp[0];
}
close GCF;

open IN, "$ARGV[1]" or die $!;
open BGC_OUT, ">$ARGV[2]";
my %info;
while(<IN>){
chomp;
my @tmp = split("\t");
my $min_value = min @tmp;
push @{$info{$gcf{$.}}}, $min_value;
my $line = join("\t", $.,$min_value);
print BGC_OUT "$line\n";
}
close IN;
close BGC_OUT;

open GCF_OUT, ">$ARGV[3]";
print GCF_OUT "GCF\tDistance to Big-Fam\tNovetly\n";
foreach my $g (sort {$a<=>$b} keys %info){
my $sum = 0;
my $num = $#{$info{$g}}+1;
foreach my $value (@{$info{$g}}){
$sum += $value;
}
my $average = $sum/$num;
my $threshold = $ARGV[4];
my $tag;
($average<=$threshold)?($tag = "No"):($tag = "Yes");
my $line = join("\t", $g,$average,$tag);
print GCF_OUT "$line\n";
}
close GCF_OUT;


