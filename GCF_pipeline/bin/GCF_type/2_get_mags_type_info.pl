use strict;

if (@ARGV != 2){
print "[ERROR] missing input file\n";
print "[EXAMPLE] perl $0 <BGC merge info file> <output file>\n";
exit;
}

my ($in,$out) = @ARGV;

my %info;
my @types;
open IN, $in or die $!;
while(<IN>){
chomp;
next unless /^\d/;
my @tmp = split("\t");
push @{$info{$tmp[2]}{$tmp[3]}{gcf}}, $tmp[0] unless grep /^$tmp[0]$/, @{$info{$tmp[2]}{$tmp[3]}{gcf}};
push @{$info{$tmp[2]}{$tmp[3]}{bgc}}, $tmp[1] unless grep /^$tmp[1]$/, @{$info{$tmp[2]}{$tmp[3]}{bgc}};
push @types, $tmp[3] unless grep /^$tmp[3]$/, @types;
}
close IN;

open OUT, ">$out";
my $header = join("\t", "Phylum","BGC_number",@types,"GCF_number",@types);
print OUT "$header\n";
foreach my $phy (sort keys %info){
my $b_c = 0;
my $g_c = 0;
my @b_numbers;
my @g_numbers;
foreach my $type (@types){
my $b_count;
my $g_count;
($info{$phy}{$type}{bgc})?($b_count = $#{$info{$phy}{$type}{bgc}}+1):($b_count = 0);
($info{$phy}{$type}{gcf})?($g_count = $#{$info{$phy}{$type}{gcf}}+1):($g_count = 0);
$b_c += $b_count;
$g_c += $g_count;
push @b_numbers, $b_count;
push @g_numbers, $g_count;
}
my $line = join("\t", $phy,$b_c,@b_numbers,$g_c,@g_numbers);
print OUT "$line\n";
}
close OUT;
