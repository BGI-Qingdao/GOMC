use strict;

if (@ARGV != 4){
print "[ERROR] missing input file\n";
print "[EXAMPLE] perl $0 <GCF_phylum_order file> <GCF_info file> <GCF_dist file> <Output>\n";
exit;
}
my ($order,$gcf_info,$gcf_dist,$out) = @ARGV;

my @phylums;
open ORDER, $order or die $!;
while(<ORDER>){
chomp;
push @phylums, $_;
}
close ORDER;

my %fun;
open GCF_INFO, $gcf_info or die $!;
while(<GCF_INFO>){
chomp;
next unless /^\d/;
my @tmp = split("\t");
my $phylum = (split(";", $tmp[-1]))[-1];
$fun{$tmp[0]} = $phylum;
#push @phylums, $phylum unless grep /^$phylum$/, @phylums;
}
close GCF_INFO;

open GCF_DIST, $gcf_dist or die $!;
my %info;
while(<GCF_DIST>){
chomp;
next unless /^\d/;
my @tmp = split("\t");
my $value = get_value($tmp[1]);
push @{$info{$value}{$fun{$tmp[0]}}}, $tmp[0];
}
close GCF_DIST;

open OUT, ">$out";
print OUT "Cos_dist\tPhylum\tBGC_number\n";
my @values = qw/0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1/;
foreach my $value (@values){
unless ($info{$value}){
foreach my $f (@phylums){
@{$info{$value}{$f}}=();
}
}
foreach my $f (@phylums){
unless ($info{$value}{$f}){
@{$info{$value}{$f}}=();
}
my $number;

(@{$info{$value}{$f}})?($number=$#{$info{$value}{$f}}+1):($number=0);
my $line = join("\t", $value,$f,$number);
print OUT "$line\n";
}
}
close OUT;

sub get_value{
my ($in) = @_;
my $value;
$value = 0 if $in == 0;
$value = 0.1 if $in>0&&$in<=0.1;
$value = 0.2 if $in>0.1&&$in<=0.2;
$value = 0.3 if $in>0.2&&$in<=0.3;
$value = 0.4 if $in>0.3&&$in<=0.4;
$value = 0.5 if $in>0.4&&$in<=0.5;
$value = 0.6 if $in>0.5&&$in<=0.6;
$value = 0.7 if $in>0.6&&$in<=0.7;
$value = 0.8 if $in>0.7&&$in<=0.8;
$value = 0.9 if $in>0.8&&$in<=0.9;
$value = 1 if $in>0.9&&$in<=1;
return $value;
}

