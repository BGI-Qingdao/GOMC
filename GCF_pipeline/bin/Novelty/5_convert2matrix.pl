use strict;
if (@ARGV != 3){
print "[ERROR] missing input file\n";
print "[EXAMPLE] perl $0 <order file> <plot file> <converted file>\n";
exit;
}

my ($order_file,$plot_in,$convert_out) = @ARGV;

my @orders;
open ORDER, $order_file or die $!;
while(<ORDER>){
chomp;
push @orders, $_;
}
close ORDER;

my %info;
my @values;
open IN, $plot_in or die $!;
while(<IN>){
chomp;
next unless /^\d/;
my @tmp = split("\t");
$info{$tmp[0]}{$tmp[1]} = $tmp[2];
push @values, $tmp[0] unless grep /^$tmp[0]$/, @values;
}
close IN;

open OUT, ">$convert_out";
my $header = join("\t", "dist",@orders);
print OUT "$header\n";
foreach my $value (@values){
my $line = "$value";
foreach my $second (@orders){
$line.="\t$info{$value}{$second}";
}
print OUT "$line\n";
}
close OUT;
